/*  thread_pool.c -- Implementation for the thread pool.

    Copyright (c) 2013-2016 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

    Minified by: Ziwei Xue

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
 * This file implements a thread pool for multi-threading applications.
 * It consists of two distinct interfaces: thread pools an thread job queues.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling except to pick tasks from
 * queues that have room to store the result.
 *
 * Upon completion, the return value from the function pointer is
 * added to back to the queue if the result is required.  We may have
 * multiple queues in use for the one pool.
 *
 * To see example usage, please look at the #ifdef TEST_MAIN code in
 * thread_pool.c.
 */

#include <stdlib.h>
#include <inttypes.h>
#include <signal.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>
#include "thread_pool.h"

/** Minimum stack size for threads.  Required for some rANS codecs
 * that use over 2Mbytes of stack for encoder / decoder state
**/
#define MIN_THREAD_STACK (3 * 1024 * 1024)


tpool_t* tpool_init(int n) 
{
    int t_idx = 0;
    size_t stack_size = 0;
    pthread_attr_t pattr;
    bool pattr_init_done;
    tpool_t* p = malloc(sizeof(*p));
    if (!p) {return NULL;}

    p->tsize = n;
    p->njobs = 0;
    p->nwaiting = 0;
    p->shutdown = false;
    p->q_head = NULL;
    p->t_stack = NULL;
    p->n_count = 0;
    p->n_running = 0;
    p->t = malloc(n * sizeof(tpool_worker_t));
    if (!p->t) {free(p); return NULL;}
    p->t_stack = malloc(n * sizeof(int));
    if (!p->t_stack) {free(p);free(p->t); return NULL;}

    p->t_stack_top = -1;

    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&p->tpool_mu, &attr);
    pthread_mutexattr_destroy(&attr);  

    pthread_mutex_lock(&p->tpool_mu);

    if (pthread_attr_init(&pattr) < 0) goto cleanup;
    pattr_init_done = true;
    if (pthread_attr_getstacksize(&pattr, &stack_size) < 0) goto cleanup;
    if (stack_size < MIN_THREAD_STACK) {
        if (pthread_attr_setstacksize(&pattr, MIN_THREAD_STACK) < 0) goto cleanup;
    }
    pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);
    /* initialization of worker threads */
    for (t_idx = 0; t_idx < n; t_idx++) {
        tpool_worker_t *w = &p->t[t_idx];
        p->t_stack[t_idx] = 0;
        w->p = p;
        w->idx = t_idx;
        pthread_cond_init(&w->pending_c, NULL);
        if (0 != pthread_create(&w->tid, &pattr, tpool_worker, w)) goto cleanup;
    }

    pthread_mutex_unlock(&p->tpool_mu);
    pthread_attr_destroy(&pattr);
    return p;
cleanup:
    p->shutdown = true;
    for (int j = 0; j < t_idx; j++) {
        pthread_join(p->t[j].tid, NULL);
        pthread_cond_destroy(&p->t[j].pending_c);
    }
    pthread_mutex_destroy(&p->tpool_mu);
    if (pattr_init_done) pthread_attr_destroy(&pattr);
    free(p->t_stack);
    free(p->t);
    free(p);
    return NULL;
}

void tpool_destroy(tpool_t *p)
{
    pthread_mutex_lock(&p->tpool_mu);
    p->shutdown = true;
    /* broadcast shutdown message to all worker threads */
    for (int i = 0; i < p->tsize; i++) {
        pthread_cond_signal(&p->t[i].pending_c);
    }
    pthread_mutex_unlock(&p->tpool_mu);
    for (int i = 0; i < p->tsize; i++) {
        pthread_join(p->t[i].tid, NULL);
    }
    pthread_mutex_destroy(&p->tpool_mu);
    
    for (int i = 0; i < p->tsize; i++) {
        pthread_cond_destroy(&p->t[i].pending_c);
    }
    if (p->t_stack) free(p->t_stack);
    free(p->t);
    free(p);
}

tpool_process_t *tpool_process_init(tpool_t *p, int qsize, bool in_only)
{
    tpool_process_t *q = malloc(sizeof(*q));
    if (!q) return NULL;

    pthread_cond_init(&q->output_avail_c,   NULL);
    pthread_cond_init(&q->input_not_full_c, NULL);
    pthread_cond_init(&q->input_empty_c,    NULL);
    pthread_cond_init(&q->none_processing_c,NULL);

    q->p = p;
    q->job_head = NULL;
    q->job_tail = NULL; 
    q->result_head = NULL;
    q->result_tail = NULL;
    q->next_serial = 0;
    q->curr_serial = 0;
    q->no_more_input = 0;
    q->n_job     = 0;
    q->n_result    = 0;
    q->n_processing= 0;
    q->qsize       = qsize;
    q->in_only     = in_only;
    q->shutdown    = 0;
    q->wake_dispatch = 0;
    q->ref_count   = 1;

    q->next        = NULL;
    q->prev        = NULL;
    tpool_process_attach(p,q);
    return q;    
}

void tpool_process_destroy(tpool_process_t *q)
{
    if (!q) return;
    pthread_mutex_lock(&q->p->tpool_mu);
    q->no_more_input = true;
    pthread_mutex_unlock(&q->p->tpool_mu);

    tpool_process_reset(q,0);

    pthread_mutex_lock(&q->p->tpool_mu);
    
    tpool_process_detach(q->p, q);
    
    tpool_process_shutdown(q);
    
    if (--q->ref_count > 0) {
        pthread_mutex_unlock(&q->p->tpool_mu);
        return;
    }

    pthread_cond_destroy(&q->output_avail_c);
    pthread_cond_destroy(&q->input_not_full_c);
    pthread_cond_destroy(&q->input_empty_c);
    pthread_cond_destroy(&q->none_processing_c);
    pthread_mutex_unlock(&q->p->tpool_mu);
    free(q);
    return;
}

void tpool_process_attach(tpool_t *p, tpool_process_t *q)
{
    pthread_mutex_lock(&p->tpool_mu);
    /* insert the process into the process linked-list */
    if (p->q_head) {
        q->next = p->q_head;
        q->prev = p->q_head->prev;
        p->q_head->prev->next = q;
        p->q_head->prev = q;
        p->q_head = q;
    } else {
        p->q_head = q;
        q->next = q;
        q->prev = q;
    }
    pthread_mutex_unlock(&p->tpool_mu);
}

void tpool_process_detach(tpool_t *p, tpool_process_t *q) 
{
    pthread_mutex_lock(&p->tpool_mu);
    if (!p->q_head || !q->next || !q->prev ) return;
    tpool_process_t *curr = p->q_head, *first = curr;
    do {
        /* find and detach */
        if (curr == q) {
            curr->prev->next = curr->next;
            curr->next->prev = curr->prev;
            p->q_head = curr->next;
            q->next = q->prev = NULL;

            /* one left */
            if (p->q_head == q) {
                p->q_head = NULL;
            }
            break;
        }
        curr = curr->next;
    } while (curr != first);
    pthread_mutex_unlock(&p->tpool_mu);
}

int tpool_process_flush(tpool_process_t *q)
{
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    for (int i = 0; i < p->tsize; i++) {
        if (p->t_stack[i])
            pthread_cond_signal(&p->t[i].pending_c);
    }
    if (q->qsize < q->n_result + q->n_job + q->n_processing)
        q->qsize = q->n_result + q->n_job + q->n_processing;
    if (q->shutdown) {
        while (q->n_processing)
            pthread_cond_wait(&q->none_processing_c, &p->tpool_mu);
    }
    /* Wait for n_job and n_processing to hit zero */
    while (!q->shutdown && (q->n_job || q->n_processing)) {
        struct timeval now;
        struct timespec timeout;

        while (q->n_job && !q->shutdown) {
            gettimeofday(&now, NULL);
            timeout.tv_sec = now.tv_sec + 1;
            timeout.tv_nsec = now.tv_usec * 1000;
            pthread_cond_timedwait(&q->input_empty_c, &p->tpool_mu, &timeout);
        }

        while (q->n_processing > 0) {
            gettimeofday(&now, NULL);
            timeout.tv_sec = now.tv_sec + 1;
            timeout.tv_nsec = now.tv_usec * 1000;
            pthread_cond_wait(&q->none_processing_c, &p->tpool_mu); 
        }
        if (q->shutdown) break;
    }
    pthread_mutex_unlock(&p->tpool_mu);
    return 0;
}

static void *tpool_worker(void *arg)
{
    tpool_worker_t *w = (tpool_worker_t *)arg;
    tpool_t *p = w->p;
    tpool_job_t *j;

    pthread_mutex_lock(&p->tpool_mu);
    while (!p->shutdown) {
        assert(p->q_head == 0 || (p->q_head->prev && p->q_head->next));
        bool work2do = 0;
        tpool_process_t *first = p->q_head, *q = first;
        do {
            if (q && q->job_head && q->qsize - q->n_result > p->tsize - p->nwaiting && !q->shutdown) {
                work2do = true;
                break;
            }
            if (q) q = q->next;
        } while (q && q != first);
        if (!work2do) {
            /* no queues available, wait */
            p->nwaiting++;
            /* push this thread to the top of the waiting stack */
            if (p->t_stack_top == -1 || p->t_stack_top > w->idx)
                p->t_stack_top = w->idx;
            p->t_stack[w->idx] = 1;
            pthread_cond_wait(&w->pending_c, &p->tpool_mu);


            /* Finish and find new stack top */
            p->t_stack[w->idx] = 0;
            int i;
            p->t_stack_top = -1;
            for (int i = 0; i < p->tsize; ++i) {
                if (p->t_stack[i]) {
                    p->t_stack_top = i;
                    break;
                }
            }
            p->nwaiting--;
            continue;
        }
        q->ref_count++;
        while (q->job_head && q->qsize - q->n_result > q->n_processing) {
            if (p->shutdown) goto shutdown;
            if (q->shutdown) break;
            j = q->job_head;
            assert(j->p == p);
            if (!(q->job_head = j->next)) q->job_tail = NULL;

            q->n_processing++;
            if (q->n_job-- >= q->qsize) pthread_cond_broadcast(&q->input_not_full_c);

            if (q->n_job == 0) pthread_cond_broadcast(&q->input_empty_c);
            p->njobs--;
            pthread_mutex_unlock(&p->tpool_mu);
            if (tpool_add_result(j, j->func(j->arg)) < 0) goto err;
            free(j);
            pthread_mutex_lock(&p->tpool_mu);
        }
        if (--q->ref_count == 0) {
            tpool_process_destroy(q);
        } else {
            if (p->q_head) p->q_head = p->q_head->next;
        }
    }
shutdown:
    pthread_mutex_unlock(&p->tpool_mu);
    return NULL;
err:
    pthread_mutex_lock(&p->tpool_mu);
    tpool_process_t* first = p->q_head, *q = first;
    if (q) {
            do {
                tpool_process_shutdown(q);
                q->shutdown = 2; // signify error.
                q = q->next;
            } while (q != first);
    }   
    return NULL;
}

int tpool_process_reset(tpool_process_t *q, bool _free)
{
    tpool_job_t *j, *jn, *jhead;
    tpool_result_t *r, *rn, *rhead;

    pthread_mutex_lock(&q->p->tpool_mu);
    q->next_serial = INT_MAX;

    /* clear all inputs */
    jhead = q->job_head;
    q->job_head = q->job_tail = NULL;
    q->n_job = 0;

    /* clear all results */
    rhead = q->result_head;
    q->result_head = q->result_tail = NULL;
    q->n_result = 0;
    /* unlock here as the process linked-list is detached */
    pthread_mutex_unlock(&q->p->tpool_mu);
    /* release memory using the cleanup function */
    for (j = jhead; j; j = jn) {
        jn = j->next;
        if (j->job_cleanup) j->job_cleanup(j->arg);
        free(j);
    }

    for (r = rhead; r; r = rn) {
        rn = r->next;
        if (r->result_cleanup) r->result_cleanup(r->data);
        if (r && _free) {
            free(r->data);
            r->data = NULL;
            free(r);
        }
    }
    /* Wait for any jobs being processed to complete */
    if (tpool_process_flush(q) != 0) return -1;
    /* clean any new reuslt */
    pthread_mutex_lock(&q->p->tpool_mu);
    rhead = q->result_head;
    q->result_head = q->result_tail = NULL;
    q->n_result = 0;
    q->next_serial = q->curr_serial = 0;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&q->p->tpool_mu);
    /* again we clean the unwanted result */
    for (r = rhead; r; r = rn) {
        rn = r->next;
        if (r->result_cleanup) r->result_cleanup(r->data);
        if (r && _free) {
            free(r->data);
            r->data = NULL;
            free(r);
        }
    }
    return 0;

}

int tpool_dispatch(tpool_t *p, 
                   tpool_process_t *q, 
                   void *(*func)(void *arg),
                   void *arg,
                   void (*job_cleanup)(void *arg),
                   void (*result_cleanup)(void *data),
                   bool nonblock)
{
    tpool_job_t *j;
    pthread_mutex_lock(&p->tpool_mu);
    if ((q->no_more_input || q->n_job >= q->qsize  && nonblock)) {
        pthread_mutex_unlock(&p->tpool_mu);;
        return -1;
    }
    if (!(j = malloc(sizeof(*j)))) {
        pthread_mutex_unlock(&p->tpool_mu);
        return -1;
    }
    j->func = func;
    j->arg = arg;
    j->job_cleanup = job_cleanup;
    j->result_cleanup = result_cleanup;
    j->next = NULL;
    j->p = p;
    j->q = q;
    j->serial = q->curr_serial++;
    if (!nonblock) {
        while (q->no_more_input || q->n_job >= q->qsize && !q->shutdown && !q->wake_dispatch) {
            pthread_cond_wait(&q->input_not_full_c, &q->p->tpool_mu);
        }
        if (q->no_more_input || q->shutdown) {
            free(j);
            pthread_mutex_unlock(&p->tpool_mu);
            return -1;
        }
        if (q->wake_dispatch) {
            q->wake_dispatch = 0;
        }
    }
    p->njobs++;
    q->n_job++;
    if (q->job_tail) {
        q->job_tail->next = j;
        q->job_tail = j;
    } else {
        q->job_head = q->job_tail = j;
    }

    if (!q->shutdown) wake_next_worker(q,true);
    pthread_mutex_unlock(&p->tpool_mu);
    return 0;
}

void tpool_wake_dispatch(tpool_process_t *q)
{
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    q->wake_dispatch = 1;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&p->tpool_mu);
}

static void wake_next_worker(tpool_process_t *q, bool locked)
{
    if (!q) return;
    tpool_t *p = q->p;
    if (!locked) {
        pthread_mutex_lock(&p->tpool_mu);
    }
    assert(q->prev && q->next); /* the process has been attached */
    p->q_head = q;
    assert(p->njobs == q->n_job);

    int running = p->tsize - p->nwaiting;
    int sig = p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting && (q->n_processing < q->qsize - q->n_result);

    if (sig) {
        pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
    }
    if (!locked) {
        pthread_mutex_unlock(&p->tpool_mu);
    }
}

static int tpool_add_result(tpool_job_t *j, void *data)
{
    tpool_t *p = j->p;
    tpool_process_t *q = j->q;
    tpool_result_t *r;
    pthread_mutex_lock(&p->tpool_mu);
    
    if (--q->n_processing == 0) {
        int r;
        if ((r = pthread_cond_signal(&q->none_processing_c)))
            printf("output_avail_c signal failed with code %d\n", r);
    }
    if (q->in_only) {
        pthread_mutex_unlock(&p->tpool_mu);
        return 0;
    }
    
    r = malloc(sizeof(*r));
    if (!r) {
        pthread_mutex_unlock(&p->tpool_mu);
        tpool_process_shutdown(q);
        return -1;
    }

    r->next = NULL;
    r->data = data;
    r->result_cleanup = j->result_cleanup;
    r->serial = j->serial;

    q->n_result++;
    if (q->result_tail) {
        q->result_tail->next = r;
        q->result_tail = r;
    } else {
        q->result_head = q->result_tail = r;
    }
    
    if (r->serial == q->next_serial) {
        int r;
        
        if ((r = pthread_cond_broadcast(&q->output_avail_c)))
            printf("output_avail_c broadcast failed with code %d\n", r);
    }
    pthread_mutex_unlock(&p->tpool_mu);
    return 0;
}

tpool_result_t *tpool_next_result(tpool_process_t *q)
{
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    tpool_result_t *r, *last;
    if (q->shutdown) return NULL;

    /* searching for ordered output by serial number */
    
    for (last = NULL, r = q->result_head; r; last = r, r = r->next) {
        if (r->serial == q->next_serial)
            break;
    }

    if (r) {
        /* remove r from out linked-list */
        if (q->result_head == r) q->result_head = r->next;
        else last->next = r->next;
        
        if (q->result_tail == r) q->result_tail = last;
        if (!q->result_head) q->result_tail = NULL;

        q->next_serial++;
        q->n_result--;

        if (q->qsize && q->n_result < q->qsize) {
            if (q->n_job < q->qsize)
                pthread_cond_signal(&q->input_not_full_c);
            if (!q->shutdown)
                wake_next_worker(q,1);
        }
    }
    pthread_mutex_unlock(&p->tpool_mu);
    return r;
}

tpool_result_t *tpool_next_result_wait(tpool_process_t *q)
{
    tpool_result_t *r;
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    /* non-timed wait 
    if (q->n_job || q->n_processing || q->n_result)
        pthread_cond_wait(&q->output_avail_c, &p->tpool_mu);
    */
    while (!(r = tpool_next_result(q))) {
        struct timeval now;
        struct timespec timeout;
        gettimeofday(&now, NULL);
        timeout.tv_sec = now.tv_sec + 10;
        timeout.tv_nsec = now.tv_usec * 1000;
        q->ref_count++;
        if (q->shutdown) {
            int rc = --q->ref_count;
            pthread_mutex_unlock(&p->tpool_mu);
            if (rc == 0) tpool_process_destroy(q);
            return NULL;
        }     
        pthread_cond_timedwait(&q->output_avail_c, &q->p->tpool_mu, &timeout);
        q->ref_count--;
    }
    pthread_mutex_unlock(&p->tpool_mu);
    return r;
}

bool tpool_process_empty(tpool_process_t *q)
{
    bool empty;
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    empty = q->n_job == 0 && q->n_processing == 0 && q->n_result == 0;
    pthread_mutex_unlock(&p->tpool_mu);
    return empty;
}

static void tpool_process_shutdown(tpool_process_t *q)
{
    tpool_t *p = q->p;
    pthread_mutex_lock(&p->tpool_mu);
    q->shutdown = true;
    pthread_cond_broadcast(&q->output_avail_c);
    pthread_cond_broadcast(&q->input_not_full_c);
    pthread_cond_broadcast(&q->input_empty_c);
    pthread_cond_broadcast(&q->none_processing_c);
    pthread_mutex_unlock(&p->tpool_mu);
}


bool tpool_process_is_shutdown(tpool_process_t *q)
{
    tpool_t *p = q->p;
    bool shutdown;
    pthread_mutex_lock(&p->tpool_mu);
    shutdown = q->shutdown;
    pthread_mutex_unlock(&p->tpool_mu);
    return shutdown;
}

void tpool_delete_result(tpool_result_t *r, int free_data) 
{
    if (!r)
        return;

    if (free_data && r->data)
        free(r->data);

    free(r);
}

void *tpool_result_data(tpool_result_t *r) 
{
    return r->data;
}

