/*  thread_pool.h -- API for the thread pool.

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

#ifndef THREAD_POOL_H
#define THREAD_POOL_H
#define _XOPEN_SOURCE 600
#include <pthread.h>
#include <stdint.h>
#include <ucontext.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tpool tpool_t;
typedef struct tpool_process tpool_process_t;
typedef struct tpool_worker tpool_worker_t;
typedef struct tpool_job tpool_job_t;
typedef struct tpool_result tpool_result_t;

/**
 * < Implementation of the thread pool >
 */

/**
 * @brief tpool data structure.
 */
struct tpool {
    int nwaiting; /* How many workers waiting for new jobs */
    int njobs;    /* Howmany jobs are waiting in all queues */
    int shutdown;

    tpool_process_t *q_head;

    /* threads */
    int tsize; /* maximum number of jobs */
    tpool_worker_t *t;
    /* array of worker IDs free */
    int *t_stack, t_stack_top;

    pthread_mutex_t tpool_mu;
    /* Tracking of average number of running jobs */
    int n_count, n_running;

    long long total_time, wait_time;
};

/**
 * @brief 
 * 
 */
struct tpool_process {
    tpool_t *p;
    /* The input job forms a linked-list */
    tpool_job_t *job_head;
    tpool_job_t *job_tail;
    /* The output result forms a linked-list */
    tpool_result_t *result_head;
    tpool_result_t *result_tail;

    int qsize;    /* max size of IO queues */
    uint64_t next_serial; /* next serial for output */
    uint64_t curr_serial; /* current serial (next input) */

    bool no_more_input; /* disable dispatching of more jobs */
    int n_job; /* no. items in input queue; was njobs */
    int n_result; /* no. items in output queue; was nresult */
    int n_processing; /* no. items being proessed */

    bool shutdown; /* true if the thread pool is destroyed */
    bool in_only;  /* if true, don't queue result up */
    bool wake_dispatch; /* unblocks waiting dispatchers */

    int ref_count; /* reference count. used to track destruction */

    pthread_cond_t output_avail_c;   /* Signalled on each new output */
    pthread_cond_t input_not_full_c; /* Input queue is no longer full */
    pthread_cond_t input_empty_c;    /* Input queue has become empty */
    pthread_cond_t none_processing_c;/* n_processing has hit zero */

    tpool_process_t *next, *prev; /* the process forms a circular linked-list */
};

/**
 * @brief The minimal wrapper for a thread
 * 
 */
struct tpool_worker {
    tpool_t *p;
    int idx;
    pthread_t tid;
    pthread_cond_t pending_c;
};

/**
 * @brief 
 * 
 */
struct tpool_job {
    void *(*func)(void *arg);
    void *arg;
    void (*job_cleanup)(void *arg); /* callback function for job */
    void (*result_cleanup)(void *data); /* callback function for result */
    tpool_job_t *next;

    tpool_t *p;
    tpool_process_t *q;
    uint64_t serial; /* serial number for ordering */
};

struct tpool_result {
    void (*result_cleanup)(void *data);
    void *data;
    uint64_t serial; /* serial number for ordering */
    tpool_result_t *next;
};


tpool_process_t *tpool_process_init(tpool_t *p, int qsize, bool in_only);
void tpool_process_destroy(tpool_process_t *q);
tpool_t* tpool_init(int n);
void tpool_destroy(tpool_t *p);
static void *tpool_worker(void *arg);
static int tpool_add_result(tpool_job_t *j, void *data);
void tpool_process_attach(tpool_t *p, tpool_process_t *q);
void tpool_process_detach(tpool_t *p, tpool_process_t *q);
int tpool_process_reset(tpool_process_t *q, bool _free);
int tpool_process_flush(tpool_process_t *q);
void tpool_delete_result(tpool_result_t *r, int free_data);
static void tpool_process_shutdown(tpool_process_t *q);
static void wake_next_worker(tpool_process_t *q, bool locked);
int tpool_dispatch(tpool_t *p, 
                   tpool_process_t *q, 
                   void *(*func)(void *arg),
                   void *arg,
                   void (*job_cleanup)(void *arg),
                   void (*result_cleanup)(void *data),
                   bool nonblock);
tpool_result_t *tpool_next_result(tpool_process_t *q);
tpool_result_t *tpool_next_result_wait(tpool_process_t *q);
bool tpool_process_empty(tpool_process_t *q);
static void tpool_process_shutdown(tpool_process_t *q);
bool tpool_process_is_shutdown(tpool_process_t *q);
void hts_tpool_delete_result(tpool_result_t *r, int free_data);
void *hts_tpool_result_data(tpool_result_t *r);

/**
 * </ Implementation of the thread pool >
 */

#ifdef __cplusplus
}
#endif
#endif

