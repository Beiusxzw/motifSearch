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


void *doit_square_u(void *arg) {

    int job = *(int *)arg;
    int *res;
    res = malloc(sizeof(*res));
    job = (job<0) ? -job*job : job*job;
    *res = (job<0) ? -job*job : job*job;
    for (int i = 0; i < 100000000; ++i) {
        *res += 1;
    }
    free(arg);
    return res;
}


int main(int argc, char const *argv[])
{
    /* Test */
    pthread_setconcurrency(2);
    tpool_t *p = tpool_init(8);
    tpool_process_t *q = tpool_process_init(p, 16, false);
    tpool_result_t *r;
    for (int i = 0; i < 100; ++i) {
        int *ip = malloc(sizeof(*ip));
        *ip = i;
        int blk;
        do {
            blk = tpool_dispatch(p, q, doit_square_u, ip, NULL, NULL, true);
            if (blk == -1) {
                usleep(10000);
            }
            while (r = tpool_next_result(q)) {
                printf("RESULT: %d\n", *(int *)r->data);
                tpool_delete_result(r, true);
            }
        } while (blk == -1);
    }
    tpool_process_flush(q);
    while (r = tpool_next_result(q)) {
        printf("RESULT: %d\n", *(int *)r->data);
        tpool_delete_result(r, true);
    }
    tpool_process_destroy(q);
    tpool_destroy(p);
    pthread_exit(NULL);
}

