#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include "thread_pool.h"
#include "motifSearch.h"
#include "fasta.h"

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX_THREADS  sysconf(_SC_NPROCESSORS_ONLN)
#define MAX_CHROM 100

void version()
{
    printf("Snow's motifSearch version 0.0.2\n");
    printf("Usage:\n");
    printf("\t-f/--fasta\tfasta file\n");
    printf("\t-m/--motif\tmotif string\n");
    printf("\t-p/--nthreads\tnumber of threads\n");
}

void usage()
{
    version();
    printf("");
}

void *test(void *arg) {
    struct par_arg *parg = (struct par_arg *) arg;
    printf("entry %s\n", parg->entry->name);
}


int main(int argc, char const *argv[])
{
    static bool verbose_flag;
    int n_threads = 0;
    char *file_path = NULL;
    char *motif = NULL;
    FastaIndex *fi;
    int c;
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutex_t pt_mu;
    pthread_mutex_init(&pt_mu, &attr);

    while (1)
    {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"verbose", no_argument, &verbose_flag, 1},
                {"brief", no_argument, &verbose_flag, 0},
                /* These options don’t set a flag.
             We distinguish them by their indices. */
                {"fasta", required_argument, 0, 'f'},
                {"motif", required_argument, 0, 'm'},
                {"nthreads", optional_argument, 0, 'p'},
                {"help", no_argument, NULL, 'h'},
                {"version", no_argument, NULL, 'v'},
                {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long(argc, argv, "f:hm:p:v", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;

        case 'f':
            file_path = optarg;
            break;

        case 'h':
            usage();
            exit(0);

        case 'm':
            motif = optarg;
            break;

        case 'p':
            n_threads = MIN(strtol(optarg, NULL, 10), MAX_THREADS) ;
            break;

        case 'v':
            version();
            exit(0);

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            fatal("Argument capture failed\n");
        }
    }
      
    n_threads = n_threads ? n_threads : MAX_THREADS;
    char *pattern[MAX_PATTERN_LEN];
    int num = parse_motif_pattern(motif, &pattern);


    pthread_setconcurrency(2);
    tpool_t *p = tpool_init(n_threads);
    tpool_process_t *q = tpool_process_init(p, 16, true);

    char *index_file_path = malloc(strlen(file_path)+4);
    strcpy(index_file_path, file_path);
    strcpy(index_file_path + strlen(file_path), ".fai");

    if (!(fi = readFastaIndex(index_file_path, 0))) {
        printf("No index file found. Generating index file...\n");
        fi = writeFastaIndex(file_path, 0, true);
    }
    int num_chrom = kv_size(fi->sequence_names);
    char *chrom_name;
    FastaIndexEntry *entry;
    kh_foreach(fi->name_field, chrom_name, entry, {
            int blk;
            struct par_arg *arg = malloc(sizeof(struct par_arg));
            /* full header or not ?*/
            char buf[100];
            memset(buf, 0, 100);
            for (int i = 0; i < strlen(entry->name); ++i) {
                if (entry->name[i] != ' ') {
                    buf[i] = entry->name[i];
                }
            }
            arg->chrom = buf;
            arg->file_path = file_path;
            arg->entry = entry;
            arg->n_patterns = num;
            arg->pattern = pattern;
            arg->pt_mu = &pt_mu;
            arg->n_threads = n_threads;
            do {
                blk = tpool_dispatch(p, q, search_fasta_par, (void *)arg, NULL, free_par_arg, true);
                if (blk == -1) {
                    usleep(10000);
                }
            } while (blk == -1);
    })

    tpool_process_flush(q);
    tpool_process_destroy(q);
    tpool_destroy(p);
    pthread_exit(NULL);
}
