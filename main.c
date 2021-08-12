#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include "thread_pool.h"
#include "motifSearch.h"

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX_THREADS  sysconf(_SC_NPROCESSORS_ONLN)
#define MAX_CHROM 100

void version()
{
    printf("Snow's motifSearch version 0.0.1\n");
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

int main(int argc, char const *argv[])
{
    static bool verbose_flag;
    int n_threads = 0;
    char *file_path = NULL;
    char *motif = NULL;
    int c;
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
    struct chrom_seq* chrom_seqs[MAX_CHROM];

    if (n_threads > 1) {
        tpool_t *p = tpool_init(n_threads);
        tpool_process_t *q = tpool_process_init(p, 16, true);
        int num_chrom = read_fasta(file_path, chrom_seqs);
        
        for (int i = 0; i < num_chrom; ++i) {
            int blk;
            struct par_arg *arg = malloc(sizeof(struct par_arg));
            arg->chrom_seq = chrom_seqs[i];
            arg->motif_len = strlen(motif);
            arg->n_patterns = num;
            arg->pattern = pattern;
            do {
                blk = tpool_dispatch(p, q, &search_fasta_par, arg, NULL, &free_par_arg, true);
                if (blk == -1) {
                    usleep(10000);
                }
            } while (blk == -1);
        }
        tpool_process_flush(q);
        tpool_process_destroy(q);
        tpool_destroy(p);
    } else {
        search_fasta(file_path, pattern, num, strlen(motif));
        for (int i = 0; i < num; ++i)
        {
            free(pattern[i]);
        }
        return 0;
    }
}
