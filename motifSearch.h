#ifndef _MOTIF_SEARCH
#define _MOTIF_SEARCH

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "kseq.h"
#include "utils.h"
#include "fasta.h"
#include "./ahocorasick/include/ahocorasick.h"

#define MAX_MOTIF_LEN 64
#define MAX_PATTERN_LEN 512

KSEQ_INIT(FILE*, read)

struct pt_info {
	int motif_len;
	char* chrom;
	pthread_mutex_t *mu;
	char *seq;
};

struct par_arg {
	char* chrom;
    char* file_path;
	FastaIndexEntry *entry;
	const char** pattern;
	int n_patterns;
	int motif_len;
	pthread_mutex_t *pt_mu;
	int n_threads;
};

void search_motif(struct ahocorasick *aho, const char* seq, const char* chrom, int motif_len, bool uselock, pthread_mutex_t *mu);
void init_ahocorasick(struct ahocorasick *aho, const char** pattern, int n_patterns);
void search_fasta(const char** file_path, const char** pattern, int n_patterns, int motif_len);
int parse_motif_pattern(char* motif, char** pattern);
void search_fasta_par(void *arg);
void search_fasta_par_test(void *arg);
void free_par_arg(void *arg);

#endif