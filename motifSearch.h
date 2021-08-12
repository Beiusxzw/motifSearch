#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"
#include "utils.h"
#include "./ahocorasick/include/ahocorasick.h"

#define MAX_MOTIF_LEN 25
#define MAX_PATTERN_LEN 80

KSEQ_INIT(FILE*, read)

struct motif_info{
	int motif_len;
	char* chrom;
};

struct chrom_seq {
    char* chrom;
    char* seq;
};

struct par_arg {
	struct chrom_seq *chrom_seq;
	const char** pattern;
	int n_patterns;
	int motif_len
};

void search_motif(struct ahocorasick *aho, const char* seq, const char* chrom, int motif_len);
void init_ahocorasick(struct ahocorasick *aho, const char** pattern, int n_patterns);
void search_fasta(const char** file_path, const char** pattern, int n_patterns, int motif_len);
int parse_motif_pattern(char* motif, char** pattern);
int read_fasta(const char** file_path, struct chrom_seq **chrom_seqs);
void search_fasta_par(void *arg);
void free_par_arg(void *arg);