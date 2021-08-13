#include <stdbool.h>
#include "motifSearch.h"

/* Another array to help us do complement of char */
static char ntCompTable[256];
static bool inittedCompTable = false;

void initNtCompTable(void)
{
    ntCompTable[' '] = ' ';
    ntCompTable['-'] = '-';
    ntCompTable['='] = '=';
    ntCompTable['a'] = 't';
    ntCompTable['c'] = 'g';
    ntCompTable['g'] = 'c';
    ntCompTable['t'] = 'a';
    ntCompTable['u'] = 'a';
    ntCompTable['n'] = 'n';
    ntCompTable['-'] = '-';
    ntCompTable['.'] = '.';
    ntCompTable['A'] = 'T';
    ntCompTable['C'] = 'G';
    ntCompTable['G'] = 'C';
    ntCompTable['T'] = 'A';
    ntCompTable['U'] = 'A';
    ntCompTable['N'] = 'N';
    ntCompTable['R'] = 'Y';
    ntCompTable['Y'] = 'R';
    ntCompTable['M'] = 'K';
    ntCompTable['K'] = 'M';
    ntCompTable['S'] = 'S';
    ntCompTable['W'] = 'W';
    ntCompTable['V'] = 'B';
    ntCompTable['H'] = 'D';
    ntCompTable['D'] = 'H';
    ntCompTable['B'] = 'V';
    ntCompTable['X'] = 'N';
    ntCompTable['r'] = 'y';
    ntCompTable['y'] = 'r';
    ntCompTable['s'] = 's';
    ntCompTable['w'] = 'w';
    ntCompTable['m'] = 'k';
    ntCompTable['k'] = 'm';
    ntCompTable['v'] = 'b';
    ntCompTable['h'] = 'd';
    ntCompTable['d'] = 'h';
    ntCompTable['b'] = 'v';
    ntCompTable['x'] = 'n';
    inittedCompTable = true;
}

static char ntMixedCaseChars[256];

void initNtMixedCaseChars(void)
{
    static bool initted = false;

    if (!initted)
    {
        ntMixedCaseChars['a'] = 'a';
        ntMixedCaseChars['A'] = 'A';
        ntMixedCaseChars['c'] = 'c';
        ntMixedCaseChars['C'] = 'C';
        ntMixedCaseChars['g'] = 'g';
        ntMixedCaseChars['G'] = 'G';
        ntMixedCaseChars['t'] = 't';
        ntMixedCaseChars['T'] = 'T';
        ntMixedCaseChars['n'] = 'n';
        ntMixedCaseChars['N'] = 'N';
        ntMixedCaseChars['u'] = 'u';
        ntMixedCaseChars['U'] = 'U';
        ntMixedCaseChars['-'] = 'n';
        initted = true;
    }
}

/* Reverse the order of the bytes. */
void reverseBytes(char *bytes, long length)
{
    long halfLen = (length >> 1);
    char *end = bytes + length;
    char c;
    while (--halfLen >= 0)
    {
        c = *bytes;
        *bytes++ = *--end;
        *end = c;
    }
}

void complement(char *dna, long length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        dna[i] = ntCompTable[(int)dna[i]];
    }
}

char *reverseComplement(char *dna, long length)
{
	char* new_dna = malloc(strlen(dna));
	memcpy(new_dna, dna, strlen(dna));
    reverseBytes(new_dna, length);
    complement(new_dna, length);
    return new_dna;
}

/* A little array to help us decide if a character is a 
 * nucleotide, and if so convert it to lower case. */
char ntChars[256];
static bool initted = false;
static void initNtChars()
{

    if (!initted)
    {
        memset(ntChars, 0, sizeof(ntChars));
        ntChars['a'] = ntChars['A'] = 'a';
        ntChars['c'] = ntChars['C'] = 'c';
        ntChars['g'] = ntChars['G'] = 'g';
        ntChars['t'] = ntChars['T'] = 't';
        ntChars['u'] = ntChars['U'] = 'u';
        initted = true;
    }
}

bool is_valid_dna(char *poly, int size)
{
    if (!initted)
        initNtChars();
    int i;
    int dnaCount = 0;

    for (i = 0; i < size; ++i)
    {
        if (!ntChars[(int)poly[i]])
            return false;
    }
    return true;
}

bool is_valid_ambiguity_codes(char *poly, int size)
{
    if (!initted)
        initNtChars();
    int i;
    int dnaCount = 0;

    for (i = 0; i < size; ++i)
    {
        if (!ntCompTable[(int)poly[i]])
            return false;
    }
    return true;
}

void upper_str(char *__s, size_t __size)
{
    for (size_t i = 0; i < __size; i++)
    {
        __s[i] = toupper(__s[i]);
    }
}

void aho_callback(void *arg, struct aho_match_t *m)
{
	struct motif_info *t = (struct tmp *)arg;
	printf("%s\t%lu\t%lu\n", t->chrom, m->pos, m->pos + t->motif_len);
}

void search_motif(struct ahocorasick *aho, const char* seq, const char* chrom, int motif_len)
{
	struct motif_info *arg;
	if (!(arg = malloc(sizeof(struct motif_info)))) {
		fatal("Memory allocation failed");
	}
	arg->chrom = chrom;
	arg->motif_len = motif_len;
	aho_register_match_callback(aho, &aho_callback, (void *)arg);
	aho_findtext(aho, seq, strlen(seq));
}

void init_ahocorasick(struct ahocorasick *aho, const char** pattern, int n_patterns)
{
	aho_init(aho);
    /** 
    for (int i = 0; i < n_patterns; ++i)
    {
        printf("# Pattern: %s\n", pattern[i]);
    }
    **/
	for (int i = 0; i < n_patterns; i++)
	{
		aho_add_match_text(aho, pattern[i], strlen(pattern[i]));
	}
	aho_create_trie(aho);
}

void search_fasta(const char** file_path, const char** pattern, int n_patterns, int motif_len)
{
	kseq_t *seq;
	FILE* fp;
    int n = 0, slen = 0, qlen = 0;
    if (!(fp = fopen(file_path, "r"))) {
		fatal("File open failed\n");
	}
    seq = kseq_init(fileno(fp));
	struct ahocorasick aho;
	init_ahocorasick(&aho, pattern, n_patterns);
	while (kseq_read(seq) >= 0)
    {
        upper_str(seq->seq.s, strlen(seq->seq.s));
		search_motif(&aho, seq->seq.s, seq->name.s, motif_len);
        ++n;
    }
	aho_destroy(&aho);
	kseq_destroy(seq);
    fclose(fp);
}

int read_fasta(const char** file_path, struct chrom_seq **chrom_seqs)
{
	kseq_t *seq;
	FILE* fp;
    int n = 0, slen = 0, qlen = 0;
    if (!(fp = fopen(file_path, "r"))) {
		fatal("File open failed\n");
	}
    seq = kseq_init(fileno(fp));
    while (kseq_read(seq) >= 0)
    {
        upper_str(seq->seq.s, strlen(seq->seq.s));
        chrom_seqs[n] = malloc(sizeof(struct chrom_seq));
        chrom_seqs[n]->seq = malloc(seq->seq.l+1);
        chrom_seqs[n]->chrom = malloc(seq->name.l+1);
		memcpy(chrom_seqs[n]->seq, seq->seq.s, seq->seq.l);
        memcpy(chrom_seqs[n]->chrom, seq->name.s, seq->name.l);
        ++n;
    }
    kseq_destroy(seq);
    fclose(fp);
    return n;
}

void search_fasta_par(void *arg)
{
    struct par_arg *parg = (struct par_arg *)arg;
    char* pattern = parg->pattern;
    int n_patterns = parg->n_patterns;
    struct chrom_seq *chrom_seq = parg->chrom_seq;
    int motif_len = parg->motif_len;
    struct ahocorasick aho;
	init_ahocorasick(&aho, pattern, n_patterns);
    search_motif(&aho, chrom_seq->seq, chrom_seq->chrom, motif_len);
    aho_destroy(&aho);
}

void free_par_arg(void *arg)
{
    struct par_arg *parg = (struct par_arg *)arg;
    free(parg->chrom_seq->chrom);
    free(parg->chrom_seq->seq);
    free(parg->chrom_seq);
    free(parg);
}

char* parse_iupac(char c)
{
	char* new_char = malloc(5);
	switch (c) {
		case 'M':
			return "AC";
		case 'R':
			return "AG";
		case 'W':
			return "AT";
		case 'S':
			return "CG";
		case 'Y':
			return "CT";
		case 'K':
			return "GT";
		case 'V':
			return "ACG";
		case 'H':
			return "ACT";
		case 'D':
			return "AGT";
		case 'B':
			return "CGT";
		case 'N':
			return "AGCT";
		default:
			free(new_char);
			new_char = NULL;
	}
	return new_char;
}

int parse_motif_pattern_help(char* motif, char** pattern, int* index)
{
	int ind = *index;
	if (is_valid_dna(motif, strlen(motif))) {
		bool exist = false;
		for (int n = 0; n < ind; n++) {
			if (strncmp(pattern[n], motif, strlen(motif)) == 0)
			{
				exist = true;
			}
		}
		// simple DNA motif, no ambiguity codes.
		if (!exist) {
			pattern[ind] = malloc(strlen(motif)+1);
			memcpy(pattern[ind], motif, strlen(motif));
			pattern[++ind] = reverseComplement(motif, strlen(motif));
			*index = ind+1;
		}
	} else if (is_valid_ambiguity_codes(motif, strlen(motif))) {
		char* tmp = malloc(sizeof(motif)+1);
		memcpy(tmp, motif, strlen(motif));
		for (int i = 0; i < strlen(motif); ++i) {
			char* amb = parse_iupac(tmp[i]);
			if (amb) {
				for (int j = 0; j < strlen(amb); ++j) {
					tmp[i] = amb[j];
					parse_motif_pattern_help(tmp, pattern, index);
				}
			}
		}
		free(tmp);
	} else {
		fatal("Invalid ambiguity codes!");
	}
}

int parse_motif_pattern(char* motif, char** pattern)
{
	int ind = 0;
	initNtChars();
	initNtCompTable();
	parse_motif_pattern_help(motif, pattern, &ind);
	return ind;
}