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

char* parse_iupac(char c)
{
	char* new_char = malloc(5);
	switch (c) {
		case 'M':
			memcpy(new_char, "AC", 2); break;
		case 'R':
			memcpy(new_char, "AG", 2); break;
		case 'W':
			memcpy(new_char, "AT", 2); break;
		case 'S':
			memcpy(new_char, "CG", 2); break;
		case 'Y':
			memcpy(new_char, "CT", 2); break;
		case 'K':
			memcpy(new_char, "GT", 2); break;
		case 'V':
			memcpy(new_char, "ACG", 2); break;
		case 'H':
			memcpy(new_char, "ACT", 2); break;
		case 'D':
			memcpy(new_char, "AGT", 2); break;
		case 'B':
			memcpy(new_char, "CGT", 2); break;
		case 'N':
			memcpy(new_char, "AGCT", 2); break;
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
			*index = ind;
		}
	} else if (is_valid_ambiguity_codes(motif, strlen(motif))) {
		char* tmp = malloc(sizeof(motif)+1);
		memcpy(tmp, motif, strlen(motif));
		for (int i = 0; i < strlen(motif); ++i) {
			char* amb = parse_iupac(motif[i]);
			if (amb) {
				for (int j = 0; j < strlen(amb); ++j) {
					tmp[i] = amb[j];
					parse_motif_pattern_help(tmp, pattern, index);
				}
				free(amb);
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
	return ind+1;
}