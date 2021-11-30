// ****************************************
// C Interface for bedtools' Fasta Index
// ----------------------------------------

#ifndef _FASTA_H
#define _FASTA_H

#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdbool.h>
#include "khash.h"
#include "kvec.h"

KHASH_MAP_INIT_STR(str_hash_t, const void *);
typedef kvec_t(char *) stringVec;
typedef khash_t(str_hash_t) stringMap;

struct fmm {
    void *mm;
    size_t fs;
};

typedef struct FastaIndexEntry {
    char* name;         // Name of the fasta sequence
    int64_t length;     // length of the sequences bytes stored
    int64_t offset;     // Offset to the begining of the file
    int64_t line_blen;  // line length in bytes, sequence characters
    int64_t line_len;   // line length including newline
    bool full_header;
} FastaIndexEntry;

FastaIndexEntry *fastaIndexEntryInit();
FastaIndexEntry *fastaIndexEntryInitData(char* name, uint64_t length, uint64_t offset, uint64_t line_blen, uint64_t line_len, bool full_header);
void fastaIndexEntryDestory(FastaIndexEntry* entry);


typedef struct FastaIndex {
    stringVec sequence_names; // For sorting the FastaIndexEntries
    stringMap *name_field;    // HashMap of <name, FastaIndexEntry>
    bool full_header;         // Whether to remove the '>' of the name
} FastaIndex;

FastaIndex *FastaIndexInit();
void fastaIndexDestory(FastaIndex *fi);

stringVec split(const char *_s, char *delim);
void *writeFastaIndex(char* fasta_file_path, bool full_header, bool return_index);
void entryToIndex(FastaIndexEntry *entry, FastaIndex *fi);
void indexToFile(FastaIndex *fi, FILE* op);
void entryToFile(FastaIndexEntry *entry, FILE* op);
void indexToStdout(FastaIndex *fi);
void *readFastaByMmap(char* fasta_file_path);
char *getFastaSequenceMmap(void *filemm, FastaIndex *fi, char *seq_name);
char *getFastaSequenceMmap2(void *filemm, FastaIndexEntry *entry);
FastaIndex *readFastaIndex(char* index_file_path, bool full_header) ;



#endif