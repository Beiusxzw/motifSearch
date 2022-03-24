// ****************************************
// C Interface for bedtools' Fasta Index
// ----------------------------------------

#include "fasta.h"
#include "utils.h"


FastaIndex *fastaIndexInit()
{
    FastaIndex *fi = (FastaIndex *) malloc (sizeof(FastaIndex));
    fi->name_field = kh_init(str_hash_t);
    kv_init(fi->sequence_names);
    return fi;
}

FastaIndexEntry *fastaIndexEntryInit()
{
    FastaIndexEntry *entry = (FastaIndexEntry *) malloc (sizeof(FastaIndexEntry));
    memset(entry, 0, sizeof(FastaIndexEntry));
    entry->offset = -1;
    return entry;
}

FastaIndexEntry *fastaIndexEntryInitData(char* name, uint64_t length, uint64_t offset, uint64_t line_blen, uint64_t line_len, bool full_header) {
    FastaIndexEntry *entry = (FastaIndexEntry *) malloc (sizeof(FastaIndexEntry));
    if (name[strlen(name)-1]=='\n') {
        name[strlen(name)-1]='\0';
    }
    entry->name = name;
    entry->length = length;
    entry->line_len = line_len;
    entry->line_blen = line_blen;
    entry->offset = offset;
    entry->full_header = full_header;
    return entry;
}

void fastaIndexEntryEmpty(FastaIndexEntry* entry)
{
    free(entry->name);
    memset(entry, 0, sizeof(FastaIndexEntry));
    entry->offset = -1;
}

void fastaIndexEntryDestory(FastaIndexEntry* entry)
{
    free(entry);
}

void fastaIndexDestory(FastaIndex *fi)
{
    FastaIndexEntry *v;
    kh_foreach_value(fi->name_field, v, {
        free(v);
    });
    kh_destroy(str_hash_t, fi->name_field);
    kv_destroy(fi->sequence_names);
    free(fi);
}

stringVec split(const char *_s, char *delim)
{
    uint8_t i;
    stringVec fields;
    kv_init(fields);
    char* tok;
    char* s = malloc(strlen(_s));
    memcpy(s, _s, strlen(_s));
    tok = strtok(s, delim);
    while (tok != NULL) {
        if (tok[strlen(tok)-1] == '\n') tok[strlen(tok)-1] = '\0';
        kv_push(char *, fields, tok);
        tok = strtok(NULL, delim);
    }
    return fields;
}

void *writeFastaIndex(char* fasta_file_path, bool full_header, bool return_index)
{
    int fasta_file_path_len = strlen(fasta_file_path);
    char *index_file_path = malloc(strlen(fasta_file_path)+4);
    strcpy(index_file_path, fasta_file_path);
    strcpy(index_file_path + fasta_file_path_len, ".fai");
    FastaIndexEntry *entry = fastaIndexEntryInit();
    entry->full_header = full_header;
    FastaIndex *fi = fastaIndexInit();
    char *line = NULL;
    int64_t line_length;
    int64_t offset = 0;
    int64_t line_number;
    bool mismatched_line_len = false;
    bool empty_line = false;
    FILE* fp;
    FILE* op;
    size_t bufsize = 0;
    if (!(fp = fopen(fasta_file_path, "r"))) fatalf("Error: could not open fasta file %s\n", fasta_file_path);
    if (!(op = fopen(index_file_path, "w+"))) fatalf("Error: could not open fasta index file for writing %zu\n", index_file_path);
    while ((line_length = getline(&line, &bufsize, fp)) != -1) {
        line_length -= 1;
        ++line_number;
        if (line[0] == ';') {
            // fasta comment and skip
        } else if (line[0] == '+') {
            // fasta quality header
            line_length = getline(&line, &bufsize, fp);
            line_length -= 1;
            offset += line_length + 1;
            line_length = getline(&line, &bufsize, fp);
            line_length -= 1;
        } else if (line[0] == '>' || line[0] == '@') {
            // if we aren't on the first entry, push the last sequence into the index
            if (entry->name != "" && entry->name != NULL) {
                mismatched_line_len = false;
                empty_line = false;
                entryToIndex(entry, fi);
                fastaIndexEntryEmpty(entry);
            }
            
            char *new_name = malloc(line_length+1);
            strcpy(new_name, line + 1);
            entry->name = new_name;
            if (entry->name[strlen(entry->name)-1]=='\n') {
                entry->name[strlen(entry->name)-1]='\0';
            }
        } else {
            // assume we have a sequence file
            if (entry->offset == -1) {
                entry->offset = offset;
            }
            entry->length += line_length;
            if (entry->line_len) {
                if (mismatched_line_len || empty_line) {
                    if (line_length == 0) {
                        empty_line = true;
                    }
                    else {
                        if (empty_line) {
                            fatal("Error: found an empty line\n");
                        } else {
                            fatal("Error: mismatched line length\n");
                        }
                        fatal("Error: in generating index file\n");
                    }
                }
                if (entry->line_len != line_length + 1) {
                    mismatched_line_len = true;
                    if (line_length == 0) {
                        empty_line = true;
                    }
                } 
            } else {
                entry->line_len = line_length + 1;
            }
            entry->line_blen = entry->line_len - 1;
        }
        offset += line_length + 1;
    }
    entryToIndex(entry, fi);
    fastaIndexEntryEmpty(entry);
    indexToFile(fi, op);
    fclose(op);
    if (return_index) {
        return fi;
    } else {
        fastaIndexDestory(fi);
    }
}

void entryToIndex(FastaIndexEntry *entry, FastaIndex *fi)
{
    char *_name = malloc(strlen(entry->name)+1);
    memset(_name, 0, strlen(entry->name)+1);
    strcpy(_name, entry->name);
    FastaIndexEntry *new_entry = fastaIndexEntryInitData(_name, entry->length, entry->offset, entry->line_blen, entry->line_len, entry->full_header);
    int ret;
    khiter_t i = kh_put(str_hash_t, fi->name_field, new_entry->name, &ret);
    if (!ret) kh_del(str_hash_t, fi->name_field, new_entry->name);
    kh_value(fi->name_field, i) = new_entry;
}

void indexToFile(FastaIndex *fi, FILE* op)
{
    // TODO: The index has not been sorted!
    FastaIndexEntry *entry;
    kh_foreach_value(fi->name_field, entry, {
        entryToFile(entry, op);
    })
}

void indexToStdout(FastaIndex *fi)
{
    FastaIndexEntry *entry;
    // TODO: The index has not been sorted!
    kh_foreach_value(fi->name_field, entry, {
        printf("%s\t%lld\t%lld\t%lld\t%lld\n", entry->name, entry->length, entry->offset, entry->line_blen, entry->line_len);
    })
}

void entryToFile(FastaIndexEntry *entry, FILE* op)
{
    if (entry->full_header) {
        fprintf(op, "%s\t%lld\t%lld\t%lld\t%lld\n",entry->name, entry->length, entry->offset, entry->line_blen, entry->line_len);
    } else {
        char buf[100];
        memset(buf, 0, 100);
        for (int i = 0; i < strlen(entry->name); ++i) {
            if (entry->name[i] != ' ') {
                buf[i] = entry->name[i];
            }
        }
        fprintf(op, "%s\t%lld\t%lld\t%lld\t%lld\n", buf, entry->length, entry->offset, entry->line_blen, entry->line_len);
    }
}

FastaIndex *readFastaIndex(char* index_file_path, bool full_header) 
{
    FastaIndex *fi = fastaIndexInit();
    char* line = NULL;
    uint64_t line_num;
    uint64_t line_len;
    FILE* fp;
    int ret;
    if (!(fp = fopen(index_file_path, "rb"))) {
        printf("Error: could not open fasta index file %s\n", index_file_path);
        return NULL;
    }
    while (ret = getline(&line, &line_len, fp) != -1) {
        ++line_num;
        stringVec fields = split(line, "\t");
        if (kv_size(fields) == 5) {
            char* name = kv_A(fields, 0);
            uint64_t end;
            uint64_t length;
            uint64_t line_blen;
            uint64_t line_len;
            FastaIndexEntry *entry = fastaIndexEntryInitData(
                name, 
                strtoll(kv_A(fields, 1), &length, 10), 
                strtoll(kv_A(fields, 2), &end, 10), 
                strtoll(kv_A(fields, 3), &line_blen, 10), 
                strtoll(kv_A(fields, 4), &line_len, 10), full_header);
            kv_push(char*, fi->sequence_names, name);
            entryToIndex(entry, fi);
        } else {
            fatalf("Warning: malformed fasta index file %s\n", index_file_path);
        }
    }
    return fi;
}

void *readFastaByMmap(char* fasta_file_path)
{
    FILE* fp;
    if (!(fp = fopen(fasta_file_path, "r"))) fatalf("Error: could not open fasta file %s\n", fasta_file_path);
    int fd = fileno(fp); // Get file discriptor for mmap
    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        fatal("Failed to stat the file\n");
    }
    size_t filesize= sb.st_size;
    void *filemm = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    struct fmm *ret = malloc(sizeof(struct fmm));
    ret->mm = filemm;
    ret->fs = filesize;
    return ret;
}

char *getFastaSequenceMmap(void *filemm, FastaIndex *fi, char *seq_name)
{
    khiter_t k = kh_get(str_hash_t, fi->name_field, seq_name);
    if (k == kh_end(fi->name_field)) fatalf("%s not found in index", seq_name);
    FastaIndexEntry *entry = kh_val(fi->name_field, k);
    int64_t newlines_in_sequence = entry->length / entry->line_blen;
    int64_t seqlen = newlines_in_sequence + entry->length;
    char* seq = (char *)calloc(seqlen + 1, 1);
    memcpy(seq, (char *) filemm + entry->offset, seqlen);
    seq[seqlen] = '\0'; // end of string char
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    if (*pbegin == '\n' || *pbegin == '\0') {
        seq = seq + 1;
    }
    if (*pend == '\n') {
        *pend = '\0';
    }
    char* ret = (char *)calloc(seqlen + 1, 1);
    size_t t = 0;
    for (size_t i = 0; i < seqlen; ++i) {
        if (seq[i] != '\n') {
            ret[t] = seq[i];
            t++;
        }
    }
    free(seq);
    return ret;
}

char *getFastaSequenceMmap2(void *filemm, FastaIndexEntry *entry)
{
    int64_t newlines_in_sequence = entry->length / entry->line_blen;
    int64_t seqlen = newlines_in_sequence + entry->length;
    char* seq = (char *)calloc(seqlen + 1, 1);
    memcpy(seq, (char *) filemm + entry->offset, seqlen);
    seq[seqlen] = '\0'; // end of string char
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    if (*pbegin == '\n' || *pbegin == '\0') {
        seq = seq + 1;
    }
    if (*pend == '\n') {
        *pend = '\0';
    }
    char* ret = (char *)calloc(seqlen + 1, 1);
    size_t t = 0;
    for (size_t i = 0; i < seqlen; ++i) {
        if (seq[i] != '\n') {
            ret[t] = seq[i];
            t++;
        }
    }
    free(seq);
    return ret;
}


/*
int main(int argc, char const *argv[])
{

    FastaIndex* fi = readFastaIndex("/Users/snow/Documents/refData/fasta/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa.fai", 0);
    struct fmm *m = readFastaByMmap("/Users/snow/Documents/refData/fasta/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa");
    printf("%s\n", getFastaSequenceMmap(m->mm, fi, "KI270392.1"));
    munmap(m->mm, m->fs);
    free(m);
    fastaIndexDestory(fi);


    // writeFastaIndex("/Users/snowxue/Desktop/Homo_sapiens.GRCh38.dna.primary_assembly.fa", 0);  

}

*/