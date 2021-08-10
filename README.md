# motifSearch
Snow's motif search in C (A really simplified command, expecting to be **5×** faster than homer's scanMotifGenome in single-thread mode)


Motif searching supporting [IUPAC ambiguous codes](https://droog.gs.washington.edu/mdecode/images/iupac.html).


## Usage 

 ```
 motifSearch -f <FASTA> -m <MOTIF> -p <THREAD> > <OUTPUT-BED>
-f/--fasta      fasta file
-m/--motif      motif string
-p/--nthreads   number of threads
```

To re-compile, `gcc main.c motifSearch.c ahocorasick/src/*c -w -o motifSearch`

## TODO

- [ ] stats compare to other methods
- [ ] multi-threaded searching
- [ ] mismatch handling
