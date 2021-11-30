# motifSearch
Current version: 0.0.2 == multi-threaded searching
Snow's motif search in C (A simple command-line-tool, expecting to be faster than homer's scanMotifGenome)


Motif searching supporting [IUPAC ambiguous codes](https://droog.gs.washington.edu/mdecode/images/iupac.html).


## Usage 

```
motifSearch -f <FASTA> -m <MOTIF> -p <THREAD> > <OUTPUT-BED>
-f/--fasta      fasta file
-m/--motif      motif string
-p/--nthreads   number of threads
```

To compile, `make && make clean`

## TODO

- [ ] stats compare to other methods
- [ ] mismatch handling
