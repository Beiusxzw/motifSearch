CFLAGS=	 -g -Wall -Wc++-compat -w #-Wextra
INCLUDES=
OBJS=
PROG= motifSearch
PROG_EXTRA= supp
LIBS=	 -lm -lz -lpthread
HEADERS := $(wildcard *.h) $(wildcard $(AHOCORASICK_DIR)/includes/*.h)
AHOCORASICK_DIR= ./ahocorasick/src
SHARED_CS= motifSearch.c thread_pool.c

ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif


ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread
endif

.PHONY:all extra clean depend
.SUFFIXES:.c .o

.c.o: $(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(HEADERS) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

motifSearch:main.o fasta.o motifSearch.o thread_pool.o ahocorasick.a
	$(CC) $(CFLAGS) main.o fasta.o motifSearch.o thread_pool.o ahocorasick.a -o $@ -L. $(LIBS)

clean:
	rm -fr *.o a.out $(PROG_EXTRA) *~ *.a *.dSYM build dist mappy*.so mappy.c python/mappy.c mappy.egg*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

main.o:main.c $(SHARED_CS) $(HEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
fasta.o: fasta.c $(SHARED_CS) $(HEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
motifSearch.o: motifSearch.c $(SHARED_CS) $(HEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
thread_pool.o: thread_pool.c $(SHARED_CS) $(HEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


ahocorasick.a: aho_queue.o aho_trie.o ahocorasick.o
		$(AR) -csru $@  aho_queue.o aho_trie.o ahocorasick.o
aho_queue.o: $(AHOCORASICK_DIR)/aho_queue.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
aho_trie.o: $(AHOCORASICK_DIR)/aho_trie.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
ahocorasick.o: $(AHOCORASICK_DIR)/ahocorasick.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@