#include <errno.h>
#include <setjmp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>

/* Optimization hints */
#if defined __GNUC__ || defined __llvm__
#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif 

/* Defines a unique identifier of type const void*. */
#define concat(x,y) x##y
#define void_uid(name) \
    static const int concat(name, ___) = 0;\
    const void *name = & concat(name, ___);

/**
 * @brief Takes a pointer to a member variable and computes pointer to the 
 * structure that contains it.
 * @param ptr ptr to the member
 * @param type type of the structure
 * @param member member name
 */
#define parent_struct(ptr, type, member) \
    (ptr ? ((type*) (((char*) ptr) - offsetof(type, member))) : NULL)

#define fatal(s) do { \
    fprintf(stderr, "\033[1;31m%s\033[0m\n", s); \
    exit(1); \
} while (0);


#define fatalf(fmt, args...) { \
    fprintf(stderr, "\033[0;31m"); \
    fprintf(stderr, fmt, ## args); \
    fprintf(stderr, "\033[0m\n"); \
    exit(1); \
} while (0);