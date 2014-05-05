#ifndef MISMATCHDEF_INCLUDED
#define MISMATCHDEF_INCLUDED


#define T Mismatch_T
typedef struct T *T;

struct T {
  char nt;
  int shift;			/* Used to record shifts */
  char mapq;
  char quality;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  T next;	    /* Used for linking similar mismatches together */
};

#undef T
#endif

