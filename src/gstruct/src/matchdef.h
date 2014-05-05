#ifndef MATCHDEF_INCLUDED
#define MATCHDEF_INCLUDED


#define T Match_T
typedef struct T *T;

struct T {
  int shift;			/* Used to record shifts */
  int mapq;
  char quality;
  long int count;
};


#undef T
#endif
