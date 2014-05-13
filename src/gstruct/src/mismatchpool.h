/* $Id: mismatchpool.h 135603 2014-05-08 20:10:58Z twu $ */
#ifndef MISMATCHPOOL_INCLUDED
#define MISMATCHPOOL_INCLUDED

typedef struct Mismatchpool_T *Mismatchpool_T;

#include "mismatchdef.h"
#include "list.h"

#define T Mismatchpool_T

extern void
Mismatchpool_free (T *old);
extern void
Mismatchpool_free_memory (T this);
extern void
Mismatchpool_report_memory (T this);
extern T
Mismatchpool_new (void);
extern void
Mismatchpool_reset (T this);
extern List_T
Mismatchpool_push (List_T list, T this, char nt, int shift, int mapq, char quality, int xs);
extern List_T
Mismatchpool_pop (List_T list, Mismatch_T *x);

#undef T
#endif


