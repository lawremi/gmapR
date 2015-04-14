/* $Id: matchpool.h 135603 2014-05-08 20:10:58Z twu $ */
#ifndef MATCHPOOL_INCLUDED
#define MATCHPOOL_INCLUDED

typedef struct Matchpool_T *Matchpool_T;

#include "matchdef.h"
#include "list.h"

#define T Matchpool_T

extern void
Matchpool_free (T *old);
extern void
Matchpool_free_memory (T this);
extern void
Matchpool_report_memory (T this);
extern T
Matchpool_new (void);
extern void
Matchpool_reset (T this);
extern List_T
Matchpool_push (List_T list, T this, int shift, int mapq, char quality, int xs);
extern List_T
Matchpool_pop (List_T list, Match_T *x);

#undef T
#endif


