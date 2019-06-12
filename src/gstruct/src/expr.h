/* $Id: expr.h 206180 2017-05-11 20:34:12Z twu $ */
#ifndef EXPR_INCLUDED
#define EXPR_INCLUDED
#include "intlist.h"
#include "genomicpos.h"

#define T Expr_T
typedef struct T *T;

extern void
Expr_free (T *old);
extern T
Expr_new (Genomicpos_T chrpos, long int total, Intlist_T acclist, Intlist_T exonlist, Intlist_T ntlist);

extern Genomicpos_T
Expr_chrpos (T this);
extern long int
Expr_total (T this);
extern Intlist_T
Expr_acclist (T this);
extern Intlist_T
Expr_exonlist (T this);
extern Intlist_T
Expr_ntlist (T this);


#undef T
#endif
