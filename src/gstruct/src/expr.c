static char rcsid[] = "$Id: expr.c 206180 2017-05-11 20:34:12Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "expr.h"
#include <stdio.h>		/* For sprintf */
#include <stdlib.h>
#include "mem.h"

#define T Expr_T

struct T {
  Genomicpos_T chrpos;
  long int total;
  Intlist_T acclist;
  Intlist_T exonlist;
  Intlist_T ntlist;
};


void
Expr_free (T *old) {
  Intlist_free(&(*old)->acclist);
  Intlist_free(&(*old)->exonlist);
  Intlist_free(&(*old)->ntlist);

  FREE(*old);
  return;
}


T
Expr_new (Genomicpos_T chrpos, long int total, Intlist_T acclist, Intlist_T exonlist, Intlist_T ntlist) {
  T new = (T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->total = total;
  new->acclist = acclist;
  new->exonlist = exonlist;
  new->ntlist = ntlist;

  return new;
}

Genomicpos_T
Expr_chrpos (T this) {
  return this->chrpos;
}

long int
Expr_total (T this) {
  return this->total;
}

Intlist_T
Expr_acclist (T this) {
  return this->acclist;
}

Intlist_T
Expr_exonlist (T this) {
  return this->exonlist;
}

Intlist_T
Expr_ntlist (T this) {
  return this->ntlist;
}

