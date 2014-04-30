/* $Id: genomicpos.h 133073 2014-04-11 20:21:07Z twu $ */
#ifndef GENOMICPOS_INCLUDED
#define GENOMICPOS_INCLUDED
#include "types.h"

/* A genomic position, typically 3 billion or less, requiring 32 bits
   or 4 bytes */
typedef UINT4 Genomicpos_T;

/* Needed for dynprog procedures from GMAP */
typedef UINT4 Univcoord_T;
typedef UINT4 Chrpos_T;

#define T Genomicpos_T

extern char *
Genomicpos_commafmt (T N);
extern int
Genomicpos_compare (const void *a, const void *b);

#undef T
#endif
