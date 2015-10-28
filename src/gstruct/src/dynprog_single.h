/* $Id: dynprog_single.h 143470 2014-08-05 23:16:21Z twu $ */
#ifndef DYNPROG_SINGLE_INCLUDED
#define DYNPROG_SINGLE_INCLUDED

#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "iit-read.h"
#include "types.h"
#include "dynprog.h"

#define T Dynprog_T


extern char *
Dynprog_single_gap (int *finalscore, int *finalc, char **md_string, T dynprog,
		    char *rsequence, char *gsequence, char *gsequence_alt,
		    int *nindels, char **delstrings, int rlength, int glength,
		    int extraband_single);

extern int
Dynprog_single_nogap_simple (int *nmismatches, char *rsequence, char *gsequence, char *gsequence_alt, int seqlength);

extern char *
Dynprog_single_nogap (int *finalscore, int *nmismatches, int *finalc, char **md_string, T dynprog,
		      char *rsequence, char *gsequence, char *gsequence_alt,
		      int *nindels, char **delstrings, int rlength, int glength,
		      int extraband_single);

#undef T
#endif

