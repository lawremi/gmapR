/* $Id: dynprog_single.h 202212 2017-01-07 02:07:41Z twu $ */
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
Dynprog_end5_gap (int *finalscore, T dynprog, char *rev_rsequence, char *rev_gsequence,
		  int rlength, int glength);
extern int
Dynprog_end3_gap (int *finalscore, T dynprog, char *rsequence, char *gsequence,
		  int rlength, int glength);

extern int
Dynprog_single_nogap_simple (int *nmismatches, char *rsequence, char *gsequence, char *gsequence_alt, int seqlength);

extern char *
Dynprog_single_nogap (int *finalscore, int *nmismatches, int *finalc, char **md_string, T dynprog,
		      char *rsequence, char *gsequence, char *gsequence_alt,
		      int *nindels, char **delstrings, int rlength, int glength,
		      int extraband_single);

#undef T
#endif

