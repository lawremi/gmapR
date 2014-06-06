/* $Id: dynprog_single.h 137605 2014-05-30 00:33:41Z twu $ */
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
		    bool jump_late_p, int extraband_single);

#undef T
#endif

