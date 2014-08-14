/* $Id: cappaths.h 52921 2011-11-21 21:27:25Z twu $ */
#ifndef CAPPATHS_INCLUDED
#define CAPPATHS_INCLUDED
#include "genomicpos.h"
#include "iit-read.h"

extern void
Cappaths_setup ();

extern Genomicpos_T
Cappaths_solve_genestart (Genomicpos_T pathstart,
			  long int *tally_matches, long int *tally_mismatches,
			  IIT_T end_exons_iit, char *chr, Genomicpos_T chrlength,
			  bool forwardp);
extern Genomicpos_T
Cappaths_solve_geneend (Genomicpos_T pathend,
			long int *tally_matches, long int *tally_mismatches,
			IIT_T end_exons_iit, char *chr, Genomicpos_T chrlength,
			bool forwardp);

#endif


