#ifndef DYNPROG_NOGAP_INCLUDED
#define DYNPROG_NOGAP_INCLUDED

#include "dynprog.h"

#define T Dynprog_T

extern Score8_T **
Dynprog_nogap_8 (T this, char *rsequence, char *gsequence, char *gsequence_alt,
		 int rlength, int glength, Mismatchtype_T mismatchtype, int uband);

extern Score16_T **
Dynprog_nogap_16 (T this, char *rsequence, char *gsequence, char *gsequence_alt,
		  int rlength, int glength, Mismatchtype_T mismatchtype, int uband);

extern int
Dynprog_nogap_find_best_endpoint_to_queryend_8 (int *bestr, int *bestc, Score8_T **matrix,
						int rlength, int glength, int lband, int uband);
extern int
Dynprog_nogap_find_best_endpoint_to_queryend_16 (int *bestr, int *bestc, Score16_T **matrix,
                                                 int rlength, int glength, int lband, int uband);

#undef T
#endif

