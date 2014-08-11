/* $Id: translation.h 124200 2014-01-22 22:51:51Z twu $ */
#ifndef TRANSLATION_INCLUDED
#define TRANSLATION_INCLUDED
#include <stdio.h>
#include "bool.h"

#define T Translation_T
typedef struct T *T;

extern char
Translation_get_codon (char a, char b, char c);

extern char *
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 char *genome, int startpos, int endpos, int npairs,
			 bool backwardp, bool revcompp, bool fulllengthp, int cds_startpos);
#undef T
#endif




 
