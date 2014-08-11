#ifndef CIGAR_INCLUDED
#define CIGAR_INCLUDED
#include "types.h"
#include "dynprog.h"
#include "list.h"

extern char *
Dynprog_tokens_string (List_T tokens);
extern void
Dynprog_tokens_free (List_T *tokens);


extern char *
Dynprog_cigar_8 (char **md_string, int *finalc,
		 Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		 int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		 int *nindels, int queryoffset, int genomeoffset, bool revp,
		 Univcoord_T chroffset, Univcoord_T chrhigh);

extern char *
Dynprog_cigar_nogap_8 (int *nmismatches, char **md_string, int *finalc,
		       int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		       int *nindels, int queryoffset, int genomeoffset, bool revp,
		       Univcoord_T chroffset, Univcoord_T chrhigh);

extern char *
Dynprog_cigar_16 (char **md_string, int *finalc,
		  Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		  int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		  int *nindels, int queryoffset, int genomeoffset, bool revp,
		  Univcoord_T chroffset, Univcoord_T chrhigh);

extern char *
Dynprog_cigar_nogap_16 (int *nmismatches, char **md_string, int *finalc,
			int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
			int *nindels, int queryoffset, int genomeoffset, bool revp,
			Univcoord_T chroffset, Univcoord_T chrhigh);

extern char *
Dynprog_cigar_std (char **md_string, int *finalc,
		   Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
		   int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		   int *nindels, int queryoffset, int genomeoffset, bool revp,
		   Univcoord_T chroffset, Univcoord_T chrhigh);

#endif
