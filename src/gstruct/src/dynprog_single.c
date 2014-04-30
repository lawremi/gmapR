static char rcsid[] = "$Id: dynprog_single.c 134424 2014-04-25 22:23:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_single.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include <ctype.h>		/* For tolower */
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif


#include "bool.h"
#include "except.h"
#include "assert.h"
#include "mem.h"
#include "listdef.h"
#include "complement.h"
#include "dynprog_simd.h"


/* Tests whether get_genomic_nt == genomicseg in compute_scores procedures */
/* #define EXTRACT_GENOMICSEG 1 */


/* Prints parameters and results */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* print_vector */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

#define ENDSEQUENCE_PVALUE 0.001 /* Have stricter threshold for making end exons */

#define SINGLE_OPEN_HIGHQ -12	/* was -10 */
#define SINGLE_OPEN_MEDQ -8
#define SINGLE_OPEN_LOWQ -4

#define SINGLE_EXTEND_HIGHQ -3	/* was -3 */
#define SINGLE_EXTEND_MEDQ -2
#define SINGLE_EXTEND_LOWQ -1



#define T Dynprog_T


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static int
find_best_endpoint_to_queryend_indels_8 (int *bestr, int *bestc, Score8_T **matrix,
					 int rlength, int glength, int lband, int uband,
					 bool jump_late_p) {
  Score8_T bestscore = NEG_INFINITY_8;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }
  }

  return (int) bestscore;
}
#endif


static int
find_best_endpoint_to_queryend_indels_16 (int *bestr, int *bestc, Score16_T **matrix,
					  int rlength, int glength, int lband, int uband,
					  bool jump_late_p) {
  Score16_T bestscore = NEG_INFINITY_16;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }
  }

  return (int) bestscore;
}


static int
find_best_endpoint_to_queryend_indels_std (int *bestr, int *bestc, Score32_T **matrix, 
					   int rlength, int glength, int lband, int uband,
					   bool jump_late_p) {
  Score32_T bestscore = NEG_INFINITY_32;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }
  }

  return (int) bestscore;
}




#if 0
static List_T
single_gap_simple (int *finalscore, int *nmatches, int *nmismatches,
		   char *rsequence, char *rsequenceuc, int rlength, char *gsequence, char *gsequence_alt,
		   int roffset, int goffset,
		   int mismatchtype, int dynprogindex) {
  int score;
  List_T pairs = NULL;
  int r;
  int querycoord, genomecoord;
  int c1, c1_uc, c2, c2_alt;
  Pairdistance_T **pairdistance_array_type;

  debug(printf("Starting single_gap_simple\n"));
  pairdistance_array_type = pairdistance_array[mismatchtype];

  *finalscore = 0;
  *nmatches = *nmismatches = 0;

  /* Push from left to right, so we don't need to do List_reverse() later */
  for (r = 1; r <= rlength; r++) {
    querycoord = genomecoord = r-1;

    c1 = rsequence[querycoord];
    c1_uc = rsequenceuc[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u\n",goffset,genomecoord));

    } else if (c1_uc == c2 || c1_uc == c2_alt) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
    } else {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
  }

  if (*nmismatches > 1) {
    return (List_T) NULL;
  } else {
    return pairs;
  }
}
#endif


char *
Dynprog_single_gap (int *finalscore, int *finalc, char **md_string, T dynprog,
		    char *rsequence, char *gsequence, char *gsequence_alt,
		    char *nindels, char *deletion_string, int rlength, int glength,
		    bool jump_late_p, int extraband_single) {
  char *cigar;
  int bestr, bestc;

  Mismatchtype_T mismatchtype;
  int lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;

  Score16_T **matrix16;
  Direction16_T **directions16_nogap, **directions16_Egap, **directions16_Fgap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif

  mismatchtype = HIGHQ;
  open = SINGLE_OPEN_HIGHQ;
  extend = SINGLE_EXTEND_HIGHQ;

  /* Rlength: maxlookback+MAXPEELBACK.  Glength +EXTRAMATERIAL */
  debug(printf("Aligning single gap middle with wideband = %d and extraband %d\n",widebandp,extraband_single));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif
  debug(printf("\n"));

#if 0
  if (rlength > dynprog->max_rlength || glength > dynprog->max_glength) {
    debug(printf("rlength %d or glength %d is too long.  Returning NULL\n",rlength,glength));
    return NEG_INFINITY_32;
  }
#endif

  debug(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
  
  /* Have to set widebandp to be true */
  Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_single,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    matrix8 = Dynprog_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
			     rsequence,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
			     goffset,chroffset,chrhigh,watsonp,
#endif
			     mismatchtype,open,extend,
			     lband,uband,jump_late_p,/*revp*/false);

    *finalscore = find_best_endpoint_to_queryend_indels_8(&bestr,&bestc,matrix8,rlength,glength,
							  lband,uband,jump_late_p);

    cigar = Dynprog_cigar_8(&(*finalc),directions8_nogap,directions8_Egap,directions8_Fgap,
			    bestr,bestc,rsequence,gsequence,gsequence,nindels,
			    /*queryoffset*/0,/*genomeoffset*/0,
			    /*revp*/false,/*chroffset*/0,/*chrhigh*/0);
    
  } else {
    matrix16 = Dynprog_simd_16(&directions16_nogap,&directions16_Egap,&directions16_Fgap,dynprog,
			       rsequence,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
			       goffset,chroffset,chrhigh,watsonp,
#endif
			       mismatchtype,open,extend,
			       lband,uband,jump_late_p,/*revp*/false);
    *finalscore = find_best_endpoint_to_queryend_indels_16(&bestr,&bestc,matrix16,rlength,glength,
							   lband,uband,jump_late_p);

    cigar = Dynprog_cigar_16(&(*finalc),directions16_nogap,directions16_Egap,directions16_Fgap,
			     bestr,bestc,rsequence,gsequence,gsequence,nindels,
			     /*queryoffset*/0,/*genomeoffset*/0,
			     /*revp*/false,/*chroffset*/0,/*chrhigh*/0);
  }

#else

  matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			    rsequence,gsequence,gsequence_alt,rlength,glength,
			    goffset,chroffset,chrhigh,watsonp,mismatchtype,open,extend,
			    lband,uband,jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
  *finalscore = find_best_endpoint_to_queryend_indels_std(&bestr,&bestc,matrix8,rlength,glength,
							  lband,uband,jump_late_p);

  cigar = Dynprog_cigar_std(&(*finalc),directions_nogap,directions_Egap,directions_Fgap,
			    bestr,bestc,rsequence,gsequence,gsequence,nindels,
			    /*queryoffset*/0,/*genomeoffset*/0,
			    /*revp*/false,/*chroffset*/0,/*chrhigh*/0);

#endif

  /*
  Directions_free(directions);
  Matrix_free(matrix);
  */

  debug(printf("End of dynprog single gap\n"));

  return cigar;
}



