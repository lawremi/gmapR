static char rcsid[] = "$Id: bamtally.c 178481 2015-11-09 20:55:28Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bamtally.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif

#include "except.h"
#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "complement.h"
#include "list.h"
#include "iit-read.h"
#include "iit-write.h"
#include "interval.h"
#include "genome.h"
#include "table.h"
#include "ucharlist.h"

#include "tally.h"		/* Includes calls to matchdef.h and mismatchdef.h */
#include "translation.h"


/* Specific to BAM */
#include "samflags.h"		/* For flags */
#include "samread.h"
#include "bamread.h"
#include "parserange.h"

#define ARRAY_THRESHOLD 20

static bool totals_only_p = false;


/* Alloc and block control structure */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Revise per read */
/* Can define debug1(x) as: if (linei > XX) {x;} */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Bamtally_iit */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Print block */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif



typedef enum {FIRST, SECOND} Pairend_T;


#ifdef HAVE_SSE4_1
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}

static void
print_vector_32_dec (__m128i x) {
  printf("%d %d %d %d\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}

static void
print_vector_16_dec (__m128i x) {
  printf("%d %d %d %d %d %d %d %d ",
	   _mm_extract_epi16(x,0),_mm_extract_epi16(x,1),_mm_extract_epi16(x,2),_mm_extract_epi16(x,3),
	   _mm_extract_epi16(x,4),_mm_extract_epi16(x,5),_mm_extract_epi16(x,6),_mm_extract_epi16(x,7));
  return;
}

static void
print_vector_8_dec (__m128i x) {
  printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	   _mm_extract_epi8(x,0),_mm_extract_epi8(x,1),_mm_extract_epi8(x,2),_mm_extract_epi8(x,3),
	   _mm_extract_epi8(x,4),_mm_extract_epi8(x,5),_mm_extract_epi8(x,6),_mm_extract_epi8(x,7),
	   _mm_extract_epi8(x,8),_mm_extract_epi8(x,9),_mm_extract_epi8(x,10),_mm_extract_epi8(x,11),
	   _mm_extract_epi8(x,12),_mm_extract_epi8(x,13),_mm_extract_epi8(x,14),_mm_extract_epi8(x,15));
  return;
}

#else
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 (_mm_extract_epi16(x,1) << 16) | (_mm_extract_epi16(x,0) & 0x0000FFFF),
	 (_mm_extract_epi16(x,3) << 16) | (_mm_extract_epi16(x,2) & 0x0000FFFF),
	 (_mm_extract_epi16(x,5) << 16) | (_mm_extract_epi16(x,4) & 0x0000FFFF),
	 (_mm_extract_epi16(x,7) << 16) | (_mm_extract_epi16(x,6) & 0x0000FFFF));
  return;
}

static void
print_vector_dec (__m128i x) {
  printf("%u %u %u %u\n",
	 (_mm_extract_epi16(x,1) << 16) | (_mm_extract_epi16(x,0) & 0x0000FFFF),
	 (_mm_extract_epi16(x,3) << 16) | (_mm_extract_epi16(x,2) & 0x0000FFFF),
	 (_mm_extract_epi16(x,5) << 16) | (_mm_extract_epi16(x,4) & 0x0000FFFF),
	 (_mm_extract_epi16(x,7) << 16) | (_mm_extract_epi16(x,6) & 0x0000FFFF));
  return;
}
#endif


/************************************************************************
 *   Procedures for Mismatch_T
 ************************************************************************/

static Match_T
find_match_byshift (List_T matches, int shift) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->shift == shift) {
      return match;
    }
  }
  return (Match_T) NULL;
}

static Match_T
find_match_bynm (List_T matches, int nm) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->nm == nm) {
      return match;
    }
  }
  return (Match_T) NULL;
}

static Match_T
find_match_byxs (List_T matches, int xs) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->xs == xs) {
      return match;
    }
  }
  return (Match_T) NULL;
}

/* Go -1 to -readlength, then +readlength to +1 */
static int
Match_shift_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}

static int
Match_nm_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->nm < y->nm) {
    return -1;
  } else if (y->nm < x->nm) {
    return +1;
  } else {
    return 0;
  }
}

static int
Match_xs_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->xs < y->xs) {
    return -1;
  } else if (y->xs < x->xs) {
    return +1;
  } else {
    return 0;
  }
}


/************************************************************************
 *   Procedures for Mismatch_T
 ************************************************************************/

static int
Mismatch_chain_length (Mismatch_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}


static int
Mismatch_count_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return 0;
  }
}


/* Go -1 to -readlength, then +readlength to +1 */
static int
Mismatch_shift_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}

static int
Mismatch_nm_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->nm > y->nm) {
    return -1;
  } else if (x->nm < y->nm) {
    return +1;
  } else {
    return 0;
  }
}

static int
Mismatch_xs_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->xs > y->xs) {
    return -1;
  } else if (x->xs < y->xs) {
    return +1;
  } else {
    return 0;
  }
}


static Mismatch_T
find_mismatch_byshift (List_T mismatches, char nt, int shift) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->shift == shift) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_bynm (List_T mismatches, char nt, int nm) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->nm == nm) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_byxs (List_T mismatches, char nt, int xs) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->xs == xs) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_nt (List_T mismatches, char nt) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}


static int
Insertion_chain_length (Insertion_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}

/* Go -1 to -readlength, then +readlength to +1 */
static int
Insertion_shift_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}

static int
Insertion_nm_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->nm < y->nm) {
    return -1;
  } else if (y->nm < x->nm) {
    return +1;
  } else {
    return 0;
  }
}

static int
Insertion_xs_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->xs < y->xs) {
    return -1;
  } else if (y->xs < x->xs) {
    return +1;
  } else {
    return 0;
  }
}


static int
Deletion_chain_length (Deletion_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}

/* Go -1 to -readlength, then +readlength to +1 */
static int
Deletion_shift_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}

static int
Deletion_nm_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->nm < y->nm) {
    return -1;
  } else if (y->nm < x->nm) {
    return +1;
  } else {
    return 0;
  }
}

static int
Deletion_xs_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->xs < y->xs) {
    return -1;
  } else if (y->xs < x->xs) {
    return +1;
  } else {
    return 0;
  }
}


/************************************************************************
 *   Probability calculations
 ************************************************************************/

#define NGENOTYPES 10

static char *genotype_string[NGENOTYPES] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
static bool genotype_heterozygousp[NGENOTYPES] = {false,true,true,true,false,true,true,false,true,false};


static double
allele_logprob[3][MAX_QUALITY_SCORE+1] = {
  /* Not in genotype: log(1/3*10^(-Q/10)) */
  {-1.098612,
   -1.328871, -1.559129, -1.789388, -2.019646, -2.249905,
   -2.480163, -2.710422, -2.940680, -3.170939, -3.401197,
   -3.631456, -3.861714, -4.091973, -4.322231, -4.552490,
   -4.782748, -5.013007, -5.243265, -5.473524, -5.703782,
   -5.934041, -6.164299, -6.394558, -6.624817, -6.855075,
   -7.085334, -7.315592, -7.545851, -7.776109, -8.006368,
   -8.236626, -8.466885, -8.697143, -8.927402, -9.157660,
   -9.387919, -9.618177, -9.848436, -10.078694, -10.308953},

  /* In heterozygote: log(1/2 - 1/3*10^(-Q/10)) */
  {-1.7917595,
   -1.4472174, -1.2389754, -1.0998002, -1.0015828, -0.9299061,
   -0.8764201, -0.8358837, -0.8048159, -0.7808079, -0.7621401,
   -0.7475561, -0.7361213, -0.7271306, -0.7200462, -0.7144544,
   -0.7100349, -0.7065382, -0.7037694, -0.7015754, -0.6998362,
   -0.6984568, -0.6973624, -0.6964940, -0.6958048, -0.6952576,
   -0.6948232, -0.6944782, -0.6942043, -0.6939868, -0.6938141,
   -0.6936769, -0.6935679, -0.6934814, -0.6934126, -0.6933580,
   -0.6933147, -0.6932802, -0.6932528, -0.6932311, -0.6932138},

  /* In homozygote: log(1 - 10^(-Q/10)) */
  {/*-Inf*/-1.58,
   -1.5814737534, -0.9968430440, -0.6955244713, -0.5076758737, -0.3801304081,
   -0.2892681872, -0.2225515160, -0.1725565729, -0.1345519603, -0.1053605157,
   -0.0827653027, -0.0651741732, -0.0514182742, -0.0406248442, -0.0321335740,
   -0.0254397275, -0.0201543648, -0.0159758692, -0.0126691702, -0.0100503359,
   -0.0079749983, -0.0063295629, -0.0050244739, -0.0039890173, -0.0031672882,
   -0.0025150465, -0.0019972555, -0.0015861505, -0.0012597185, -0.0010005003,
   -0.0007946439, -0.0006311565, -0.0005013129, -0.0003981864, -0.0003162778,
   -0.0002512202, -0.0001995461, -0.0001585019, -0.0001259005, -0.0001000050}
};


#if 0
static void
make_quality_scores_constant (int quality_score_constant) {
  int i, Q;
  double value;

  for (i = 0; i < 3; i++) {
    value = allele_logprob[i][quality_score_constant];
    for (Q = 0; Q <= MAX_QUALITY_SCORE; Q++) {
      allele_logprob[i][Q] = value;
    }
  }

  return;
}
#endif


static int
genotype_compatibility[4][NGENOTYPES] = {
  
  {/*AA*/ 2, /*AC*/ 1, /*AG*/ 1, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 0,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* C */
  {/*AA*/ 0, /*AC*/ 1, /*AG*/ 0, /*AT*/ 0,
   /*CC*/ 2, /*CG*/ 1, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* G */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 1, /*AT*/ 0,
   /*CC*/ 0, /*CG*/ 1, /*CT*/ 0,
   /*GG*/ 2, /*GT*/ 1,
   /*TT*/ 0},

  /* T */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 0, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 1,
   /*TT*/ 2}
};


#if 0
static int
compute_probs_list (double *probs, double *loglik, char refnt, int minimum_quality_score,
		    int quality_score_adj) {
  int bestg;
  double total, adj_loglik[NGENOTYPES], maxlik;
  List_T p;
  Match_T match;
  Mismatch_T mismatch;
  int Q;
  int g;

  for (g = 0; g < NGENOTYPES; g++) {
    loglik[g] = 0.0;
  }

  for (p = list_matches_byquality; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    Q = minimum_quality_score - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (refnt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }

  for (p = mismatches_byquality; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    Q = mismatch->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (mismatch->nt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }
  
  /* Raise all loglikelihoods to maximum, to avoid underflow when taking exp() */
  maxlik = loglik[0];
  bestg = 0;
  for (g = 1; g < NGENOTYPES; g++) {
    if (loglik[g] > maxlik) {
      maxlik = loglik[g];
      bestg = g;
    } else if (loglik[g] == maxlik) {
      bestg = -1;
    }
  }

  for (g = 0; g < NGENOTYPES; g++) {
    adj_loglik[g] = loglik[g] - maxlik;
  }

  total = 0.0;
  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] = exp(adj_loglik[g]);
    total += probs[g];
  }

  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] /= total;
  }


  return bestg;
}



static int
compute_probs_array (double *probs, double *loglik, char refnt,
		     int *matches_byquality, int n_matches_byquality,
		     List_T mismatches_byquality, int quality_score_adj) {
  int bestg;
  int quality;
  int count;
  double total, adj_loglik[NGENOTYPES], maxlik;
  List_T p;
  Mismatch_T mismatch;
  int Q;
  int g;

  for (g = 0; g < NGENOTYPES; g++) {
    loglik[g] = 0.0;
  }

  for (quality = 0; quality < n_matches_byquality; quality++) {
    if ((count = matches_byquality[quality]) > 0) {
      Q = quality - quality_score_adj;
      if (Q < 1) {
	Q = 1;
      } else if (Q > MAX_QUALITY_SCORE) {
	Q = MAX_QUALITY_SCORE;
      }

      switch (refnt) {
      case 'A':
	for (g = 0; g < NGENOTYPES; g++) {
	  loglik[g] += count * allele_logprob[genotype_compatibility[0][g]][Q];
	}
	break;

      case 'C':
	for (g = 0; g < NGENOTYPES; g++) {
	  loglik[g] += count * allele_logprob[genotype_compatibility[1][g]][Q];
	}
	break;

      case 'G':
	for (g = 0; g < NGENOTYPES; g++) {
	  loglik[g] += count * allele_logprob[genotype_compatibility[2][g]][Q];
	}
	break;

      case 'T':
	for (g = 0; g < NGENOTYPES; g++) {
	  loglik[g] += count * allele_logprob[genotype_compatibility[3][g]][Q];
	}
	break;
      }
    }
  }

  for (p = mismatches_byquality; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    Q = mismatch->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (mismatch->nt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }
  
  /* Raise all loglikelihoods to maximum, to avoid underflow when taking exp() */
  maxlik = loglik[0];
  bestg = 0;
  for (g = 1; g < NGENOTYPES; g++) {
    if (loglik[g] > maxlik) {
      maxlik = loglik[g];
      bestg = g;
    } else if (loglik[g] == maxlik) {
      bestg = -1;
    }
  }

  for (g = 0; g < NGENOTYPES; g++) {
    adj_loglik[g] = loglik[g] - maxlik;
  }

  total = 0.0;
  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] = exp(adj_loglik[g]);
    total += probs[g];
  }

  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] /= total;
  }

  return bestg;
}
#endif




static bool
pass_variant_filter_p (long int nmatches, long int delcounts_plus, long int delcounts_minus,
		       List_T mismatches_by_shift, int min_depth, int variant_strands) {
  long int depth;
  bool plus_strand_p[4], minus_strand_p[4];
  List_T ptr;
  Mismatch_T mismatch;

  depth = nmatches + delcounts_plus + delcounts_minus;
  ptr = mismatches_by_shift;
  while (ptr != NULL) {
    mismatch = (Mismatch_T) List_head(ptr);
    depth += mismatch->count;
    ptr = List_next(ptr);
  }
  if (depth < min_depth) {
    return false;
  }


  if (variant_strands == 0) {
    return true;
  } else if (variant_strands == 1) {
    /* Not sure if we should look at insertions_byshift also */
    if (mismatches_by_shift == NULL && delcounts_plus == 0 && delcounts_minus == 0) {
      return false;
    } else {
      return true;
    }

  } else {
    if (delcounts_plus > 0 && delcounts_minus > 0) {
      return true;
    }

    plus_strand_p[0] = plus_strand_p[1] = plus_strand_p[2] = plus_strand_p[3] = false;
    minus_strand_p[0] = minus_strand_p[1] = minus_strand_p[2] = minus_strand_p[3] = false;

    /* Look for two strands on a given mismatch allele */
    for (ptr = mismatches_by_shift; ptr != NULL; ptr = List_next(ptr)) {
      mismatch = (Mismatch_T) List_head(ptr);
      if (mismatch->shift > 0) {
	if (delcounts_minus > 0) {
	  return true;
	}
	switch (mismatch->nt) {
	case 'A':
	  if (minus_strand_p[/*A*/0] == true) {
	    return true;
	  } else {
	    plus_strand_p[/*A*/0] = true;
	  }
	  break;
	case 'C':
	  if (minus_strand_p[/*C*/1] == true) {
	    return true;
	  } else {
	    plus_strand_p[/*C*/1] = true;
	  }
	  break;
	case 'G': 
	  if (minus_strand_p[/*G*/2] == true) {
	    return true;
	  } else {
	    plus_strand_p[/*G*/2] = true;
	  }
	  break;
	case 'T':
	  if (minus_strand_p[/*T*/3] == true) {
	    return true;
	  } else {
	    plus_strand_p[/*T*/3] = true;
	  }
	  break;
	}

      } else if (mismatch->shift < 0) {
	if (delcounts_plus > 0) {
	  return true;
	}
	switch (mismatch->nt) {
	case 'A':
	  if (plus_strand_p[/*A*/0] == true) {
	    return true;
	  } else {
	    minus_strand_p[/*A*/0] = true;
	  }
	  break;
	case 'C':
	  if (plus_strand_p[/*C*/1] == true) {
	    return true;
	  } else {
	    minus_strand_p[/*C*/1] = true;
	  }
	  break;
	case 'G':
	  if (plus_strand_p[/*G*/2] == true) {
	    return true;
	  } else {
	    minus_strand_p[/*G*/2] = true;
	  }
	  break;
	case 'T':
	  if (plus_strand_p[/*T*/3] == true) {
	    return true;
	  } else {
	    minus_strand_p[/*T*/3] = true;
	  }
	  break;
	}
      }
    }

    return false;
  }
}


static void
print_allele_counts_simple (FILE *fp, Tally_T this, Genome_T genome, Genomicpos_T chroffset,
			    Genomicpos_T chrpos) {
  List_T p;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;

  if (this->insertions_byshift != NULL) {
    fprintf(fp,"\n");
    fprintf(fp,"^");
    for (p = this->insertions_byshift; p != NULL; p = List_next(p)) {
      ins = (Insertion_T) List_head(p);
      fprintf(fp,"%s %ld ref:%ld",ins->segment,ins->count,this->n_fromleft_plus+this->n_fromleft_minus);
    }
    fprintf(fp,"\n");
  }

  if (this->deletions_byshift != NULL) {
    fprintf(fp,"\n");
    fprintf(fp,"_");
    for (p = this->deletions_byshift; p != NULL; p = List_next(p)) {
      del = (Deletion_T) List_head(p);
      fprintf(fp,"%s %ld ref:%ld",del->segment,del->count,this->n_fromleft_plus+this->n_fromleft_minus);
    }
    fprintf(fp,"\n");
  }

  if (this->nmatches + this->delcounts_plus + this->delcounts_minus == 0) {
    fprintf(fp,"%c0",Genome_get_char(genome,chroffset+chrpos-1U));
  } else {
    fprintf(fp,"%c%d",this->refnt,this->nmatches);
  }

  for (p = this->mismatches_byshift; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    fprintf(fp," %c%ld",mismatch->nt,mismatch->count);
  }

  if (this->delcounts_plus + this->delcounts_minus > 0) {
    fprintf(fp," _%ld",this->delcounts_plus + this->delcounts_minus);
  }

  return;
}


static bool
pass_difference_filter_p (double *probs, double *loglik, Tally_T this,
			  Genome_T genome, char *printchr, Genomicpos_T chroffset, Genomicpos_T chrpos,
			  int quality_score_adj, bool genomic_diff_p) {
  int bestg;

  if (genomic_diff_p == false) {
    return true;
  } else {
#if 0
    if (this->use_array_p == false) {
      bestg = compute_probs_list(probs,loglik,this->refnt,this->nmatches,
				 this->mismatches_byquality,quality_score_adj);
    } else {
      bestg = compute_probs_array(probs,loglik,this->refnt,this->matches_byquality,this->n_matches_byquality,
				  this->mismatches_byquality,quality_score_adj);
    }
#endif

    if (0 && bestg < 0) {
      return false;
    } else if (genotype_heterozygousp[bestg] == true) {
      if (this->refnt == 'A') {
	if (genotype_compatibility[0][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  printchr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'C') {
	if (genotype_compatibility[1][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  printchr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'G') {
	if (genotype_compatibility[2][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  printchr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'T') {
	if (genotype_compatibility[3][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  printchr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else {
	fprintf(stderr,"Reference nt not ACGT\n");
	return false;
      }

    } else if (this->refnt == 'A') {
      return (bestg == 0) ? false : true;
    } else if (this->refnt == 'C') {
      return (bestg == 4) ? false : true;
    } else if (this->refnt == 'G') {
      return (bestg == 7) ? false : true;
    } else if (this->refnt == 'T') {
      return (bestg == 9) ? false : true;
    } else {
      fprintf(stderr,"Reference nt not ACGT\n");
      return false;
    }
  }
}



static Genomicpos_T
print_runlength (Tally_T *block_tallies, Genomicpos_T *exonstart, Genomicpos_T lastpos,
		 Genomicpos_T blockstart, Genomicpos_T blockptr, 
		 char *printchr) {
  Tally_T this;
  int blocki, lasti;
  Genomicpos_T chrpos;

  lasti = blockptr - blockstart;

  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (this->nmatches > 0 || this->delcounts_plus + this->delcounts_minus > 0 ||
	this->mismatches_byshift != NULL || this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
      chrpos = blockstart + blocki;
      if (chrpos > lastpos + 1U) {
	if (*exonstart == 0U) {
	  /* Skip printing */
	} else {
	  printf(">1 %s:%u..%u\n",printchr,*exonstart,lastpos);
	}
	*exonstart = chrpos;
      }
      lastpos = chrpos;
    }
  }

  return lastpos;
}


static List_T
make_mismatches_unique (List_T mismatches
#ifdef USE_MISMATCHPOOL
			, Mismatchpool_T mismatchpool
#endif
			) {
  List_T unique_mismatches = NULL, ptr;
  Mismatch_T mismatch, mismatch0;

  for (ptr = mismatches; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    if ((mismatch0 = find_mismatch_nt(unique_mismatches,mismatch->nt)) != NULL) {
      mismatch0->count += mismatch->count;

      /* Insert mismatch into list */
      mismatch->next = mismatch0->next;
      mismatch0->next = mismatch;

      mismatch0->shift += 1; /* Used here as nshifts */
    } else {
#ifdef USE_MISMATCHPOOL
      unique_mismatches = Mismatchpool_push(unique_mismatches,mismatchpool,
					    mismatch->nt,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
      mismatch0 = (Mismatch_T) List_head(unique_mismatches);
      mismatch0->count = mismatch->count;
      mismatch0->next = mismatch;
#else
      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
      mismatch0->count = mismatch->count;
      mismatch0->next = mismatch;
      unique_mismatches = List_push(unique_mismatches,mismatch0);
#endif
    }
  }

  return unique_mismatches;
}

static List_T
make_mismatches_unique_signed (List_T mismatches
#ifdef USE_MISMATCHPOOL
	, Mismatchpool_T mismatchpool
#endif
	) {
  List_T unique_mismatches = NULL, ptr;
  Mismatch_T mismatch, mismatch0;

  for (ptr = mismatches; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    if ((mismatch0 = find_mismatch_nt(unique_mismatches,mismatch->nt)) != NULL) {
      mismatch0->count += mismatch->count;
      if (mismatch->shift > 0) {
	mismatch0->count_plus += mismatch->count;
      } else {
	mismatch0->count_minus += mismatch->count;
      }

      /* Insert mismatch into list */
      mismatch->next = mismatch0->next;
      mismatch0->next = mismatch;

      mismatch0->shift += 1; /* Used here as nshifts */
    } else {
#ifdef USE_MISMATCHPOOL
      unique_mismatches = Mismatchpool_push(unique_mismatches,mismatchpool,
					    mismatch->nt,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
      mismatch0 = (Mismatch_T) List_head(unique_mismatches);
      mismatch0->count = mismatch->count;
      if (mismatch->shift > 0) {
	mismatch0->count_plus = mismatch->count;
	mismatch0->count_minus = 0;
      } else {
	mismatch0->count_minus = mismatch->count;
	mismatch0->count_plus = 0;
      }
      mismatch0->next = mismatch;
#else
      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
      mismatch0->count = mismatch->count;
      if (mismatch->shift > 0) {
	mismatch0->count_plus = mismatch->count;
	mismatch0->count_minus = 0;
      } else {
	mismatch0->count_minus = mismatch->count;
	mismatch0->count_plus = 0;
      }
      mismatch0->next = mismatch;
      unique_mismatches = List_push(unique_mismatches,mismatch0);
#endif
    }
  }

  return unique_mismatches;
}



static long int
block_total (Tally_T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr) {
  long int total;
  Tally_T this;
  int blocki, lasti;
  List_T ptr;
  Mismatch_T mismatch;

  lasti = blockptr - blockstart;

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    total += this->nmatches + this->delcounts_plus + this->delcounts_minus;
    for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
      mismatch = (Mismatch_T) List_head(ptr);
      total += mismatch->count;
    }
    Tally_clear(this);
  }

  return total;
}


static void
print_zeroes (Genomicpos_T start, Genomicpos_T end, char *printchr, int blocksize, 
	      Genome_T genome, Genomicpos_T chroffset, bool blockp) {
  Genomicpos_T chrpos, chrpos0, blockstart, blockend;
  
  if (start < end) {
    for (chrpos = start; chrpos + blocksize < end; chrpos += blocksize) {
      blockstart = chrpos;
      blockend = chrpos + blocksize;
      if (blockp == true) {
	printf(">%ld %s:%u..%u\n",/*total*/0,printchr,blockstart,blockend-1U);
      }
      for (chrpos0 = blockstart; chrpos0 < blockend; chrpos0++) {
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos0);
	}
	printf("%c0\n",Genome_get_char(genome,chroffset+chrpos0-1U));
      }
    }

    if (chrpos < end) {
      blockstart = chrpos;
      blockend = end;
      if (blockp == true) {
	printf(">%ld %s:%u..%u\n",/*total*/0,printchr,blockstart,blockend-1U);
      }
      for (chrpos0 = blockstart; chrpos0 < blockend; chrpos0++) {
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos0);
	}
	printf("%c0\n",Genome_get_char(genome,chroffset+chrpos0-1U));
      }
    }
  }

  return;
}


/************************************************************************
 *  Gene coordinates
 ************************************************************************/


static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_inplace (char *sequence, Genomicpos_T length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


typedef struct Exon_T *Exon_T;

struct Exon_T {
  Genomicpos_T exonstart;
  Genomicpos_T exonend;
  int length;
  char *ntstring;
};


static void
Exon_free (Exon_T *old) {
  FREE((*old)->ntstring);
  FREE(*old);
  return;
}

static Exon_T
Exon_new (Genomicpos_T exonstart, Genomicpos_T exonend, char *chr, Genomicpos_T chroffset,
	  int sign, Genome_T genome) {
  Exon_T new = (Exon_T) MALLOC(sizeof(*new));
  Genomicpos_T exonlow;

  new->exonstart = exonstart;
  new->exonend = exonend;

  if (sign > 0) {
    exonlow = exonstart;
    new->length = exonend - exonstart + 1;
    new->ntstring = (char *) CALLOC(new->length + 1,sizeof(char));
    Genome_fill_buffer_simple(genome,chroffset+exonlow-1U,new->length,new->ntstring);
  } else {
    exonlow = exonend;
    new->length = exonstart - exonend + 1;
    new->ntstring = (char *) CALLOC(new->length + 1,sizeof(char));
    Genome_fill_buffer_simple(genome,chroffset+exonlow-1U,new->length,new->ntstring);
    make_complement_inplace(new->ntstring,new->length);
  }

  return new;
}


static Exon_T *
get_exons (int *nexons, char *annot, char *chr, Genomicpos_T chroffset, int sign, Genome_T genome) {
  Exon_T *array;
  List_T exons;
  Genomicpos_T exonstart, exonend;
  char *p;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;

  exons = (List_T) NULL;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      exons = List_push(exons,(void *) Exon_new(exonstart,exonend,chr,chroffset,sign,genome));
    }

    /* Advance to next exon */
    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  exons = List_reverse(exons);
  *nexons = List_length(exons);
  array = (Exon_T *) List_to_array(exons,NULL);
  List_free(&exons);

  return array;
}


typedef struct Gene_T *Gene_T;
struct Gene_T {
  int exoni;

  char *acc;
  char *genename;
  Exon_T *exons;
  int *cum_exonlength;
  int nexons;

  int translation_start;
  int translation_end;
};


static void
Gene_free (Gene_T *old) {
  int i;

  FREE((*old)->acc);
  FREE((*old)->genename);
  FREE((*old)->cum_exonlength);
  for (i = 0; i < (*old)->nexons; i++) {
    Exon_free(&((*old)->exons[i]));
  }
  FREE((*old)->exons);

  FREE(*old);
  return;
}

static Gene_T
Gene_new (char *acc, char *genename, Exon_T *exons, int nexons, int sign, Genomicpos_T startcoord) {
  Gene_T new = (Gene_T) MALLOC(sizeof(*new));
  int exoni;

  new->acc = acc;
  new->genename = genename;
  new->exons = exons;
  new->nexons = nexons;
  new->translation_start = 0;
  new->translation_end = -1;

  new->cum_exonlength = (int *) MALLOC(nexons * sizeof(int));
  new->cum_exonlength[0] = 0;
  if (sign > 0) {
    for (exoni = 1; exoni < nexons; exoni++) {
      new->cum_exonlength[exoni] = new->cum_exonlength[exoni-1] + (exons[exoni-1]->exonend - exons[exoni-1]->exonstart + 1);
    }
    new->exoni = 0;

  } else {
    for (exoni = 1; exoni < nexons; exoni++) {
      new->cum_exonlength[exoni] = new->cum_exonlength[exoni-1] + (exons[exoni-1]->exonstart - exons[exoni-1]->exonend + 1);
    }
    new->exoni = nexons - 1;

  }

  return new;
}


typedef struct Chrpos_pair_T *Chrpos_pair_T;
struct Chrpos_pair_T {
  Genomicpos_T frame0pos;
  Genomicpos_T frame1pos;
};

static void
Chrpos_pair_free (Chrpos_pair_T *old) {
  FREE(*old);
  return;
}

static Chrpos_pair_T
Chrpos_pair_new (Genomicpos_T frame0pos, Genomicpos_T frame1pos) {
  Chrpos_pair_T new = (Chrpos_pair_T) MALLOC(sizeof(*new));

  new->frame0pos = frame0pos;
  new->frame1pos = frame1pos;
  return new;
}

static void
Chrpos_pair_print (Chrpos_pair_T this) {
  printf(" %u,%u",this->frame0pos,this->frame1pos);
  return;
}

static int
Chrpos_pair_cmp (const void *a, const void *b) {
  Chrpos_pair_T x = * (Chrpos_pair_T *) a;
  Chrpos_pair_T y = * (Chrpos_pair_T *) b;

  if (x->frame0pos < y->frame0pos) {
    return -1;
  } else if (y->frame0pos < x->frame0pos) {
    return +1;
  } else if (x->frame1pos < y->frame1pos) {
    return -1;
  } else if (y->frame1pos < x->frame1pos) {
    return +1;
  } else {
    return 0;
  }
}

static List_T
Chrpos_pair_uniq (List_T pairs) {
  List_T unique = NULL;
  Chrpos_pair_T *array;
  int npairs, i, j, k;

  if (pairs == NULL) {
    return (List_T) NULL;

  } else {
    array = (Chrpos_pair_T *) List_to_array_n(&npairs,pairs);
    qsort(array,npairs,sizeof(Chrpos_pair_T),Chrpos_pair_cmp);
    List_free(&pairs);

    i = 0;
    while (i < npairs) {
      unique = List_push(unique,(void *) array[i]);
      j = i + 1;
      while (j < npairs && !Chrpos_pair_cmp(&(array[j]),&array[i])) {
	j++;
      }
      for (k = i + 1; k < j; k++) {
	Chrpos_pair_free(&(array[k]));
      }

      i = j;
    }

    FREE(array);
    return unique;
  }
}
    

static char *
concatenate_exons (int *length, Exon_T *exons, int nexons) {
  char *gene_sequence;
  int i;

  *length = 0;
  for (i = 0; i < nexons; i++) {
    *length += exons[i]->length;
  }

  gene_sequence = MALLOC((*length + 1) * sizeof(char));
  *length = 0;
  for (i = 0; i < nexons; i++) {
    strcpy(&(gene_sequence[*length]),exons[i]->ntstring);
    *length += exons[i]->length;
  }
  gene_sequence[*length] = '\0';
  return gene_sequence;
}


static char *codon_table[64] = 
  {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
   "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
   "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
   "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

static char *aa_abbrev_table[64] = 
  {"Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", "Ile", "Ile", "Met", "Ile",
   "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro", "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu",
   "Glu", "Asp", "Glu", "Asp", "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
   "***", "Tyr", "***", "Tyr", "Ser", "Ser", "Ser", "Ser", "***", "Cys", "Trp", "Cys", "Leu", "Phe", "Leu", "Phe"};

static char aa_table[65] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";


static Ucharlist_T
push_char (int *nbytes, Ucharlist_T list, unsigned char x) {
#ifdef DEBUG2
  if (x == 0) {
    printf("%d: Pushing terminating char 0\n",*nbytes);
  } else {
    printf("%d: Pushing char %c\n",*nbytes,x);
  }
#endif
  list = Ucharlist_push(list,x);
  *nbytes += 1;
  return list;
}

static Ucharlist_T
push_codon (int *nbytes, Ucharlist_T list, unsigned char x) {
#ifdef DEBUG2
  if (x == 255) {
    printf("%d: Pushing codon terminator %d\n",*nbytes,x);
  } else {
    printf("%d: Pushing codon %d (%s)\n",*nbytes,x,codon_table[x]);
  }
#endif
  list = Ucharlist_push(list,x);
  *nbytes += 1;
  return list;
}

static Ucharlist_T
push_string (int *nbytes, Ucharlist_T list, char *string) {
  int length, i;

  debug2(printf("%d: Pushing string %s\n",*nbytes,string));
  length = strlen(string);
  for (i = 0; i < length; i++) {
    list = Ucharlist_push(list,string[i]);
  }
  list = Ucharlist_push(list,'\0');
  *nbytes += (length + 1);
  return list;
}

static Ucharlist_T
push_int (int *nbytes, Ucharlist_T list, int value) {
  unsigned char x;

  debug2(printf("%d: Pushing int %d\n",*nbytes,value));
  x = value & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);
  
  *nbytes += 4;
  return list;
}



typedef struct Cell_T *Cell_T;
struct Cell_T {
  int value;
  int tally;
};

static void
Cell_free (Cell_T *old) {
  FREE(*old);
  return;
}

static Cell_T
Cell_new (int value, int tally) {
  Cell_T new = (Cell_T) MALLOC(sizeof(*new));

  new->value = value;
  new->tally = tally;
  return new;
}


static void
print_shift_list (Intlist_T shift_list) {
  int *shift_tally;
  int shift, max_shift = -1000000, min_shift = 1000000;
  Intlist_T p;
  bool firstp = true;

  if (shift_list != NULL) {
    for (p = shift_list; p != NULL; p = Intlist_next(p)) {
      shift = Intlist_head(p);
      if (shift > max_shift) {
	max_shift = shift;
      }
      if (shift < min_shift) {
	min_shift = shift;
      }
    }

    shift_tally = (int *) CALLOC(max_shift - min_shift + 1,sizeof(int));

    for (p = shift_list; p != NULL; p = Intlist_next(p)) {
      shift = Intlist_head(p);
      shift_tally[shift - min_shift] += 1;
    }
  
    printf("(");

    if (max_shift < 0) {
      shift = max_shift;
    } else {
      shift = -1;
    }
    for ( ; shift >= min_shift; shift--) {
      if (shift_tally[shift - min_shift] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%d@%d",shift_tally[shift - min_shift],shift);
      }
    }

    for (shift = max_shift; shift >= 0 && shift >= min_shift; shift--) {
      if (shift_tally[shift - min_shift] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%d@%d",shift_tally[shift - min_shift],shift);
      }
    }
    printf(")");

    FREE(shift_tally);
  }

  return;
}

static Ucharlist_T
iit_shift_list (Ucharlist_T bytes, int *nbytes, Intlist_T shift_list) {
  int *shift_tally;
  int shift, max_shift = -1000000, min_shift = 1000000;
  int length;
  Intlist_T p;

  if (shift_list != NULL) {
    for (p = shift_list; p != NULL; p = Intlist_next(p)) {
      shift = Intlist_head(p);
      if (shift > max_shift) {
	max_shift = shift;
      }
      if (shift < min_shift) {
	min_shift = shift;
      }
    }

    shift_tally = (int *) CALLOC(max_shift - min_shift + 1,sizeof(int));

    for (p = shift_list; p != NULL; p = Intlist_next(p)) {
      shift = Intlist_head(p);
      shift_tally[shift - min_shift] += 1;
    }
  
    length = 0;
    for (shift = min_shift; shift <= max_shift; shift++) {
      if (shift_tally[shift - min_shift] > 0) {
	length++;
      }
    }
    bytes = push_int(&(*nbytes),bytes,length);

    if (max_shift < 0) {
      shift = max_shift;
    } else {
      shift = -1;
    }
    for ( ; shift >= min_shift; shift--) {
      if (shift_tally[shift - min_shift] > 0) {
	bytes = push_int(&(*nbytes),bytes,shift);
	bytes = push_int(&(*nbytes),bytes,shift_tally[shift - min_shift]);
      }
    }

    for (shift = max_shift; shift >= 0 && shift >= min_shift; shift--) {
      if (shift_tally[shift - min_shift] > 0) {
	bytes = push_int(&(*nbytes),bytes,shift);
	bytes = push_int(&(*nbytes),bytes,shift_tally[shift - min_shift]);
      }
    }

    FREE(shift_tally);
  }

  return bytes;
}


static void
print_nm_list (Intlist_T nm_list) {
  int *nm_tally;
  int nm, max_nm = 0, min_nm = 100;
  Intlist_T p;
  bool firstp = true;

  if (nm_list != NULL) {
    for (p = nm_list; p != NULL; p = Intlist_next(p)) {
      nm = Intlist_head(p);
      if (nm > max_nm) {
	max_nm = nm;
      }
      if (nm < min_nm) {
	min_nm = nm;
      }
    }

    nm_tally = (int *) CALLOC(max_nm - min_nm + 1,sizeof(int));

    for (p = nm_list; p != NULL; p = Intlist_next(p)) {
      nm = Intlist_head(p);
      nm_tally[nm - min_nm] += 1;
    }
  
    printf("(");

    for (nm = min_nm; nm <= max_nm; nm++) {
      if (nm_tally[nm - min_nm] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%dNM%d",nm_tally[nm - min_nm],nm);
      }
    }
    printf(")");

    FREE(nm_tally);
  }

  return;
}

static Ucharlist_T
iit_nm_list (Ucharlist_T bytes, int *nbytes, Intlist_T nm_list) {
  int *nm_tally;
  int nm, max_nm = 0, min_nm = 100;
  int length;
  Intlist_T p;

  if (nm_list != NULL) {
    for (p = nm_list; p != NULL; p = Intlist_next(p)) {
      nm = Intlist_head(p);
      if (nm > max_nm) {
	max_nm = nm;
      }
      if (nm < min_nm) {
	min_nm = nm;
      }
    }

    nm_tally = (int *) CALLOC(max_nm - min_nm + 1,sizeof(int));

    for (p = nm_list; p != NULL; p = Intlist_next(p)) {
      nm = Intlist_head(p);
      nm_tally[nm - min_nm] += 1;
    }
  
    length = 0;
    for (nm = min_nm; nm <= max_nm; nm++) {
      if (nm_tally[nm - min_nm] > 0) {
	length++;
      }
    }
    bytes = push_int(&(*nbytes),bytes,length);

    for (nm = min_nm; nm <= max_nm; nm++) {
      if (nm_tally[nm - min_nm] > 0) {
	bytes = push_int(&(*nbytes),bytes,nm);
	bytes = push_int(&(*nbytes),bytes,nm_tally[nm - min_nm]);
      }
    }

    FREE(nm_tally);
  }

  return bytes;
}


#if 0
static void
print_quality_list (Intlist_T quality_list) {
  int quality_tally[MAX_QUALITY_SCORE+1];
  int quality, max_quality = 0, min_quality = MAX_QUALITY_SCORE;
  Intlist_T p;
  bool firstp = true;

  if (quality_list != NULL) {
    memset(quality_tally,0,(MAX_QUALITY_SCORE+1)*sizeof(int));
    for (p = quality_list; p != NULL; p = Intlist_next(p)) {
      quality = Intlist_head(p);
      if (quality > max_quality) {
	max_quality = quality;
      }
      if (quality < min_quality) {
	min_quality = quality;
      }
      quality_tally[quality] += 1;
    }
  
    printf("(");
    for (quality = max_quality; quality >= min_quality; quality--) {
      if (quality_tally[quality] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%dQ%d",quality_tally[quality],quality);
      }
    }
    printf(")");
  }

  return;
}
#endif

#if 0
static Ucharlist_T
iit_quality_list (Ucharlist_T bytes, int *nbytes, Intlist_T quality_list, int quality_score_adj) {
  int length;
  int quality_tally[MAX_QUALITY_SCORE+1];
  int quality, max_quality = 0, min_quality = MAX_QUALITY_SCORE;
  Intlist_T p;

  if (quality_list != NULL) {
    memset(quality_tally,0,(MAX_QUALITY_SCORE+1)*sizeof(int));
    for (p = quality_list; p != NULL; p = Intlist_next(p)) {
      quality = Intlist_head(p);
      if (quality > max_quality) {
	max_quality = quality;
      }
      if (quality < min_quality) {
	min_quality = quality;
      }
      quality_tally[quality] += 1;
    }
  
    length = 0;
    for (quality = max_quality; quality >= min_quality; quality--) {
      if (quality_tally[quality] > 0) {
	length += 1;
      }
    }
    bytes = push_int(&(*nbytes),bytes,length);


    for (quality = max_quality; quality >= min_quality; quality--) {
      if (quality_tally[quality] > 0) {
	bytes = push_int(&(*nbytes),bytes,quality - quality_score_adj);
	bytes = push_int(&(*nbytes),bytes,quality_tally[quality]);
      }
    }
  }

  return bytes;
}
#endif


#if 0
static void
print_mapq_list (Intlist_T mapq_list) {
  int mapq_tally[MAX_MAPQ_SCORE+1];
  int mapq, max_mapq = 0, min_mapq = MAX_MAPQ_SCORE;
  Intlist_T p;
  bool firstp = true;

  if (mapq_list != NULL) {
    memset(mapq_tally,0,(MAX_MAPQ_SCORE+1)*sizeof(int));
    for (p = mapq_list; p != NULL; p = Intlist_next(p)) {
      mapq = Intlist_head(p);
      if (mapq > max_mapq) {
	max_mapq = mapq;
      }
      if (mapq < min_mapq) {
	min_mapq = mapq;
      }
      mapq_tally[mapq] += 1;
    }
  
    printf("(");
    for (mapq = max_mapq; mapq >= min_mapq; mapq--) {
      if (mapq_tally[mapq] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%dMAPQ%d",mapq_tally[mapq],mapq);
      }
    }
    printf(")");
  }

  return;
}
#endif


static void
print_xs_list (Intlist_T xs_list) {
  int xs_tally[3];
  int xs;
  Intlist_T p;
  bool firstp = true;

  if (xs_list != NULL) {
    memset(xs_tally,0,3*sizeof(int));
    for (p = xs_list; p != NULL; p = Intlist_next(p)) {
      xs = Intlist_head(p);
      xs_tally[xs] += 1;
    }
  
    printf("(");
    for (xs = 2; xs >= 0; xs--) {
      if (xs_tally[xs] > 0) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  printf(",");
	}
	printf("%dXS",xs_tally[xs]);
	switch (xs) {
	case 0: printf("0"); break;
	case 1: printf("+"); break;
	case 2: printf("-"); break;
	default: abort();
	}
      }
    }
    printf(")");
  }

  return;
}

static Ucharlist_T
iit_xs_list (Ucharlist_T bytes, int *nbytes, Intlist_T xs_list) {
  int length;
  int xs_tally[3];
  int xs;
  Intlist_T p;

  if (xs_list != NULL) {
    memset(xs_tally,0,3*sizeof(int));
    for (p = xs_list; p != NULL; p = Intlist_next(p)) {
      xs = Intlist_head(p);
      xs_tally[xs] += 1;
    }
  
    length = 0;
    for (xs = 2; xs >= 0; xs--) {
      if (xs_tally[xs] > 0) {
	length += 1;
      }
    }
    bytes = push_int(&(*nbytes),bytes,length);

    for (xs = 2; xs >= 0; xs--) {
      if (xs_tally[xs] > 0) {
	bytes = push_int(&(*nbytes),bytes,xs);
	bytes = push_int(&(*nbytes),bytes,xs_tally[xs]);
      }
    }
  }

  return bytes;
}



typedef struct Codon_T *Codon_T;
struct Codon_T {
  char codoni;
  int tally_plus;
  int tally_minus;
  Intlist_T shift_list;
  Intlist_T nm_list;
  Intlist_T xs_list;
};

static void
Codon_free (Codon_T *old) {
  Intlist_free(&(*old)->xs_list);
  Intlist_free(&(*old)->nm_list);
  Intlist_free(&(*old)->shift_list);
  FREE(*old);
  return;
}

static Codon_T
Codon_new (char codoni, int tally_plus, int tally_minus, Intlist_T shift_list, Intlist_T nm_list, Intlist_T xs_list) {
  Codon_T new = (Codon_T) MALLOC(sizeof(*new));

  new->codoni = codoni;
  new->tally_plus = tally_plus;
  new->tally_minus = tally_minus;
  new->shift_list = shift_list;
  new->nm_list = nm_list;
  new->xs_list = xs_list;

  return new;
}


static int
Codon_cmp (const void *a, const void *b) {
  Codon_T x = * (Codon_T *) a;
  Codon_T y = * (Codon_T *) b;

  if (x->tally_plus + x->tally_minus > y->tally_plus + y->tally_minus) {
    return -1;
  } else if (y->tally_plus + y->tally_minus > x->tally_plus + x->tally_minus) {
    return +1;
  } else if (x->codoni < y->codoni) {
    return -1;
  } else if (y->codoni < x->codoni) {
    return +1;
  } else {
    return 0;
  }
}


static Ucharlist_T
process_codons_plus (Ucharlist_T bytes, int *nbytes, Tally_T tally0, Tally_T tally1, Tally_T tally2,
		     Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		     Genome_T genome, Genomicpos_T chroffset, bool signed_counts_p, bool print_cycles_p,
		     bool print_nm_scores_p, bool print_xs_scores_p, int quality_score_adj,
		     Tally_outputtype_T output_type) {
  int total_plus = 0, total_minus = 0;
  int codon_tally_plus[64], codon_tally_minus[64], ref_codon_tally_plus = 0, ref_codon_tally_minus = 0;
  Readevid_T *readevid0, *readevid1, *readevid2;
  int nreads0, nreads1, nreads2;
  int i = 0, j = 0, k = 0;
  unsigned int linei0, linei1, linei2, highest;
  char codoni, ref_codoni;
  int shift, nm, xs;
  List_T alt_codons;
  Intlist_T codon_shift[64], codon_nm[64], codon_xs[64],
    ref_codon_shift = NULL, ref_codon_nm = NULL, ref_codon_xs = NULL;
  Codon_T *array;
  int ncodons;


  readevid0 = (Readevid_T *) List_to_array_n(&nreads0,tally0->readevidence);
  readevid1 = (Readevid_T *) List_to_array_n(&nreads1,tally1->readevidence);
  readevid2 = (Readevid_T *) List_to_array_n(&nreads2,tally2->readevidence);

  if ((ref_codoni = Tally_codoni_plus(tally0,tally1,tally2,chrpos0,chrpos1,chrpos2,genome,chroffset)) < 0) {
    /* Skip, because genome has one or more Ns here */

  } else if (nreads0 == 0 || nreads1 == 0 || nreads2 == 0) {
    /* Lack of read evidence, due to indel or non-overlapping reads */
    if (output_type == OUTPUT_IIT) {
      debug2(printf("Total signed aa counts:\n"));
      bytes = push_int(&(*nbytes),bytes,/*total_plus*/0);
      bytes = push_int(&(*nbytes),bytes,/*total_minus*/0);

      /* Reference aa and counts */
      bytes = push_codon(&(*nbytes),bytes,ref_codoni);
      bytes = push_int(&(*nbytes),bytes,/*ref_codon_tally_plus*/0);
      bytes = push_int(&(*nbytes),bytes,/*ref_codon_tally_minus*/0);

    } else {
      printf(" %c[%s]",aa_table[ref_codoni],codon_table[ref_codoni]);
      if (signed_counts_p == false) {
	printf("0");
      } else {
	printf("0|0");
      }
    }

  } else {
    memset(codon_tally_plus,0,64*sizeof(int));
    memset(codon_tally_minus,0,64*sizeof(int));
    memset(codon_shift,0,64*sizeof(Intlist_T));
    memset(codon_nm,0,64*sizeof(Intlist_T));
    memset(codon_xs,0,64*sizeof(Intlist_T));

#if 0
    /* No need to sort, since they should be in order from highest readi to lowest */
    qsort(readevid0,nreads0,sizeof(Readevid_T),Readevid_cmp);
    qsort(readevid1,nreads1,sizeof(Readevid_T),Readevid_cmp);
    qsort(readevid2,nreads2,sizeof(Readevid_T),Readevid_cmp);
#endif

    while (i < nreads0 && j < nreads1 && k < nreads2) {
      linei0 = Readevid_linei(readevid0[i]);
      linei1 = Readevid_linei(readevid1[j]);
      linei2 = Readevid_linei(readevid2[k]);

      highest = linei0;
      if (linei1 > highest) highest = linei1;
      if (linei2 > highest) highest = linei2;

      if (linei0 == highest && linei1 == highest && linei2 == highest) {
	if ((codoni = Readevid_codoni_plus(&shift,&nm,&xs,readevid0[i],readevid1[j],readevid2[k])) < 0) {
	  /* Skip */
	} else {
	  if (shift > 0) {
	    codon_tally_plus[codoni] += 1;
	    total_plus += 1;
	  } else {
	    codon_tally_minus[codoni] += 1;
	    total_minus += 1;
	  }
	  codon_shift[codoni] = Intlist_push(codon_shift[codoni],shift);
	  codon_nm[codoni] = Intlist_push(codon_nm[codoni],nm);
	  codon_xs[codoni] = Intlist_push(codon_xs[codoni],xs);
	}
	/* printf(" %u:%c%c%c",highest,Readevid_nt(readevid0[i]),Readevid_nt(readevid1[j]),Readevid_nt(readevid2[k])); */
	i++; j++; k++;

      } else {
	if (linei0 == highest) i++;
	if (linei1 == highest) j++;
	if (linei2 == highest) k++;
      }
    }

    alt_codons = (List_T) NULL;
    for (codoni = 0; codoni < 64; codoni++) {
      if (codon_tally_plus[codoni] > 0 || codon_tally_minus[codoni] > 0) {
	if (codoni == ref_codoni) {
	  ref_codon_tally_plus = codon_tally_plus[codoni];
	  ref_codon_tally_minus = codon_tally_minus[codoni];
	  ref_codon_shift = codon_shift[codoni];
	  ref_codon_nm = codon_nm[codoni];
	  ref_codon_xs = codon_xs[codoni];
	} else {
	  alt_codons = List_push(alt_codons,(void *) Codon_new(codoni,codon_tally_plus[codoni],codon_tally_minus[codoni],
							       codon_shift[codoni],codon_nm[codoni],codon_xs[codoni]));
	}
      }
    }

    if (ref_codon_tally_plus + ref_codon_tally_minus > 0 || alt_codons != NULL) {
      if (output_type == OUTPUT_IIT) {
	debug2(printf("Total signed codon counts on plus strand:\n"));
	bytes = push_int(&(*nbytes),bytes,total_plus);
	bytes = push_int(&(*nbytes),bytes,total_minus);

	/* Plus-strand reference codon and counts */
	debug2(printf("Reference codon counts on plus strand:\n"));
	bytes = push_codon(&(*nbytes),bytes,ref_codoni);
	bytes = push_int(&(*nbytes),bytes,ref_codon_tally_plus);
	bytes = push_int(&(*nbytes),bytes,ref_codon_tally_minus);

      } else {
	printf(" %c[%s]",aa_table[ref_codoni],codon_table[ref_codoni]);
	if (signed_counts_p == false) {
	  printf("%d",ref_codon_tally_plus + ref_codon_tally_minus);
	} else {
	  printf("%d|%d",ref_codon_tally_plus,ref_codon_tally_minus);
	}
	if (print_cycles_p == true) {
	  print_shift_list(ref_codon_shift);
	}
	if (print_nm_scores_p == true) {
	  print_nm_list(ref_codon_nm);
	}
	if (print_xs_scores_p == true) {
	  print_xs_list(ref_codon_xs);
	}

	Intlist_free(&ref_codon_shift);
	Intlist_free(&ref_codon_nm);
	Intlist_free(&ref_codon_xs);
      }
      
      if (output_type == OUTPUT_IIT) {
	if (alt_codons != NULL) {
	  array = (Codon_T *) List_to_array_n(&ncodons,alt_codons);
	  qsort(array,ncodons,sizeof(Codon_T),Codon_cmp);

	  for (i = 0; i < ncodons; i++) {
	    debug2(printf("Counts on plus strand for alternate codon #%d:\n",i));
	    codoni = array[i]->codoni;
	    bytes = push_codon(&(*nbytes),bytes,codoni);
	    bytes = push_int(&(*nbytes),bytes,array[i]->tally_plus);
	    bytes = push_int(&(*nbytes),bytes,array[i]->tally_minus);
	  }
	}

	/* Terminates plus-strand amino acids */
	bytes = push_codon(&(*nbytes),bytes,255);

      } else {
	if (alt_codons != NULL) {
	  array = (Codon_T *) List_to_array_n(&ncodons,alt_codons);
	  qsort(array,ncodons,sizeof(Codon_T),Codon_cmp);

	  for (i = 0; i < ncodons; i++) {
	    codoni = array[i]->codoni;
	    printf(" %c[%s]",aa_table[codoni],codon_table[codoni]);
	    if (signed_counts_p == false) {
	      printf("%d",array[i]->tally_plus + array[i]->tally_minus);
	    } else {
	      printf("%d|%d",array[i]->tally_plus,array[i]->tally_minus);
	    }
	    if (print_cycles_p == true) {
	      print_shift_list(array[i]->shift_list);
	    }
	    if (print_nm_scores_p == true) {
	      print_nm_list(array[i]->nm_list);
	    }
	    if (print_xs_scores_p == true) {
	      print_xs_list(array[i]->xs_list);
	    }

	    Codon_free(&(array[i]));
	  }
	  
	  FREE(array);
	  List_free(&alt_codons);
	}
      }

      if (output_type == OUTPUT_IIT) {
	/* Output cycles, nm, and xs */
	debug2(printf("Cycles/nm/xs for plus-strand reference codon:\n"));
	if (print_cycles_p == true) {
	  bytes = iit_shift_list(bytes,&(*nbytes),ref_codon_shift);
	}
	if (print_nm_scores_p == true) {
	  bytes = iit_nm_list(bytes,&(*nbytes),ref_codon_nm);
	}
#if 0
	bytes = iit_quality_list(bytes,&(*nbytes),ref_codon_quality,quality_score_adj);
#endif
	if (print_xs_scores_p == true) {
	  bytes = iit_xs_list(bytes,&(*nbytes),ref_codon_xs);
	}

	Intlist_free(&ref_codon_shift);
	Intlist_free(&ref_codon_nm);
	Intlist_free(&ref_codon_xs);

	if (alt_codons != NULL) {
	  for (i = 0; i < ncodons; i++) {
	    debug2(printf("Cycles/nm/xs for plus-strand alternate codon #%d:\n",i));
	    if (print_cycles_p == true) {
	      bytes = iit_shift_list(bytes,&(*nbytes),array[i]->shift_list);
	    }
	    if (print_nm_scores_p == true) {
	      bytes = iit_nm_list(bytes,&(*nbytes),array[i]->nm_list);
	    }
#if 0
	    bytes = iit_quality_list(bytes,&(*nbytes),array[i]->quality_list,quality_score_adj);
#endif
	    if (print_xs_scores_p == true) {
	      bytes = iit_xs_list(bytes,&(*nbytes),array[i]->xs_list);
	    }
	    Codon_free(&(array[i]));
	  }

	  FREE(array);
	  List_free(&alt_codons);
	}
      }
    }
  }

  FREE(readevid2);
  FREE(readevid1);
  FREE(readevid0);

  return bytes;
}


static Ucharlist_T
process_codons_minus (Ucharlist_T bytes, int *nbytes, Tally_T tally0, Tally_T tally1, Tally_T tally2,
		      Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		      Genome_T genome, Genomicpos_T chroffset, bool signed_counts_p, bool print_cycles_p,
		      bool print_nm_scores_p, bool print_xs_scores_p, int quality_score_adj, Tally_outputtype_T output_type) {
  int total_plus = 0, total_minus = 0;
  int codon_tally_plus[64], codon_tally_minus[64], ref_codon_tally_plus = 0, ref_codon_tally_minus = 0;
  Readevid_T *readevid0, *readevid1, *readevid2;
  int nreads0, nreads1, nreads2;
  int i = 0, j = 0, k = 0;
  unsigned int linei0, linei1, linei2, highest;
  char codoni, ref_codoni;
  int shift, nm, xs;
  List_T alt_codons;
  Intlist_T codon_shift[64], codon_nm[64], codon_xs[64],
    ref_codon_shift = NULL, ref_codon_nm = NULL, ref_codon_xs = NULL;
  Codon_T *array;
  int ncodons;

  readevid0 = (Readevid_T *) List_to_array_n(&nreads0,tally0->readevidence);
  readevid1 = (Readevid_T *) List_to_array_n(&nreads1,tally1->readevidence);
  readevid2 = (Readevid_T *) List_to_array_n(&nreads2,tally2->readevidence);

  if ((ref_codoni = Tally_codoni_minus(tally0,tally1,tally2,chrpos0,chrpos1,chrpos2,genome,chroffset)) < 0) {
    /* Skip, because genome has one or more Ns here */

  } else if (nreads0 == 0 || nreads1 == 0 || nreads2 == 0) {
    /* Lack of read evidence, due to indel or non-overlapping reads */
    if (output_type == OUTPUT_IIT) {
      debug2(printf("Total signed aa counts:\n"));
      bytes = push_int(&(*nbytes),bytes,/*total_plus*/0);
      bytes = push_int(&(*nbytes),bytes,/*total_minus*/0);

      /* Reference aa and counts */
      bytes = push_codon(&(*nbytes),bytes,ref_codoni);
      bytes = push_int(&(*nbytes),bytes,/*ref_codon_tally_plus*/0);
      bytes = push_int(&(*nbytes),bytes,/*ref_codon_tally_minus*/0);

    } else {
      printf(" %c[%s]",aa_table[ref_codoni],codon_table[ref_codoni]);
      if (signed_counts_p == false) {
	printf("0");
      } else {
	printf("0|0");
      }
    }

  } else {
    memset(codon_tally_plus,0,64*sizeof(int));
    memset(codon_tally_minus,0,64*sizeof(int));
    memset(codon_shift,0,64*sizeof(Intlist_T));
    memset(codon_nm,0,64*sizeof(Intlist_T));
    memset(codon_xs,0,64*sizeof(Intlist_T));

#if 0
    /* No need to sort, since they should be in order from highest readi to lowest */
    qsort(readevid0,nreads0,sizeof(Readevid_T),Readevid_cmp);
    qsort(readevid1,nreads1,sizeof(Readevid_T),Readevid_cmp);
    qsort(readevid2,nreads2,sizeof(Readevid_T),Readevid_cmp);
#endif

    while (i < nreads0 && j < nreads1 && k < nreads2) {
      linei0 = Readevid_linei(readevid0[i]);
      linei1 = Readevid_linei(readevid1[j]);
      linei2 = Readevid_linei(readevid2[k]);

      highest = linei0;
      if (linei1 > highest) highest = linei1;
      if (linei2 > highest) highest = linei2;

      if (linei0 == highest && linei1 == highest && linei2 == highest) {
	if ((codoni = Readevid_codoni_minus(&shift,&nm,&xs,readevid0[i],readevid1[j],readevid2[k])) < 0) {
	  /* Skip */
	} else {
	  if (shift > 0) {
	    codon_tally_plus[codoni] += 1;
	    total_plus += 1;
	  } else {
	    codon_tally_minus[codoni] += 1;
	    total_minus += 1;
	  }
	  codon_shift[codoni] = Intlist_push(codon_shift[codoni],shift);
	  codon_nm[codoni] = Intlist_push(codon_nm[codoni],nm);
	  codon_xs[codoni] = Intlist_push(codon_xs[codoni],xs);
	}
	/* printf(" %u:%c%c%c",highest,Readevid_nt(readevid0[i]),Readevid_nt(readevid1[j]),Readevid_nt(readevid2[k])); */
	i++; j++; k++;

      } else {
	if (linei0 == highest) i++;
	if (linei1 == highest) j++;
	if (linei2 == highest) k++;
      }
    }

    alt_codons = (List_T) NULL;
    for (codoni = 0; codoni < 64; codoni++) {
      if (codon_tally_plus[codoni] > 0 || codon_tally_minus[codoni] > 0) {
	if (codoni == ref_codoni) {
	  ref_codon_tally_plus = codon_tally_plus[codoni];
	  ref_codon_tally_minus = codon_tally_minus[codoni];
	  ref_codon_shift = codon_shift[codoni];
	  ref_codon_nm = codon_nm[codoni];
	  ref_codon_xs = codon_xs[codoni];
	} else {
	  alt_codons = List_push(alt_codons,(void *) Codon_new(codoni,codon_tally_plus[codoni],codon_tally_minus[codoni],
							       codon_shift[codoni],codon_nm[codoni],codon_xs[codoni]));
	}
      }
    }


    if (ref_codon_tally_plus + ref_codon_tally_minus > 0 || alt_codons != NULL) {
      if (output_type == OUTPUT_IIT) {
	debug2(printf("Total signed codon counts on minus strand:\n"));
	bytes = push_int(&(*nbytes),bytes,total_plus);
	bytes = push_int(&(*nbytes),bytes,total_minus);

	/* Minus-strand reference codon and counts */
	debug2(printf("Reference codon counts on minus strand:\n"));
	bytes = push_codon(&(*nbytes),bytes,ref_codoni);
	bytes = push_int(&(*nbytes),bytes,ref_codon_tally_plus);
	bytes = push_int(&(*nbytes),bytes,ref_codon_tally_minus);

      } else {
	printf(" %c[%s]",aa_table[ref_codoni],codon_table[ref_codoni]);
	if (signed_counts_p == false) {
	  printf("%d",ref_codon_tally_plus + ref_codon_tally_minus);
	} else {
	  printf("%d|%d",ref_codon_tally_plus,ref_codon_tally_minus);
	}
	if (print_cycles_p == true) {
	  print_shift_list(ref_codon_shift);
	}
	if (print_nm_scores_p == true) {
	  print_nm_list(ref_codon_nm);
	}
	if (print_xs_scores_p == true) {
	  print_xs_list(ref_codon_xs);
	}

	Intlist_free(&ref_codon_shift);
	Intlist_free(&ref_codon_nm);
	Intlist_free(&ref_codon_xs);
      }

      if (output_type == OUTPUT_IIT) {
	if (alt_codons != NULL) {
	  array = (Codon_T *) List_to_array_n(&ncodons,alt_codons);
	  qsort(array,ncodons,sizeof(Codon_T),Codon_cmp);

	  for (i = 0; i < ncodons; i++) {
	    debug2(printf("Counts on minus strand for alternate codon #%d:\n",i));
	    codoni = array[i]->codoni;
	    bytes = push_codon(&(*nbytes),bytes,codoni);
	    bytes = push_int(&(*nbytes),bytes,array[i]->tally_plus);
	    bytes = push_int(&(*nbytes),bytes,array[i]->tally_minus);
	  }
	}
	  
	/* Terminates minus-strand amino acids */
	bytes = push_codon(&(*nbytes),bytes,255);

      } else {
	if (alt_codons != NULL) {
	  array = (Codon_T *) List_to_array_n(&ncodons,alt_codons);
	  qsort(array,ncodons,sizeof(Codon_T),Codon_cmp);

	  for (i = 0; i < ncodons; i++) {
	    codoni = array[i]->codoni;
	    printf(" %c[%s]",aa_table[codoni],codon_table[codoni]);
	    if (signed_counts_p == false) {
	      printf("%d",array[i]->tally_plus + array[i]->tally_minus);
	    } else {
	      printf("%d|%d",array[i]->tally_plus,array[i]->tally_minus);
	    }
	    if (print_cycles_p == true) {
	      print_shift_list(array[i]->shift_list);
	    }
	    if (print_nm_scores_p == true) {
	      print_nm_list(array[i]->nm_list);
	    }
	    if (print_xs_scores_p == true) {
	      print_xs_list(array[i]->xs_list);
	    }

	    Codon_free(&(array[i]));
	  }
	  
	  FREE(array);
	  List_free(&alt_codons);
	}
      }

      if (output_type == OUTPUT_IIT) {
	/* Output cycles, nm, and xs */
	debug2(printf("Cycles/nm/xs for minus-strand reference codon:\n"));
	if (print_cycles_p == true) {
	  bytes = iit_shift_list(bytes,&(*nbytes),ref_codon_shift);
	}
	if (print_nm_scores_p == true) {
	  bytes = iit_nm_list(bytes,&(*nbytes),ref_codon_nm);
	}
#if 0
	bytes = iit_quality_list(bytes,&(*nbytes),ref_codon_quality,quality_score_adj);
#endif
	if (print_xs_scores_p == true) {
	  bytes = iit_xs_list(bytes,&(*nbytes),ref_codon_xs);
	}

	Intlist_free(&ref_codon_shift);
	Intlist_free(&ref_codon_nm);
	Intlist_free(&ref_codon_xs);

	if (alt_codons != NULL) {
	  for (i = 0; i < ncodons; i++) {
	    debug2(printf("Cycles/mismatches/xs for minus-strand alternate codon #%d:\n",i));
	    if (print_cycles_p == true) {
	      bytes = iit_shift_list(bytes,&(*nbytes),array[i]->shift_list);
	    }
	    if (print_nm_scores_p == true) {
	      bytes = iit_nm_list(bytes,&(*nbytes),array[i]->nm_list);
	    }
#if 0
	    bytes = iit_quality_list(bytes,&(*nbytes),array[i]->quality_list,quality_score_adj);
#endif
	    if (print_xs_scores_p == true) {
	      bytes = iit_xs_list(bytes,&(*nbytes),array[i]->xs_list);
	    }
	    Codon_free(&(array[i]));
	  }

	  FREE(array);
	  List_free(&alt_codons);
	}
      }
    }
  }

  FREE(readevid2);
  FREE(readevid1);
  FREE(readevid0);

  return bytes;
}


static Ucharlist_T
report_plus_genes (Ucharlist_T bytes, int *nbytes, Tally_T tally2, Tally_T *block_tallies, Genomicpos_T blockstart,
		   Tally_T *prev_block_tallies, Genomicpos_T prev_blockstart, Genomicpos_T prev_blockptr,
		   List_T plus_genes, Genomicpos_T chrpos, Genome_T genome, Genomicpos_T chroffset,
		   bool signed_counts_p, bool print_cycles_p, bool print_nm_scores_p, bool print_xs_scores_p, int quality_score_adj,
		   Tally_outputtype_T output_type) {
  List_T chrpos_pairs = NULL;
  Chrpos_pair_T chrpos_pair;
  Tally_T tally0, tally1;
  List_T p;
  Gene_T gene;
  int exoni;
  int ntpos;

  char *gene_sequence, *aa_sequence;
  int genelength, translation_length;
  int ntframe;

  int ignore;


  for (p = plus_genes; p != NULL; p = List_next(p)) {
    gene = (Gene_T) List_head(p);
    exoni = gene->exoni;
    while (exoni < gene->nexons && chrpos > gene->exons[exoni]->exonend) {
      exoni++;
    }

    if (exoni >= gene->nexons) {
      /* Skip */
    } else if (chrpos < gene->exons[exoni]->exonstart) {
      /* Skip */
    } else {
      if (gene->translation_end < 0) {
	gene_sequence = concatenate_exons(&genelength,gene->exons,gene->nexons);
	aa_sequence = Translation_via_genomic(&gene->translation_start,&gene->translation_end,&translation_length,
					      gene_sequence,/*startpos*/0,/*endpos*/genelength,genelength,
					      /*backwardp*/false,/*revcompp*/false,/*fulllengthp*/true,/*cds_startpos*/0);
	FREE(aa_sequence);
	FREE(gene_sequence);
      }

      ntpos = gene->cum_exonlength[exoni] + (chrpos - gene->exons[exoni]->exonstart);
      if (ntpos >= gene->translation_start && ntpos <= gene->translation_end) {
	if ((ntframe = (ntpos - gene->translation_start) % 3) == 2) {
	  if (chrpos - 2 >= gene->exons[exoni]->exonstart) {
	    chrpos_pairs = List_push(chrpos_pairs,(void *) Chrpos_pair_new(chrpos-2,chrpos-1));

	  } else if (exoni - 1 < 0) {
	    fprintf(stderr,"Need to compute coding regions\n");

	  } else if (chrpos - 1 >= gene->exons[exoni]->exonstart) {
	    chrpos_pairs = List_push(chrpos_pairs,
				     (void *) Chrpos_pair_new(gene->exons[exoni-1]->exonend,chrpos-1));

	  } else {
	    chrpos_pairs = List_push(chrpos_pairs,
				     (void *) Chrpos_pair_new(gene->exons[exoni-1]->exonend-1,
							      gene->exons[exoni-1]->exonend));
	  }
	}

	/* printf(" (frame %d)",ntframe); */
      }
    }

    gene->exoni = exoni;
  }

  
  /* Process codons */
  chrpos_pairs = Chrpos_pair_uniq(chrpos_pairs);
  for (p = chrpos_pairs; p != NULL; p = List_next(p)) {
    chrpos_pair = (Chrpos_pair_T) List_head(p);
    if (chrpos_pair->frame0pos >= blockstart) {
      tally0 = block_tallies[chrpos_pair->frame0pos - blockstart];
    } else if (chrpos_pair->frame0pos >= prev_blockstart && chrpos_pair->frame0pos < prev_blockptr) {
      tally0 = prev_block_tallies[chrpos_pair->frame0pos - prev_blockstart];
    } else {
      tally0 = (Tally_T) NULL;
    }

    if (chrpos_pair->frame1pos >= blockstart) {
      tally1 = block_tallies[chrpos_pair->frame1pos - blockstart];
    } else if (chrpos_pair->frame1pos >= prev_blockstart && chrpos_pair->frame1pos < prev_blockptr) {
      tally1 = prev_block_tallies[chrpos_pair->frame1pos - prev_blockstart];
    } else {
      tally1 = (Tally_T) NULL;
    }

    if (tally0 != NULL && tally1 != NULL) {
      if (output_type == OUTPUT_IIT) {
	bytes = process_codons_plus(bytes,&(*nbytes),tally0,tally1,tally2,
				    chrpos_pair->frame0pos,chrpos_pair->frame1pos,/*chrpos2*/chrpos,
				    genome,chroffset,signed_counts_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				    quality_score_adj,output_type);
      } else {
	process_codons_plus(/*bytes*/NULL,&ignore,tally0,tally1,tally2,
			    chrpos_pair->frame0pos,chrpos_pair->frame1pos,/*chrpos2*/chrpos,
			    genome,chroffset,signed_counts_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
			    quality_score_adj,output_type);
      }
    }
  }

  for (p = chrpos_pairs; p != NULL; p = List_next(p)) {
    chrpos_pair = (Chrpos_pair_T) List_head(p);
    Chrpos_pair_free(&chrpos_pair);
  }
  List_free(&chrpos_pairs);

  return bytes;
}


static Ucharlist_T
report_minus_genes (Ucharlist_T bytes, int *nbytes, Tally_T tally2, Tally_T *block_tallies, Genomicpos_T blockstart,
		    Tally_T *prev_block_tallies, Genomicpos_T prev_blockstart, Genomicpos_T prev_blockptr,
		    List_T minus_genes, Genomicpos_T chrpos, Genome_T genome, Genomicpos_T chroffset,
		    bool signed_counts_p, bool print_cycles_p, bool print_nm_scores_p, bool print_xs_scores_p, int quality_score_adj,
		    Tally_outputtype_T output_type) {
  List_T chrpos_pairs = NULL;
  Chrpos_pair_T chrpos_pair;
  Tally_T tally0, tally1;
  List_T p;
  Gene_T gene;
  int exoni;
  int ntpos;

  char *gene_sequence, *aa_sequence;
  int genelength, translation_length;
  int ntframe;

  int ignore;


  for (p = minus_genes; p != NULL; p = List_next(p)) {
    gene = (Gene_T) List_head(p);
    exoni = gene->exoni;
    while (exoni >= 0 && chrpos > gene->exons[exoni]->exonstart) {
      exoni--;
    }

    if (exoni < 0) {
      /* Skip */
    } else if (chrpos < gene->exons[exoni]->exonend) {
      /* Skip */
    } else {
      if (gene->translation_end < 0) {
	gene_sequence = concatenate_exons(&genelength,gene->exons,gene->nexons);
	aa_sequence = Translation_via_genomic(&gene->translation_start,&gene->translation_end,&translation_length,
					      gene_sequence,/*startpos*/0,/*endpos*/genelength,genelength,
					      /*backwardp*/false,/*revcompp*/false,/*fulllengthp*/true,/*cds_startpos*/0);
	FREE(aa_sequence);
	FREE(gene_sequence);
      }

      ntpos = gene->cum_exonlength[exoni] + (gene->exons[exoni]->exonstart - chrpos);
      if (ntpos >= gene->translation_start && ntpos < gene->translation_end) { /* Needs to be < and not <= */
	if ((ntframe = (ntpos - gene->translation_start) % 3) == 0) {
	  if (chrpos - 2 >= gene->exons[exoni]->exonend) {
	    chrpos_pairs = List_push(chrpos_pairs,(void *) Chrpos_pair_new(chrpos-2,chrpos-1));

	  } else if (exoni + 1 >= gene->nexons) {
	    fprintf(stderr,"Need to compute coding regions\n");

	  } else if (chrpos - 1 >= gene->exons[exoni]->exonend) {
	    chrpos_pairs = List_push(chrpos_pairs,
				     (void *) Chrpos_pair_new(gene->exons[exoni+1]->exonstart,chrpos-1));

	  } else {
	    chrpos_pairs = List_push(chrpos_pairs,
				     (void *) Chrpos_pair_new(gene->exons[exoni+1]->exonstart-1,
							      gene->exons[exoni+1]->exonstart));
	  }
	}

	/* printf(" (frame %d)",ntframe); */
      }
    }

    gene->exoni = exoni;
  }

  /* Process codons */
  chrpos_pairs = Chrpos_pair_uniq(chrpos_pairs);
  for (p = chrpos_pairs; p != NULL; p = List_next(p)) {
    chrpos_pair = (Chrpos_pair_T) List_head(p);
    if (chrpos_pair->frame0pos >= blockstart) {
      tally0 = block_tallies[chrpos_pair->frame0pos - blockstart];
    } else if (chrpos_pair->frame0pos >= prev_blockstart && chrpos_pair->frame0pos < prev_blockptr) {
      tally0 = prev_block_tallies[chrpos_pair->frame0pos - prev_blockstart];
    } else {
      tally0 = (Tally_T) NULL;
    }

    if (chrpos_pair->frame1pos >= blockstart) {
      tally1 = block_tallies[chrpos_pair->frame1pos - blockstart];
    } else if (chrpos_pair->frame1pos >= prev_blockstart && chrpos_pair->frame1pos < prev_blockptr) {
      tally1 = prev_block_tallies[chrpos_pair->frame1pos - prev_blockstart];
    } else {
      tally1 = (Tally_T) NULL;
    }
	      
    if (tally0 != NULL && tally1 != NULL) {
      if (output_type == OUTPUT_IIT) {
	bytes = process_codons_minus(bytes,&(*nbytes),tally0,tally1,tally2,
				     chrpos_pair->frame0pos,chrpos_pair->frame1pos,/*chrpos2*/chrpos,
				     genome,chroffset,signed_counts_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				     quality_score_adj,output_type);
      } else {
	process_codons_minus(/*bytes*/NULL,&ignore,tally0,tally1,tally2,
			     chrpos_pair->frame0pos,chrpos_pair->frame1pos,/*chrpos2*/chrpos,
			     genome,chroffset,signed_counts_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
			     quality_score_adj,output_type);
      }
    }
  }

  for (p = chrpos_pairs; p != NULL; p = List_next(p)) {
    chrpos_pair = (Chrpos_pair_T) List_head(p);
    Chrpos_pair_free(&chrpos_pair);
  }
  List_free(&chrpos_pairs);

  return bytes;
}


static void
print_gene_info (List_T plus_genes, List_T minus_genes, Genomicpos_T chrpos) {
  List_T p;
  Gene_T gene;
  int exoni;
  int ntpos;
  bool firstp = true;

  /* Print gene information */
  for (p = plus_genes; p != NULL; p = List_next(p)) {
    gene = (Gene_T) List_head(p);
    exoni = gene->exoni;
#if 0
    /* Already done by report_plus_genes */
    while (exoni < gene->nexons && chrpos > gene->exons[exoni]->exonend) {
      exoni++;
    }
#endif

    if (exoni >= gene->nexons) {
      /* Skip */
    } else if (chrpos < gene->exons[exoni]->exonstart) {
      /* Skip */
    } else {
      if (firstp == true) {
	firstp = false;
      } else {
	printf("|");
      }
      ntpos = gene->cum_exonlength[exoni] + (chrpos - gene->exons[exoni]->exonstart);
      printf("+%s_%s_exon%d/%d_nt%d",gene->genename,gene->acc,exoni+1,gene->nexons,ntpos+1);
    }
  }

  for (p = minus_genes; p != NULL; p = List_next(p)) {
    gene = (Gene_T) List_head(p);
    exoni = gene->exoni;
#if 0
    /* Already done by report_minus_genes */
    while (exoni >= 0 && chrpos > gene->exons[exoni]->exonstart) {
      exoni--;
    }
#endif

    if (exoni < 0) {
      /* Skip */
    } else if (chrpos < gene->exons[exoni]->exonend) {
      /* Skip */
    } else {
      if (firstp == true) {
	firstp = false;
      } else {
	printf("|");
      }
      ntpos = gene->cum_exonlength[exoni] + (gene->exons[exoni]->exonstart - chrpos);
      printf("-%s_%s_exon%d/%d_nt%d",gene->genename,gene->acc,exoni+1,gene->nexons,ntpos+1);
    }
  }

  return;
}



static void
genes_get (List_T *plus_genes, List_T *minus_genes, IIT_T map_iit,
	   char *chr, Genomicpos_T chroffset, Genome_T genome,
	   unsigned int mincoord, unsigned int maxcoord) {
  int *matches, nmatches, delcounts_plus, delcounts_minus, index, i;
  Exon_T *exons;
  int nexons;
  char *acc, *label, *genename, *annotation, *restofheader, *p;
  int gene_namelength;
  bool alloc1p, alloc2p;
  int k;

  *plus_genes = (List_T) NULL;
  matches = IIT_get_signed(&nmatches,map_iit,chr,mincoord,maxcoord,/*sign*/+1,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    /* interval = IIT_interval(iit,index); */
    label = IIT_label(map_iit,index,&alloc1p);
    acc = (char *) CALLOC(strlen(label)+1,sizeof(char));
    strcpy(acc,label);

    annotation = IIT_annotation(&restofheader,map_iit,index,&alloc2p);

    /* Get genename from annotation */
    p = annotation;
    while (*p != '\0' && *p != '\n' && !isspace(*p)) {
      p++;
    }
    gene_namelength = p - annotation;
    genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
    strncpy(genename,annotation,gene_namelength);
    for (p = genename, k = 0; k < gene_namelength; p++, k++) {
      if (*p == '-') {
	*p = '.';
      }
    }

    exons = get_exons(&nexons,annotation,chr,chroffset,/*sign*/+1,genome);
    *plus_genes = List_push(*plus_genes,Gene_new(acc,genename,exons,nexons,/*sign*/+1,mincoord));
    if (alloc2p) {
      FREE(annotation);
    }
    if (alloc1p) {
      FREE(acc);
    }
  }
  FREE(matches);

  *minus_genes = (List_T) NULL;
  matches = IIT_get_signed(&nmatches,map_iit,chr,mincoord,maxcoord,/*sign*/-1,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    /* interval = IIT_interval(iit,index); */
    label = IIT_label(map_iit,index,&alloc1p);
    acc = (char *) CALLOC(strlen(label)+1,sizeof(char));
    strcpy(acc,label);

    annotation = IIT_annotation(&restofheader,map_iit,index,&alloc2p);

    /* Get genename from annotation */
    p = annotation;
    while (*p != '\0' && *p != '\n' && !isspace(*p)) {
      p++;
    }
    gene_namelength = p - annotation;
    genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
    strncpy(genename,annotation,gene_namelength);
    for (p = genename, k = 0; k < gene_namelength; p++, k++) {
      if (*p == '-') {
	*p = '.';
      }
    }

    exons = get_exons(&nexons,annotation,chr,chroffset,/*sign*/-1,genome);
    *minus_genes = List_push(*minus_genes,Gene_new(acc,genename,exons,nexons,/*sign*/-1,mincoord));
    if (alloc2p) {
      FREE(annotation);
    }
    if (alloc1p) {
      FREE(acc);
    }
  }
  FREE(matches);

  return;
}



static void
print_block (Tally_T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr, 
	     Tally_T *prev_block_tallies, Genomicpos_T prev_blockstart, Genomicpos_T prev_blockptr,
	     Genome_T genome, char *printchr, Genomicpos_T chroffset, IIT_T map_iit,
	     bool blockp, int quality_score_adj, int min_depth, int variant_strands,
	     bool genomic_diff_p, bool signed_counts_p, bool print_totals_p, bool print_cycles_p,
	     bool print_nm_scores_p, bool print_xs_scores_p,
	     bool want_genotypes_p, bool readlevel_p, bool print_noncovered_p) {
  Tally_T this;
  double probs[NGENOTYPES], loglik[NGENOTYPES];
  int length, i, j, g, bestg;
  List_T plus_genes = NULL, minus_genes = NULL, p;
  Gene_T gene;
  Genomicpos_T chrpos;
  int blocki, lasti;
  Match_T *match_array, match;
  Mismatch_T mismatch, mismatch0, *mm_array, *mm_subarray;
  Insertion_T ins, ins0, *ins_array_byshift, *ins_array_bynm, *ins_array_byxs, *ins_subarray;
  Deletion_T del, del0, *del_array_byshift, *del_array_bynm, *del_array_byxs, *del_subarray;
  List_T unique_mismatches_byshift, unique_mismatches_bynm, unique_mismatches_byxs, ptr;
  List_T unique_insertions_byshift, unique_insertions_bynm, unique_insertions_byxs,
    unique_deletions_byshift, unique_deletions_bynm, unique_deletions_byxs;
  long int total, total_matches_plus, total_matches_minus, total_mismatches_plus, total_mismatches_minus;
  int shift, quality, nm, xs;
  int ignore;
  bool firstp;


  lasti = blockptr - blockstart;
  debug3(printf("Printing blocki 0 to %d\n",lasti));

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (print_noncovered_p == true ||
	(pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
			       this->mismatches_byshift,min_depth,variant_strands) == true &&
	 pass_difference_filter_p(probs,loglik,this,genome,
				  printchr,chroffset,chrpos,
				  quality_score_adj,genomic_diff_p) == true)) {
      total += this->nmatches + this->delcounts_plus + this->delcounts_minus;
      for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	mismatch = (Mismatch_T) List_head(ptr);
	total += mismatch->count;
      }
    }
  }
  
  if (blockstart == blockptr) {
    /* Skip */

  } else if (total <= 0) {
    /* Total could be 0 if the block is outside chrstart..chrend or if all positions fail variant filter */
    if (print_noncovered_p == true) {
      if (blockp == true) {
	printf(">%ld %s:%u..%u\n",total,printchr,blockstart,blockptr-1U);
      }

      for (blocki = 0; blocki < lasti; blocki++) {
	chrpos = blockstart + blocki;
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos);
	}
	printf("%c0\n",Genome_get_char(genome,chroffset+chrpos-1U));
      }
    }

  } else {
    if (blockp == true) {
      printf(">%ld %s:%u..%u\n",total,printchr,blockstart,blockptr-1U);
    }

    if (map_iit != NULL) {
      genes_get(&plus_genes,&minus_genes,map_iit,printchr,chroffset,genome,
		/*mincoord*/blockstart,/*maxcoord*/blockstart+lasti-1);
    }

    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;


      /* Handle insertions */
      if (this->insertions_byshift != NULL) {
	unique_insertions_byshift = NULL;
	for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byshift,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byshift = List_push(unique_insertions_byshift,ins0);
	  }
	}
	ins_array_byshift = (Insertion_T *) List_to_array(unique_insertions_byshift,NULL);
	qsort(ins_array_byshift,List_length(unique_insertions_byshift),sizeof(Insertion_T),Insertion_count_cmp);

	unique_insertions_bynm = NULL;
	for (ptr = this->insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_bynm,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_bynm = List_push(unique_insertions_bynm,ins0);
	  }
	}
	ins_array_bynm = (Insertion_T *) List_to_array(unique_insertions_bynm,NULL);
	qsort(ins_array_bynm,List_length(unique_insertions_bynm),sizeof(Insertion_T),Insertion_count_cmp);

	unique_insertions_byxs = NULL;
	for (ptr = this->insertions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byxs,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byxs = List_push(unique_insertions_byxs,ins0);
	  }
	}
	ins_array_byxs = (Insertion_T *) List_to_array(unique_insertions_byxs,NULL);
	qsort(ins_array_byxs,List_length(unique_insertions_byxs),sizeof(Insertion_T),Insertion_count_cmp);


	printf("^");
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos);
	}

	for (i = 0; i < List_length(unique_insertions_byshift); i++) {
	  ins0 = ins_array_byshift[i];
	  if (i > 0) {
	    printf(" ");
	  }
	  if (signed_counts_p == false) {
	    printf("%s %ld",ins0->segment,ins0->count_plus+ins0->count_minus);
	  } else {
	    printf("%s %ld|%ld",ins0->segment,ins0->count_plus,ins0->count_minus);
	  }
	  if (print_cycles_p == true) {
	    if ((length = Insertion_chain_length(ins0->next)) > 0) {
	      printf("(");
	      ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	      for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
		ins_subarray[j] = ins;
	      }
	      qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_shift_cmp);

	      printf("%ld@%d",ins_subarray[0]->count,ins_subarray[0]->shift);
	      for (j = 1; j < length; j++) {
		printf(",%ld@%d",ins_subarray[j]->count,ins_subarray[j]->shift);
	      }
	      FREE(ins_subarray);
	      printf(")");
	    }
	  }
	  if (print_nm_scores_p == true) {
	    ins0 = ins_array_bynm[i];
	    if ((length = Insertion_chain_length(ins0->next)) > 0) {
	      printf("(");
	      ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	      for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
		ins_subarray[j] = ins;
	      }
	      qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_nm_cmp);

	      printf("%ldNM%d",ins_subarray[0]->count,ins_subarray[0]->nm);
	      for (j = 1; j < length; j++) {
		printf(",%ldNM%d",ins_subarray[j]->count,ins_subarray[j]->nm);
	      }
	      FREE(ins_subarray);
	      printf(")");
	    }
	  }
	  if (print_xs_scores_p == true) {
	    ins0 = ins_array_byxs[i];
	    if ((length = Insertion_chain_length(ins0->next)) > 0) {
	      printf("(");
	      ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	      for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
		ins_subarray[j] = ins;
	      }
	      qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_xs_cmp);

	      printf("%ldXS%d",ins_subarray[0]->count,ins_subarray[0]->xs);
	      for (j = 1; j < length; j++) {
		printf(",%ldXS%d",ins_subarray[j]->count,ins_subarray[j]->xs);
	      }
	      FREE(ins_subarray);
	      printf(")");
	    }
	  }
	  if (signed_counts_p == false) {
	    printf(" ref:%ld",this->n_fromleft_plus+this->n_fromleft_minus);
	  } else {
	    printf(" ref:%ld|%ld",this->n_fromleft_plus,this->n_fromleft_minus);
	  }
	}

	printf("\n");

	FREE(ins_array_byxs);
	FREE(ins_array_bynm);
	FREE(ins_array_byshift);

	for (ptr = unique_insertions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byxs);

	for (ptr = unique_insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_bynm);

	for (ptr = unique_insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byshift);
      }


      /* Handle deletions */
      if (this->deletions_byshift != NULL) {
	unique_deletions_byshift = NULL;
	for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byshift,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byshift = List_push(unique_deletions_byshift,del0);
	  }
	}
	del_array_byshift = (Deletion_T *) List_to_array(unique_deletions_byshift,NULL);
	qsort(del_array_byshift,List_length(unique_deletions_byshift),sizeof(Deletion_T),Deletion_count_cmp);

	unique_deletions_bynm = NULL;
	for (ptr = this->deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_bynm,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_bynm = List_push(unique_deletions_bynm,del0);
	  }
	}
	del_array_bynm = (Deletion_T *) List_to_array(unique_deletions_bynm,NULL);
	qsort(del_array_bynm,List_length(unique_deletions_bynm),sizeof(Deletion_T),Deletion_count_cmp);

	unique_deletions_byxs = NULL;
	for (ptr = this->deletions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byxs,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byxs = List_push(unique_deletions_byxs,del0);
	  }
	}
	del_array_byxs = (Deletion_T *) List_to_array(unique_deletions_byxs,NULL);
	qsort(del_array_byxs,List_length(unique_deletions_byxs),sizeof(Deletion_T),Deletion_count_cmp);


	printf("_");
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos);
	}

	for (i = 0; i < List_length(unique_deletions_byshift); i++) {
	  del0 = del_array_byshift[i];
	  if (i > 0) {
	    printf(" ");
	  }
	  if (signed_counts_p == false) {
	    printf("%s %ld",del0->segment,del0->count_plus+del0->count_minus);
	  } else {
	    printf("%s %ld|%ld",del0->segment,del0->count_plus,del0->count_minus);
	  }
	  if (print_cycles_p == true) {
	    if ((length = Deletion_chain_length(del0->next)) > 0) {
	      printf("(");
	      del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	      for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
		del_subarray[j] = del;
	      }
	      qsort(del_subarray,length,sizeof(Deletion_T),Deletion_shift_cmp);
	      
	      printf("%ld@%d",del_subarray[0]->count,del_subarray[0]->shift);
	      for (j = 1; j < length; j++) {
		printf(",%ld@%d",del_subarray[j]->count,del_subarray[j]->shift);
	      }
	      FREE(del_subarray);
	      printf(")");
	    }
	  }
	  if (print_nm_scores_p == true) {
	    del0 = del_array_bynm[i];
	    if ((length = Deletion_chain_length(del0->next)) > 0) {
	      printf("(");
	      del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	      for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
		del_subarray[j] = del;
	      }
	      qsort(del_subarray,length,sizeof(Deletion_T),Deletion_nm_cmp);
	      
	      printf("%ldNM%d",del_subarray[0]->count,del_subarray[0]->nm);
	      for (j = 1; j < length; j++) {
		printf(",%ldNM%d",del_subarray[j]->count,del_subarray[j]->nm);
	      }
	      FREE(del_subarray);
	      printf(")");
	    }
          }
	  if (print_xs_scores_p == true) {
	    del0 = del_array_byxs[i];
	    if ((length = Deletion_chain_length(del0->next)) > 0) {
	      printf("(");
	      del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	      for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
		del_subarray[j] = del;
	      }
	      qsort(del_subarray,length,sizeof(Deletion_T),Deletion_xs_cmp);

	      printf("%ldXS%d",del_subarray[0]->count,del_subarray[0]->xs);
	      for (j = 1; j < length; j++) {
		printf(",%ldXS%d",del_subarray[j]->count,del_subarray[j]->xs);
	      }
	      FREE(del_subarray);
	      printf(")");
	    }
          }
	  if (signed_counts_p == false) {
	    printf(" ref:%ld",this->n_fromleft_plus+this->n_fromleft_minus);
	  } else {
	    printf(" ref:%ld|%ld",this->n_fromleft_plus,this->n_fromleft_minus);
	  }
	}

	printf("\n");

	FREE(del_array_byxs);
	FREE(del_array_bynm);
	FREE(del_array_byshift);

	for (ptr = unique_deletions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byxs);

	for (ptr = unique_deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_bynm);

	for (ptr = unique_deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byshift);
      }

      /* Handle mismatches */
      if (print_noncovered_p == false && 
	  pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
				this->mismatches_byshift,min_depth,variant_strands) == false) {
	if (blockp == true) {
	  printf("\n");
	}

      } else if (print_noncovered_p == false &&
		 pass_difference_filter_p(probs,loglik,this,genome,
					  printchr,chroffset,chrpos,
					  quality_score_adj,genomic_diff_p) == false) {
	if (blockp == true) {
	  printf("\n");
	}

      } else {
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",printchr,chrpos);
	}

	if (signed_counts_p == false) {
	  total = this->nmatches + this->delcounts_plus + this->delcounts_minus;
	  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    total += mismatch->count;
	  }
	  if (print_totals_p == true) {
	    printf("%ld\t",total);
	  }

	} else {
	  total_matches_plus = total_matches_minus = 0;
	  if (this->use_array_p == false) {
	    for (ptr = this->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	      match = (Match_T) List_head(ptr);
	      if (match->shift > 0) {
		total_matches_plus += match->count;
	      } else {
		total_matches_minus += match->count;
	      }
	    }
	  } else {
	    for (i = 0; i < this->n_matches_byshift_plus; i++) {
	      total_matches_plus += this->matches_byshift_plus[i];
	    }
	    for (i = 0; i < this->n_matches_byshift_minus; i++) {
	      total_matches_minus += this->matches_byshift_minus[i];
	    }
	  }
	  assert(this->nmatches == total_matches_plus + total_matches_minus);

	  total_mismatches_plus = total_mismatches_minus = 0;
	  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    if (mismatch->shift > 0) {
	      total_mismatches_plus += mismatch->count;
	    } else {
	      total_mismatches_minus += mismatch->count;
	    }
	  }
	  if (print_totals_p == true) {
	    printf("%ld|%ld\t",
		   total_matches_plus+this->delcounts_plus+total_mismatches_plus,
		   total_matches_minus+this->delcounts_minus+total_mismatches_minus);
	  }
	}


	/* Details */
	if (readlevel_p == true || totals_only_p == true) {
	  /* Don't print details */

	} else if (this->nmatches + this->delcounts_plus + this->delcounts_minus == 0) {
	  /* Genome_T expects zero-based numbers */
	  if (signed_counts_p == false) {
	    printf("%c0",Genome_get_char(genome,chroffset+chrpos-1U));
	  } else {
	    printf("%c0|0",Genome_get_char(genome,chroffset+chrpos-1U));
	  }
	} else {
	  if (signed_counts_p == false) {
	    printf("%c%d",this->refnt,this->nmatches);
	  } else {
	    printf("%c%ld|%ld",this->refnt,total_matches_plus,total_matches_minus);
	  }

	  if (print_cycles_p == true) {
	    if (this->use_array_p == false) {
	      /* This sorts by -1 to -readlength, then +readlength to +1 */
	      if ((length = List_length(this->list_matches_byshift)) > 0) {
		printf("(");
		match_array = (Match_T *) List_to_array(this->list_matches_byshift,NULL);
		qsort(match_array,length,sizeof(Match_T),Match_shift_cmp);
		printf("%ld@%d",match_array[0]->count,match_array[0]->shift);
		for (j = 1; j < length; j++) {
		  printf(",%ld@%d",match_array[j]->count,match_array[j]->shift);
		}
		FREE(match_array);
		printf(")");
	      }
	    } else {
	      firstp = true;
	      for (shift = 1; shift < this->n_matches_byshift_minus; shift++) {
		if (this->matches_byshift_minus[shift] > 0) {
		  if (firstp == true) {
		    printf("(%ld@%d",this->matches_byshift_minus[shift],-shift);
		    firstp = false;
		  } else {
		    printf(",%ld@%d",this->matches_byshift_minus[shift],-shift);
		  }
		  this->matches_byshift_minus[shift] = 0; /* clear */
		}
	      }
	      for (shift = this->n_matches_byshift_plus - 1; shift >= 1; shift--) {
		if (this->matches_byshift_plus[shift] > 0) {
		  if (firstp == true) {
		    printf("(%ld@%d",this->matches_byshift_plus[shift],shift);
		    firstp = false;
		  } else {
		    printf(",%ld@%d",this->matches_byshift_plus[shift],shift);
		  }
		  this->matches_byshift_plus[shift] = 0; /* clear */
		}
	      }
	      if (firstp == false) {
		printf(")");
	      }
	    }
	  }

	  if (print_nm_scores_p == true) {
	    if (this->use_array_p == false) {
	      if ((length = List_length(this->list_matches_bynm)) > 0) {
		printf("(");
		match_array = (Match_T *) List_to_array(this->list_matches_bynm,NULL);
		qsort(match_array,length,sizeof(Match_T),Match_nm_cmp);
		printf("%ldNM%d",match_array[0]->count,match_array[0]->nm);
		for (j = 1; j < length; j++) {
		  printf(",%ldNM%d",match_array[j]->count,match_array[j]->nm);
		}
		FREE(match_array);
		printf(")");
	      }
	    } else {
	      firstp = true;
	      for (nm = 0; nm < this->n_matches_bynm; nm++) {
		if (this->matches_bynm[nm] > 0) {
		  if (firstp == true) {
		    printf("(%ldNM%d",this->matches_bynm[nm],nm);
		    firstp = false;
		  } else {
		    printf(",%ldNM%d",this->matches_bynm[nm],nm);
		  }
		  this->matches_bynm[nm] = 0; /* clear */
		}
	      }
	      if (firstp == false) {
		printf(")");
	      }
	    }
	  }

	  if (print_xs_scores_p == true) {
	    if (this->use_array_p == false) {
	      if ((length = List_length(this->list_matches_byxs)) > 0) {
		printf("(");
		match_array = (Match_T *) List_to_array(this->list_matches_byxs,NULL);
		qsort(match_array,length,sizeof(Match_T),Match_xs_cmp);
		printf("%ldXS",match_array[0]->count);
		switch (match_array[0]->xs) {
		case 0: printf("0"); break;
		case 1: printf("+"); break;
		case 2: printf("-"); break;
		default: abort();
		}
		for (j = 1; j < length; j++) {
		  printf(",%ldXS",match_array[j]->count);
		  switch (match_array[j]->xs) {
		  case 0: printf("0"); break;
		  case 1: printf("+"); break;
		  case 2: printf("-"); break;
		  default: abort();
		  }
		}
		FREE(match_array);
		printf(")");
	      }

	    } else {
	      firstp = true;
	      for (xs = 0; xs < this->n_matches_byxs; xs++) {
		if (this->matches_byxs[xs] > 0) {
		  if (firstp == true) {
		    printf("(%ldXS",this->matches_byxs[xs]);
		    firstp = false;
		  } else {
		    printf(",%ldXS",this->matches_byxs[xs]);
		  }
		  switch (xs) {
		  case 0: printf("0"); break;
		  case 1: printf("+"); break;
		  case 2: printf("-"); break;
		  default: abort();
		  }
		  this->matches_byxs[xs] = 0; /* clear */
		}
	      }
	      if (firstp == false) {
		printf(")");
	      }
	    }
          }
	}

	if (this->mismatches_byshift != NULL) {
#ifdef USE_MISMATCHPOOL
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift,this->mismatchpool);
	  unique_mismatches_bynm = make_mismatches_unique(this->mismatches_bynm,this->mismatchpool);
	  unique_mismatches_byxs = make_mismatches_unique(this->mismatches_byxs,this->mismatchpool);
#else
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift);
	  unique_mismatches_bynm = make_mismatches_unique(this->mismatches_bynm);
	  unique_mismatches_byxs = make_mismatches_unique(this->mismatches_byxs);
#endif

	  mm_array = (Mismatch_T *) List_to_array(unique_mismatches_byshift,NULL);
	  qsort(mm_array,List_length(unique_mismatches_byshift),sizeof(Mismatch_T),Mismatch_count_cmp);

	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    mismatch0 = mm_array[i];
	    if (signed_counts_p == false) {
	      printf(" %c%ld",mismatch0->nt,mismatch0->count);
	    } else {
	      printf(" %c%ld|%ld",mismatch0->nt,mismatch0->count_plus,mismatch0->count_minus);
	    }
	    if (print_cycles_p == true) {
	      if ((length = Mismatch_chain_length(mismatch0->next)) > 0) {
		printf("(");
		mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
		for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		  mm_subarray[j] = mismatch;
		}
		qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_shift_cmp);

		printf("%ld@%d",mm_subarray[0]->count,mm_subarray[0]->shift);
		for (j = 1; j < length; j++) {
		  printf(",%ld@%d",mm_subarray[j]->count,mm_subarray[j]->shift);
		}
		FREE(mm_subarray);
		printf(")");
	      }
	    }
	      
	    if (print_nm_scores_p == true) {
	      mismatch0 = find_mismatch_nt(unique_mismatches_bynm,mismatch0->nt);
	      if ((length = Mismatch_chain_length(mismatch0->next)) > 0) {
		printf("(");
		mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
		for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		  mm_subarray[j] = mismatch;
		}
		qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_nm_cmp);

		printf("%ldNM%d",mm_subarray[0]->count,mm_subarray[0]->nm);
		for (j = 1; j < length; j++) {
		  printf(",%ldNM%d",mm_subarray[j]->count,mm_subarray[j]->nm);
		}
		FREE(mm_subarray);
		printf(")");
	      }
	    }

	    if (print_xs_scores_p == true) {
	      mismatch0 = find_mismatch_nt(unique_mismatches_byxs,mismatch0->nt);
	      if ((length = Mismatch_chain_length(mismatch0->next)) > 0) {
		printf("(");
		mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
		for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		  mm_subarray[j] = mismatch;
		}
		qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_xs_cmp);
	      
		printf("%ldXS",mm_subarray[0]->count);
		switch (mm_subarray[0]->xs) {
		case 0: printf("0"); break;
		case 1: printf("+"); break;
		case 2: printf("-"); break;
		default: abort();
		}
		for (j = 1; j < length; j++) {
		  printf(",%ldXS",mm_subarray[j]->count);
		  switch (mm_subarray[j]->xs) {
		  case 0: printf("0"); break;
		  case 1: printf("+"); break;
		  case 2: printf("-"); break;
		  default: abort();
		  }
		}
		FREE(mm_subarray);
		printf(")");
	      }
	    }
	  }

	  FREE(mm_array);

#ifdef USE_MISMATCHPOOL
	  Mismatchpool_reset(this->mismatchpool);
#else
	  for (ptr = unique_mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byshift);

	  for (ptr = unique_mismatches_bynm; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_bynm);

	  for (ptr = unique_mismatches_byxs; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byxs);
#endif
	}

	if (this->delcounts_plus + this->delcounts_minus > 0) {
	  if (signed_counts_p == false) {
	    printf(" _%d",this->delcounts_plus + this->delcounts_minus);
	  } else {
	    printf(" _%d|%d",this->delcounts_plus,this->delcounts_minus);
	  }
	  /* No cycles or other information on deletions */
	}

	if (want_genotypes_p == true) {
#if 0
	  if (this->use_array_p == false) {
	    bestg = compute_probs_list(probs,loglik,this->refnt,this->list_matches_byquality,
				       this->mismatches_byquality,quality_score_adj);
	  } else {
	    bestg = compute_probs_array(probs,loglik,this->refnt,this->matches_byquality,this->n_matches_byquality,
					this->mismatches_byquality,quality_score_adj);
	  }
#endif

	  /* Best genotype */
	  printf("\t");
	  if (0 && bestg < 0) {
	    printf("NN");
	  } else {
	    printf("%s",genotype_string[bestg]);
	  }

	  /* Probabilities */
	  printf("\t");
	  printf("%.3g",1.0-probs[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",1.0-probs[g]);
	  }

	  /* Log-likelihoods */
	  printf("\t");
	  printf("%.3g",loglik[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",loglik[g]);
	  }
	}

	if (map_iit != NULL) {
	  printf("\t");
	  report_plus_genes(/*bytes*/NULL,&ignore,/*tally2*/this,block_tallies,blockstart,prev_block_tallies,
			    prev_blockstart,prev_blockptr,plus_genes,chrpos,genome,chroffset,signed_counts_p,
			    print_cycles_p,print_nm_scores_p,print_xs_scores_p,quality_score_adj,/*output_type*/OUTPUT_BLOCKS);
	  printf("\t");
	  report_minus_genes(/*bytes*/NULL,&ignore,/*tally2*/this,block_tallies,blockstart,prev_block_tallies,
			     prev_blockstart,prev_blockptr,minus_genes,chrpos,genome,chroffset,signed_counts_p,
			     print_cycles_p,print_nm_scores_p,print_xs_scores_p,quality_score_adj,/*output_type*/OUTPUT_BLOCKS);
	  printf("\t");
	  print_gene_info(plus_genes,minus_genes,chrpos);
	}

	printf("\n");
      }

#if 0
      /* Clear all tallies at end of block */
      Tally_clear(this);
#endif
    }

    if (map_iit != NULL) {
      for (p = plus_genes; p != NULL; p = List_next(p)) {
	gene = (Gene_T) List_head(p);
	Gene_free(&gene);
      }
      List_free(&plus_genes);

      for (p = minus_genes; p != NULL; p = List_next(p)) {
	gene = (Gene_T) List_head(p);
	Gene_free(&gene);
      }
      List_free(&minus_genes);
    }

#if 0
    /* Transfer all tallies to prev_block */
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      Tally_clear(this);
    }
#endif

  }
    
  return;
}


static void
tally_block (long int *tally_matches, long int *tally_mismatches,
	     Tally_T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr, 
	     Genome_T genome, char *printchr, Genomicpos_T chroffset, Genomicpos_T chrstart,
	     int quality_score_adj, int min_depth, int variant_strands,
	     bool genomic_diff_p, bool print_noncovered_p) {
  Tally_T this;
  double probs[NGENOTYPES], loglik[NGENOTYPES];
  Genomicpos_T chrpos;
  int blocki, lasti;
  Mismatch_T mismatch;
  List_T ptr;
  long int total;


  lasti = blockptr - blockstart;
  debug3(printf("Printing blocki 0 to %d\n",lasti));

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (print_noncovered_p == true ||
	(pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
			       this->mismatches_byshift,min_depth,variant_strands) == true &&
	 pass_difference_filter_p(probs,loglik,this,genome,
				  printchr,chroffset,chrpos,
				  quality_score_adj,genomic_diff_p) == true)) {
      total += this->nmatches + this->delcounts_plus + this->delcounts_minus;
      for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	mismatch = (Mismatch_T) List_head(ptr);
	total += mismatch->count;
      }
    }
  }
  
  if (total <= 0) {
    /* Total could be 0 if the block is outside chrstart..chrend or if all positions fail variant filter */
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      Tally_clear(this);
    }
    
  } else {
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;

      if (print_noncovered_p == false &&
	  pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
				this->mismatches_byshift,min_depth,variant_strands) == false) {
	/* Skip */

      } else if (print_noncovered_p == false &&
		 pass_difference_filter_p(probs,loglik,this,genome,
					  printchr,chroffset,chrpos,
					  quality_score_adj,genomic_diff_p) == false) {
	/* Skip */

      } else {
#if 0
	bini = (double) (chrpos - mincoord)/binstep;
	if (bini < 0) {
	  bini = 0;
	} else if (bini >= nbins) {
	  bini = nbins - 1;
	}
	binx_matches[bini] += this->nmatches + this->delcounts_plus + this->delcounts_minus;

	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    binx_mismatches[bini] += mismatch->count;
	}
#else
	tally_matches[chrpos - chrstart] += this->nmatches + this->delcounts_plus + this->delcounts_minus;

	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    tally_mismatches[chrpos - chrstart] += mismatch->count;
	}
#endif
      }

      /* Reset */
      Tally_clear(this);
    }
  }
    
  return;
}


static void
iit_block (List_T *intervallist, List_T *labellist, List_T *datalist,
	   Tally_T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr,
	   Tally_T *prev_block_tallies, Genomicpos_T prev_blockstart, Genomicpos_T prev_blockptr,
	   Genome_T genome, char *printchr, Genomicpos_T chroffset, IIT_T map_iit, int quality_score_adj,
	   int min_depth, int variant_strands, bool print_cycles_p, bool print_nm_scores_p,
	   bool print_xs_scores_p, bool print_noncovered_p) {
  Tally_T this;
  char Buffer[100], *label;
  int length, i, j, k;
  List_T plus_genes = NULL, minus_genes = NULL, p;
  Gene_T gene;
  Genomicpos_T chrpos;
  int blocki, lasti;
  Match_T *match_array, match;
  Mismatch_T mismatch, mismatch0, *mm_array, *mm_subarray;
  Insertion_T ins, ins0, *ins_array_byshift, *ins_array_bynm, *ins_array_byxs, *ins_subarray;
  Deletion_T del, del0, *del_array_byshift, *del_array_bynm, *del_array_byxs, *del_subarray;
  List_T unique_mismatches_byshift, unique_mismatches_bynm, unique_mismatches_byxs, ptr;
  List_T unique_insertions_byshift, unique_deletions_byshift, unique_insertions_bynm, unique_deletions_bynm,
    unique_insertions_byxs, unique_deletions_byxs;
  long int total, total_matches_plus, total_matches_minus, total_mismatches_plus, total_mismatches_minus;
  int ninsertions, ndeletions;
  int shift, quality, nm, xs;

  unsigned char *bytearray;
  Ucharlist_T pointers = NULL, bytes = NULL;
  int ignore, nbytes, ntotal;


  lasti = blockptr - blockstart;

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (print_noncovered_p == true ||
	pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
			      this->mismatches_byshift,min_depth,variant_strands) == true) {
      total += this->nmatches + this->delcounts_plus + this->delcounts_minus;
      for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	mismatch = (Mismatch_T) List_head(ptr);
	total += mismatch->count;
      }
    }
  }
  
  if (total <= 0) {
    /* Total could be 0 if the block is outside chrstart..chrend or if all positions fail variant filter */
#if 0
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      Tally_clear(this);
    }
#endif

  } else {
    sprintf(Buffer,"%ld",total);
    label = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
    strcpy(label,Buffer);

    *intervallist = List_push(*intervallist,(void *) Interval_new(blockstart,blockptr-1U,/*sign*/+1));
    *labellist = List_push(*labellist,(void *) label);

    if (map_iit != NULL) {
      genes_get(&plus_genes,&minus_genes,map_iit,printchr,chroffset,genome,
		/*mincoord*/blockstart,/*maxcoord*/blockstart+lasti-1);
    }

    ignore = 0;
    nbytes = 0;
    for (blocki = 0, k = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;
      debug2(printf("data for chrpos %u:\n",chrpos));

      /* Handle insertions at this position */
      debug2(printf("Pointers for insertions:\n"));
      pointers = push_int(&ignore,pointers,nbytes);
      if (this->insertions_byshift != NULL) {
	unique_insertions_byshift = NULL;
	for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byshift,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*cnounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byshift = List_push(unique_insertions_byshift,ins0);
	  }
	}
	ins_array_byshift = (Insertion_T *) List_to_array(unique_insertions_byshift,NULL);
	ninsertions = List_length(unique_insertions_byshift);
	qsort(ins_array_byshift,ninsertions,sizeof(Insertion_T),Insertion_count_cmp);


	unique_insertions_bynm = NULL;
	for (ptr = this->insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_bynm,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*cnounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_bynm = List_push(unique_insertions_bynm,ins0);
	  }
	}
	ins_array_bynm = (Insertion_T *) List_to_array(unique_insertions_bynm,NULL);
	assert(List_length(unique_insertions_bynm) == ninsertions);
	qsort(ins_array_bynm,ninsertions,sizeof(Insertion_T),Insertion_count_cmp);
	/* We hope that ins_array_byshift and ins_array_bynm have the same insertions in parallel */


	unique_insertions_byxs = NULL;
	for (ptr = this->insertions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byxs,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

#if 0
	    ins0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    ins0 = Insertion_new(chrpos,ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*cnounts*/1);
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byxs = List_push(unique_insertions_byxs,ins0);
	  }
	}
	ins_array_byxs = (Insertion_T *) List_to_array(unique_insertions_byxs,NULL);
	assert(List_length(unique_insertions_byxs) == ninsertions);
	qsort(ins_array_byxs,ninsertions,sizeof(Insertion_T),Insertion_count_cmp);



	/* Total number of different insertions at this position */
	debug2(printf("Number of insertions:\n"));
	bytes = push_int(&nbytes,bytes,ninsertions);
	for (i = 0; i < ninsertions; i++) {
	  ins0 = ins_array_byshift[i];
	  /* Counts and segment for insertion i */
	  bytes = push_int(&nbytes,bytes,ins0->count_plus);
	  bytes = push_int(&nbytes,bytes,ins0->count_minus);
	  bytes = push_int(&nbytes,bytes,this->n_fromleft_plus); /* ref count */
	  bytes = push_int(&nbytes,bytes,this->n_fromleft_minus); /* ref count */
	  bytes = push_string(&nbytes,bytes,ins0->segment);

	  /* Cycles for insertion i */
	  length = Insertion_chain_length(ins0->next);
	  ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	  for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
	    ins_subarray[j] = ins;
	  }

	  qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_shift_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->shift);
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->count);
	  }
	  FREE(ins_subarray);

	  /* NM counts for insertion i */
	  ins0 = ins_array_bynm[i];
	  length = Insertion_chain_length(ins0->next);
	  ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	  for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
	    ins_subarray[j] = ins;
	  }

	  qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_nm_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->nm);
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->count);
	  }
	  FREE(ins_subarray);

	  /* NM counts for insertion i */
	  ins0 = ins_array_byxs[i];
	  length = Insertion_chain_length(ins0->next);
	  ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	  for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
	    ins_subarray[j] = ins;
	  }

	  qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_xs_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->xs);
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->count);
	  }
	  FREE(ins_subarray);
	}

	FREE(ins_array_byxs);
	FREE(ins_array_bynm);
	FREE(ins_array_byshift);

	for (ptr = unique_insertions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byxs);

	for (ptr = unique_insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_bynm);

	for (ptr = unique_insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byshift);
      }

      /* Handle deletions at this position */
      debug2(printf("Pointers for deletions:\n"));
      pointers = push_int(&ignore,pointers,nbytes);
      if (this->deletions_byshift != NULL) {
	unique_deletions_byshift = NULL;
	for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byshift,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byshift = List_push(unique_deletions_byshift,del0);
	  }
	}
	del_array_byshift = (Deletion_T *) List_to_array(unique_deletions_byshift,NULL);
	ndeletions = List_length(unique_deletions_byshift);
	qsort(del_array_byshift,ndeletions,sizeof(Deletion_T),Deletion_count_cmp);

	unique_deletions_bynm = NULL;
	for (ptr = this->deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_bynm,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_bynm = List_push(unique_deletions_bynm,del0);
	  }
	}
	del_array_bynm = (Deletion_T *) List_to_array(unique_deletions_bynm,NULL);
	assert(List_length(unique_deletions_bynm) == ndeletions);
	qsort(del_array_bynm,ndeletions,sizeof(Deletion_T),Deletion_count_cmp);
	/* We hope that del_array_byshift and del_array_bynm have the same deletions in parallel */

	unique_deletions_byxs = NULL;
	for (ptr = this->deletions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byxs,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

#if 0
	    del0->shift += 1; /* Used here as nshifts.  Not necessary. */
#endif
	  } else {
	    del0 = Deletion_new(chrpos,del->segment,del->mlength,/*shift, used here as nshifts*/1,/*nm*/0,/*xs*/0,/*ncounts*/1);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byxs = List_push(unique_deletions_byxs,del0);
	  }
	}
	del_array_byxs = (Deletion_T *) List_to_array(unique_deletions_byxs,NULL);
	assert(List_length(unique_deletions_byxs) == ndeletions);
	qsort(del_array_byxs,ndeletions,sizeof(Deletion_T),Deletion_count_cmp);


	/* Total number of different deletions at this position */
	debug2(printf("Number of deletions:\n"));
	bytes = push_int(&nbytes,bytes,ndeletions);
	for (i = 0; i < ndeletions; i++) {
	  del0 = del_array_byshift[i];
	  /* Counts and segment for deletion i */
	  debug2(printf("plus and minus counts for deletion\n"));
	  bytes = push_int(&nbytes,bytes,del0->count_plus);
	  bytes = push_int(&nbytes,bytes,del0->count_minus);
	  debug2(printf("plus and minus counts for reference\n"));
	  bytes = push_int(&nbytes,bytes,this->n_fromleft_plus); /* ref count */
	  bytes = push_int(&nbytes,bytes,this->n_fromleft_minus); /* ref count */
	  debug2(printf("Deletion string:\n"));
	  bytes = push_string(&nbytes,bytes,del0->segment);

	  /* Cycles for deletion i */
	  length = Deletion_chain_length(del0->next);
	  del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	  for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
	    del_subarray[j] = del;
	  }

	  qsort(del_subarray,length,sizeof(Deletion_T),Deletion_shift_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->shift);
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->count);
	  }
	  FREE(del_subarray);

	  /* NM counts for deletion i */
	  del0 = del_array_bynm[i];
	  length = Deletion_chain_length(del0->next);
	  del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	  for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
	    del_subarray[j] = del;
	  }

	  qsort(del_subarray,length,sizeof(Deletion_T),Deletion_nm_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->nm);
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->count);
	  }
	  FREE(del_subarray);

	  /* XS counts for deletion i */
	  del0 = del_array_byxs[i];
	  length = Deletion_chain_length(del0->next);
	  del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	  for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
	    del_subarray[j] = del;
	  }

	  qsort(del_subarray,length,sizeof(Deletion_T),Deletion_xs_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->xs);
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->count);
	  }
	  FREE(del_subarray);
	}

	FREE(del_array_byxs);
	FREE(del_array_bynm);
	FREE(del_array_byshift);

	for (ptr = unique_deletions_byxs; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byxs);

	for (ptr = unique_deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_bynm);

	for (ptr = unique_deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byshift);
      }
      

      /* Handle allele counts at this position */
      debug2(printf("Pointers for allele counts:\n"));
      pointers = push_int(&ignore,pointers,nbytes);
      if (print_noncovered_p == true ||
	  pass_variant_filter_p(this->nmatches,this->delcounts_plus,this->delcounts_minus,
				this->mismatches_byshift,min_depth,variant_strands) == true) {
	total_matches_plus = total_matches_minus = 0;
	if (this->use_array_p == false) {
	  for (ptr = this->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    match = (Match_T) List_head(ptr);
	    if (match->shift > 0) {
	      total_matches_plus += match->count;
	    } else {
	      total_matches_minus += match->count;
	    }
	  }
	} else {
	  for (shift = 0; shift < this->n_matches_byshift_plus; shift++) {
	    total_matches_plus += this->matches_byshift_plus[shift];
	  }
	  for (shift = 0; shift < this->n_matches_byshift_minus; shift++) {
	    total_matches_minus += this->matches_byshift_minus[shift];
	  }
	}
	assert(this->nmatches == total_matches_plus + total_matches_minus);

	total_mismatches_plus = total_mismatches_minus = 0;
	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  if (mismatch->shift > 0) {
	    total_mismatches_plus += mismatch->count;
	  } else {
	    total_mismatches_minus += mismatch->count;
	  }
	}

	/* Total signed counts at this position */
	debug2(printf("Total signed counts:\n"));
	bytes = push_int(&nbytes,bytes,total_matches_plus + this->delcounts_plus + total_mismatches_plus);
	bytes = push_int(&nbytes,bytes,total_matches_minus + this->delcounts_minus + total_mismatches_minus);
	
	/* Reference nucleotide and counts */
	if (this->nmatches == 0) {
	  bytes = push_char(&nbytes,bytes,Genome_get_char(genome,chroffset+chrpos-1U));
	  bytes = push_int(&nbytes,bytes,/*plus*/0);
	  bytes = push_int(&nbytes,bytes,/*minus*/0);
	} else {
	  bytes = push_char(&nbytes,bytes,this->refnt);
	  bytes = push_int(&nbytes,bytes,total_matches_plus);
	  bytes = push_int(&nbytes,bytes,total_matches_minus);
	}

	/* Alternate nucleotide and counts */
	if (this->mismatches_byshift != NULL) {
#ifdef USE_MISMATCHPOOL
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift,this->mismatchpool);
	  unique_mismatches_bynm = make_mismatches_unique(this->mismatches_bynm,this->mismatchpool);
	  unique_mismatches_byxs = make_mismatches_unique(this->mismatches_byxs,this->mismatchpool);
#else
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift);
	  unique_mismatches_bynm = make_mismatches_unique(this->mismatches_bynm);
	  unique_mismatches_byxs = make_mismatches_unique(this->mismatches_byxs);
#endif

	  mm_array = (Mismatch_T *) List_to_array(unique_mismatches_byshift,NULL);
	  qsort(mm_array,List_length(unique_mismatches_byshift),sizeof(Mismatch_T),Mismatch_count_cmp);

	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    mismatch0 = mm_array[i];
	    bytes = push_char(&nbytes,bytes,mismatch0->nt);
	    bytes = push_int(&nbytes,bytes,mismatch0->count_plus);
	    bytes = push_int(&nbytes,bytes,mismatch0->count_minus);
	  }
	}

	if (this->delcounts_plus + this->delcounts_minus > 0) {
	  bytes = push_char(&nbytes,bytes,'_');
	  bytes = push_int(&nbytes,bytes,this->delcounts_plus);
	  bytes = push_int(&nbytes,bytes,this->delcounts_minus);
	}

	/* Terminates nucleotides */
	bytes = push_char(&nbytes,bytes,'\0');
	
	/* Cycles, nm, and optionally xs for reference */
	if (this->nmatches + this->delcounts_plus + this->delcounts_minus > 0) {
	  if (this->use_array_p == false) {
	    length = List_length(this->list_matches_byshift);
	    match_array = (Match_T *) List_to_array(this->list_matches_byshift,NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_shift_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,match_array[j]->shift);
	      bytes = push_int(&nbytes,bytes,match_array[j]->count);
	    }
	    FREE(match_array);
	  } else {
	    length = 0;
	    for (shift = 1; shift < this->n_matches_byshift_minus; shift++) {
	      if (this->matches_byshift_minus[shift] > 0) {
		length++;
	      }
	    }
	    for (shift = this->n_matches_byshift_plus - 1; shift >= 1; shift--) {
	      if (this->matches_byshift_plus[shift] > 0) {
		length++;
	      }
	    }
	    bytes = push_int(&nbytes,bytes,length);
	    for (shift = 1; shift < this->n_matches_byshift_minus; shift++) {
	      if (this->matches_byshift_minus[shift] > 0) {
		bytes = push_int(&nbytes,bytes,-shift);
		bytes = push_int(&nbytes,bytes,this->matches_byshift_minus[shift]);
		this->matches_byshift_minus[shift] = 0; /* clear */
	      }
	    }
	    for (shift = this->n_matches_byshift_plus - 1; shift >= 1; shift--) {
	      if (this->matches_byshift_plus[shift] > 0) {
		bytes = push_int(&nbytes,bytes,+shift);
		bytes = push_int(&nbytes,bytes,this->matches_byshift_plus[shift]);
		this->matches_byshift_plus[shift] = 0; /* clear */
	      }
	    }
	  }

	  if (this->use_array_p == false) {
	    length = List_length(this->list_matches_bynm);
	    match_array = (Match_T *) List_to_array(this->list_matches_bynm,NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_nm_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,match_array[j]->nm);
	      bytes = push_int(&nbytes,bytes,match_array[j]->count);
	    }
	    FREE(match_array);
	  } else {
	    length = 0;
	    for (nm = 0; nm < this->n_matches_bynm; nm++) {
	      if (this->matches_bynm[nm] > 0) {
		length++;
	      }
	    }
	    bytes = push_int(&nbytes,bytes,length);
	    for (nm = 0; nm < this->n_matches_bynm; nm++) {
	      if (this->matches_bynm[nm] > 0) {
		bytes = push_int(&nbytes,bytes,nm);
		bytes = push_int(&nbytes,bytes,this->matches_bynm[nm]);
		this->matches_bynm[nm] = 0; /* clear */
	      }
	    }
	  }

	  if (print_xs_scores_p == true) {
	    if (this->use_array_p == false) {
	      length = List_length(this->list_matches_byxs);
	      match_array = (Match_T *) List_to_array(this->list_matches_byxs,NULL);
	      qsort(match_array,length,sizeof(Match_T),Match_xs_cmp);
	      bytes = push_int(&nbytes,bytes,length);
	      for (j = 0; j < length; j++) {
		bytes = push_int(&nbytes,bytes,match_array[j]->xs);
		bytes = push_int(&nbytes,bytes,match_array[j]->count);
	      }
	      FREE(match_array);
	    } else {
	      length = 0;
	      for (xs = 0; xs < this->n_matches_byxs; xs++) {
		if (this->matches_byxs[xs] > 0) {
		  length++;
		}
	      }
	      bytes = push_int(&nbytes,bytes,length);
	      for (xs = 0; xs < this->n_matches_byxs; xs++) {
		if (this->matches_byxs[xs] > 0) {
		  bytes = push_int(&nbytes,bytes,xs);
		  bytes = push_int(&nbytes,bytes,this->matches_byxs[xs]);
		  this->matches_byxs[xs] = 0; /* clear */
		}
	      }
	    }
	  }
	}
	    
	if (this->mismatches_byshift != NULL) {
	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    /* Cycles for alternate allele i */
	    mismatch0 = mm_array[i];
	    length = Mismatch_chain_length(mismatch0->next);
	    mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	    for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
	      mm_subarray[j] = mismatch;
	    }

	    qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_shift_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->shift);
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->count);
	    }
	    FREE(mm_subarray);


	    /* NM scores for alternate allele i */
	    mismatch0 = find_mismatch_nt(unique_mismatches_bynm,mismatch0->nt);
	    length = Mismatch_chain_length(mismatch0->next);
	    mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	    for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
	      mm_subarray[j] = mismatch;
	    }
	    qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_nm_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->nm);
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->count);
	    }
	    FREE(mm_subarray);

	    if (print_xs_scores_p == true) {
	      /* XS scores for alternate allele i */
	      mismatch0 = find_mismatch_nt(unique_mismatches_byxs,mismatch0->nt);
	      length = Mismatch_chain_length(mismatch0->next);
	      mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		mm_subarray[j] = mismatch;
	      }
	      qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_xs_cmp);
	      bytes = push_int(&nbytes,bytes,length);
	      for (j = 0; j < length; j++) {
		bytes = push_int(&nbytes,bytes,mm_subarray[j]->xs);
		bytes = push_int(&nbytes,bytes,mm_subarray[j]->count);
	      }
	      FREE(mm_subarray);
	    }
	  }
	  FREE(mm_array);

#ifdef USE_MISMATCHPOOL
	  Mismatchpool_reset(this->mismatchpool);
#else
	  for (ptr = unique_mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byshift);

	  for (ptr = unique_mismatches_bynm; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_bynm);

	  for (ptr = unique_mismatches_byxs; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byxs);
#endif
	}

	if (this->delcounts_plus + this->delcounts_minus > 0) {
	  /* No cycle or other information should be output here */
	}

      }


      debug2(printf("Pointers for plus-strand amino acid counts:\n"));
      pointers = push_int(&ignore,pointers,nbytes);
      if (map_iit != NULL) {
	bytes = report_plus_genes(bytes,&nbytes,/*tally2*/this,block_tallies,blockstart,prev_block_tallies,
				  prev_blockstart,prev_blockptr,plus_genes,chrpos,genome,chroffset,/*signed_counts_p*/true,
				  /*print_cycles_p*/true,/*print_nm_scores_p*/true,print_xs_scores_p,
				  quality_score_adj,/*output_type*/OUTPUT_IIT);
      }

      debug2(printf("Pointers for minus-strand amino acid counts:\n"));
      pointers = push_int(&ignore,pointers,nbytes);
      if (map_iit != NULL) {
	bytes = report_minus_genes(bytes,&nbytes,/*tally2*/this,block_tallies,blockstart,prev_block_tallies,
				   prev_blockstart,prev_blockptr,minus_genes,chrpos,genome,chroffset,/*signed_counts_p*/true,
				   /*print_cycles_p*/true,/*print_nm_scores_p*/true,print_xs_scores_p,
				   quality_score_adj,/*output_type*/OUTPUT_IIT);
      }

#if 0
      /* Clear all tallies at end of block */
      Tally_clear(this);
#endif
    }

    if (map_iit != NULL) {
      for (p = plus_genes; p != NULL; p = List_next(p)) {
	gene = (Gene_T) List_head(p);
	Gene_free(&gene);
      }
      List_free(&plus_genes);

      for (p = minus_genes; p != NULL; p = List_next(p)) {
	gene = (Gene_T) List_head(p);
	Gene_free(&gene);
      }
      List_free(&minus_genes);
    }

#if 0
    /* Transfer all tallies to prev_block */
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      Tally_clear(this);
    }
#endif

    /* Final pointer for block */
    pointers = push_int(&ignore,pointers,nbytes);

    /* Reverse lists before appending */
    pointers = Ucharlist_append(Ucharlist_reverse(pointers),Ucharlist_reverse(bytes));
    bytearray = Ucharlist_to_array(&ntotal,pointers);
    Ucharlist_free(&pointers);
    
    *datalist = List_push(*datalist,(void *) bytearray);
  }

  return;
}


/* A modified version is in indelfix.c */
static Genomicpos_T
transfer_position (long int *grand_total, Tally_T *alloc_tallies, Tally_T *block_tallies,
		   Genomicpos_T *exonstart, Genomicpos_T lastpos,
		   Genomicpos_T *blockptr, Genomicpos_T *blockstart, Genomicpos_T *blockend,

		   Tally_T *prev_block_tallies, Genomicpos_T *prev_blockptr, Genomicpos_T *prev_blockstart, Genomicpos_T *prev_blockend,

		   long int *tally_matches, long int *tally_mismatches, 
		   List_T *intervallist, List_T *labellist, List_T *datalist,
		   Genomicpos_T chrpos, int alloci,
		   Genome_T genome, char *printchr, Genomicpos_T chroffset,
		   Genomicpos_T chrstart, Genomicpos_T chrend, IIT_T map_iit,
		   Tally_outputtype_T output_type, bool blockp, int blocksize,
		   int quality_score_adj, int min_depth, int variant_strands,
		   bool genomic_diff_p, bool signed_counts_p, bool print_totals_p, bool print_cycles_p,
		   bool print_nm_scores_p, bool print_xs_scores_p, bool want_genotypes_p,
		   bool readlevel_p, bool print_noncovered_p) {
  int blocki;

  if (chrpos < chrstart) {
    debug0(printf("    excluding chrpos %u < chrstart %u\n",chrpos,chrstart));
    Tally_clear(alloc_tallies[alloci]);
  } else if (chrpos > chrend) {
    debug0(printf("    excluding chrpos %u > chrend %u\n",chrpos,chrend));
    Tally_clear(alloc_tallies[alloci]);
  } else {
    if (chrpos >= *blockend) {
      debug0(printf("    chrpos %u >= blockend %u\n",chrpos,*blockend));
      debug0(printf("      print block from blockstart %u to blockptr %u\n",*blockstart,*blockptr));

      if (output_type == OUTPUT_RUNLENGTHS) {
	lastpos = print_runlength(block_tallies,&(*exonstart),lastpos,*blockstart,*blockptr,printchr);
      } else if (output_type == OUTPUT_TALLY) {
	tally_block(tally_matches,tally_mismatches,
		    block_tallies,*blockstart,*blockptr,genome,printchr,chroffset,chrstart,
		    quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		    print_noncovered_p);
      } else if (output_type == OUTPUT_IIT) {
	iit_block(&(*intervallist),&(*labellist),&(*datalist),
		  block_tallies,*blockstart,*blockptr,
		  prev_block_tallies,*prev_blockstart,*prev_blockptr,
		  genome,printchr,chroffset,map_iit,
		  quality_score_adj,min_depth,variant_strands,
	          print_cycles_p,print_nm_scores_p,print_xs_scores_p,print_noncovered_p);

	for (blocki = 0; blocki < *prev_blockptr - *prev_blockstart; blocki++) {
	  Tally_clear(prev_block_tallies[blocki]);
	}
	for (blocki = 0; blocki < *blockptr - *blockstart; blocki++) {
	  Tally_transfer(&(prev_block_tallies[blocki]),&(block_tallies[blocki]));
	}
	*prev_blockstart = *blockstart;
	*prev_blockptr = *blockptr;

      } else if (output_type == OUTPUT_TOTAL) {
	*grand_total += block_total(block_tallies,*blockstart,*blockptr);
      } else {
	if (print_noncovered_p == true) {
	  debug0(printf("Printing zeroes from lastpos %u to blockstart %u\n",lastpos,*blockstart));
	  print_zeroes(lastpos,*blockstart,printchr,blocksize,genome,chroffset,blockp);
	}
	print_block(block_tallies,*blockstart,*blockptr,
		    prev_block_tallies,*prev_blockstart,*prev_blockptr,
		    genome,printchr,chroffset,map_iit,
		    blockp,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		    signed_counts_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
	            want_genotypes_p,readlevel_p,print_noncovered_p);

	for (blocki = 0; blocki < *prev_blockptr - *prev_blockstart; blocki++) {
	  Tally_clear(prev_block_tallies[blocki]);
	}
	for (blocki = 0; blocki < *blockptr - *blockstart; blocki++) {
	  Tally_transfer(&(prev_block_tallies[blocki]),&(block_tallies[blocki]));
	}
	*prev_blockstart = *blockstart;
	*prev_blockptr = *blockptr;

	if (*blockptr > lastpos) {
	  lastpos = *blockptr;
	}
      }

      debug0(printf("      set blockstart to chrpos, blockend to chrpos + blocksize\n"));
      *blockstart = chrpos;
      *blockend = chrpos + blocksize;
    }

    blocki = chrpos - (*blockstart);
#if 0
    debug0(printf("      transfer position from alloci %d to blocki %d\n",alloci,blocki));
#endif

    Tally_transfer(&(block_tallies[blocki]),&(alloc_tallies[alloci]));
    block_tallies[blocki+1]->n_fromleft_plus = alloc_tallies[alloci+1]->n_fromleft_plus;
    block_tallies[blocki+1]->n_fromleft_minus = alloc_tallies[alloci+1]->n_fromleft_minus;

    /* Points just beyond maximum chrpos observed */
    if (chrpos + 1U > *blockptr) {
      *blockptr = chrpos + 1U;
#if 0
      debug0(printf("    advancing blockptr to %u\n",*blockptr));
#endif
    }
  }

  return lastpos;
}


static void
transfer_position_lh (Tally_T *alloc_tallies_low, Tally_T *alloc_tallies_high,
		      Tally_T *block_tallies_low, Tally_T *block_tallies_high,
		      Genomicpos_T *blockptr, Genomicpos_T *blockstart, Genomicpos_T *blockend,
		      long int *tally_matches_low, long int *tally_mismatches_low,
		      long int *tally_matches_high, long int *tally_mismatches_high,
		      Genomicpos_T chrpos, int alloci,
		      Genome_T genome, char *printchr, Genomicpos_T chroffset,
		      Genomicpos_T chrstart, Genomicpos_T chrend, IIT_T map_iit,
		      int blocksize, int quality_score_adj, int min_depth, int variant_strands,
		      bool genomic_diff_p, bool print_noncovered_p) {
  /* Genomicpos_T lastpos; */
  int blocki;

  if (chrpos < chrstart) {
    debug0(printf("    excluding chrpos %u < chrstart %u\n",chrpos,chrstart));
    Tally_clear(alloc_tallies_low[alloci]);
    Tally_clear(alloc_tallies_high[alloci]);
  } else if (chrpos > chrend) {
    debug0(printf("    excluding chrpos %u > chrend %u\n",chrpos,chrend));
    Tally_clear(alloc_tallies_low[alloci]);
    Tally_clear(alloc_tallies_high[alloci]);
  } else {
    if (chrpos >= *blockend) {
      debug0(printf("    chrpos %u >= blockend %u\n",chrpos,*blockend));
      debug0(printf("      print block from blockstart %u to blockptr %u\n",*blockstart,*blockptr));
      
      tally_block(tally_matches_low,tally_mismatches_low,
		  block_tallies_low,*blockstart,*blockptr,genome,printchr,chroffset,chrstart,
		  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		  print_noncovered_p);

      tally_block(tally_matches_high,tally_mismatches_high,
		  block_tallies_high,*blockstart,*blockptr,genome,printchr,chroffset,chrstart,
		  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		  print_noncovered_p);

      debug0(printf("      set blockstart to chrpos, blockend to chrpos + blocksize\n"));
      *blockstart = chrpos;
      *blockend = chrpos + blocksize;
    }

    blocki = chrpos - (*blockstart);
#if 0
    debug0(printf("      transfer position from alloci %d to blocki %d\n",alloci,blocki));
#endif

    Tally_transfer(&(block_tallies_low[blocki]),&(alloc_tallies_low[alloci]));
    Tally_transfer(&(block_tallies_high[blocki]),&(alloc_tallies_high[alloci]));
    block_tallies_low[blocki+1]->n_fromleft_plus = alloc_tallies_low[alloci+1]->n_fromleft_plus;
    block_tallies_low[blocki+1]->n_fromleft_minus = alloc_tallies_low[alloci+1]->n_fromleft_minus;
    block_tallies_high[blocki+1]->n_fromleft_plus = alloc_tallies_high[alloci+1]->n_fromleft_plus;
    block_tallies_high[blocki+1]->n_fromleft_minus = alloc_tallies_high[alloci+1]->n_fromleft_minus;

    /* Points just beyond maximum chrpos observed */
    if (chrpos + 1U > *blockptr) {
      *blockptr = chrpos + 1U;
#if 0
      debug0(printf("    advancing blockptr to %u\n",*blockptr));
#endif
    }
  }

  return;
}




static void
revise_position (char querynt, char genomicnt, int nm, int xs, int signed_shift,
		 Tally_T this, Tally_T right, bool ignore_query_Ns_p, bool readlevel_p,
		 unsigned int linei, int ncounts) {
  Match_T match;
  Mismatch_T mismatch;
  List_T ptr;
  int *newarray, oldsize;
  int max_shift_plus, max_shift_minus;

  if (right != NULL) {
    if (signed_shift > 0) {
      right->n_fromleft_plus += ncounts;
    } else {
      right->n_fromleft_minus += ncounts;
    }
  }

  if (readlevel_p == true || totals_only_p == true) {
    /* Don't store any details.  Also, nmatches captures both matches and mismatches */
    this->nmatches += ncounts;

  } else if (toupper(querynt) == toupper(genomicnt)) {
    /* Count matches, even if query quality score was low */
    this->nmatches += ncounts;
    if (this->nmatches == ARRAY_THRESHOLD) {
      /* Convert list to array */
      /* Find maximum shifts */
      max_shift_plus = max_shift_minus = 0;
      for (ptr = this->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	match = (Match_T) List_head(ptr);
	if (match->shift > 0) {
	  if (match->shift > max_shift_plus) {
	    max_shift_plus = match->shift;
	  }
	} else {
	  if (-match->shift > max_shift_minus) {
	    max_shift_minus = -match->shift;
	  }
	}
      }
      
      if (max_shift_plus < this->n_matches_byshift_plus) {
	/* Clear existing array.  Not necessary here if Tally_clear is doing the memset */
	/* memset(this->matches_byshift_plus,0,this->n_matches_byshift_plus * sizeof(int)); */
      } else {
	/* Resize array */
	oldsize = this->n_matches_byshift_plus;
	while (max_shift_plus >= this->n_matches_byshift_plus) {
	  this->n_matches_byshift_plus *= 2;
	}
	newarray = (int *) CALLOC(this->n_matches_byshift_plus,sizeof(int));
	/* memcpy(newarray,this->matches_byshift_plus,oldsize * sizeof(int)); -- Should be empty */
	FREE(this->matches_byshift_plus);
	this->matches_byshift_plus = newarray;
      }

      if (max_shift_minus < this->n_matches_byshift_minus) {
	/* Clear existing array.  Not necessary here if Tally_clear is doing the memset */
	/* memset(this->matches_byshift_minus,0,this->n_matches_byshift_minus * sizeof(int)); */
      } else {
	/* Resize array */
	oldsize = this->n_matches_byshift_minus;
	while (max_shift_minus >= this->n_matches_byshift_minus) {
	  this->n_matches_byshift_minus *= 2;
	}
	newarray = (int *) CALLOC(this->n_matches_byshift_minus,sizeof(int));
	/* memcpy(newarray,this->matches_byshift_minus,oldsize * sizeof(int)); -- Should be empty */
	FREE(this->matches_byshift_minus);
	this->matches_byshift_minus = newarray;
      }


      for (ptr = this->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	match = (Match_T) List_head(ptr);
	if (match->shift > 0) {
	  this->matches_byshift_plus[match->shift] = match->count;
	} else {
	  this->matches_byshift_minus[-match->shift] = match->count;
	}
#ifndef USE_MATCHPOOL
	Match_free(&match);
#endif
      }
#ifndef USE_MATCHPOOL
      List_free(&(this->list_matches_byshift));
      this->list_matches_byshift = (List_T) NULL;
#endif

      for (ptr = this->list_matches_bynm; ptr != NULL; ptr = List_next(ptr)) {
	match = (Match_T) List_head(ptr);
	this->matches_bynm[match->nm] = match->count;
#ifndef USE_MATCHPOOL
	Match_free(&match);
#endif
      }
#ifndef USE_MATCHPOOL
      List_free(&(this->list_matches_bynm));
      this->list_matches_bynm = (List_T) NULL;
#endif

      for (ptr = this->list_matches_byxs; ptr != NULL; ptr = List_next(ptr)) {
	match = (Match_T) List_head(ptr);
	this->matches_byxs[match->xs] = match->count;
#ifndef USE_MATCHPOOL
	Match_free(&match);
#endif
      }
#ifndef USE_MATCHPOOL
      List_free(&(this->list_matches_byxs));
      this->list_matches_byxs = (List_T) NULL;
#endif

      this->use_array_p = true;
    }

    if (this->use_array_p == false) {
      if ((match = find_match_byshift(this->list_matches_byshift,signed_shift)) != NULL) {
	match->count += ncounts;
      } else {
#ifdef USE_MATCHPOOL
        this->list_matches_byshift = Matchpool_push(this->list_matches_byshift,this->matchpool,signed_shift,nm,xs,ncounts);
#else
	this->list_matches_byshift = List_push(this->list_matches_byshift,(void *) Match_new(signed_shift,nm,xs,ncounts));
#endif
      }

      if ((match = find_match_bynm(this->list_matches_bynm,nm)) != NULL) {
	match->count += ncounts;
      } else {
#ifdef USE_MATCHPOOL
        this->list_matches_bynm = Matchpool_push(this->list_matches_bynm,this->matchpool,signed_shift,nm,xs,ncounts);
#else
        this->list_matches_bynm = List_push(this->list_matches_bynm,(void *) Match_new(signed_shift,nm,xs,ncounts));
#endif
      }

      if ((match = find_match_byxs(this->list_matches_byxs,xs)) != NULL) {
	match->count += ncounts;
      } else {
#ifdef USE_MATCHPOOL
        this->list_matches_byxs = Matchpool_push(this->list_matches_byxs,this->matchpool,signed_shift,nm,xs,ncounts);
#else
        this->list_matches_byxs = List_push(this->list_matches_byxs,(void *) Match_new(signed_shift,nm,xs,ncounts));
#endif
      }

    } else {
      if (signed_shift >= 0) {
	if (signed_shift >= this->n_matches_byshift_plus) {
	  /* Resize array */
	  oldsize = this->n_matches_byshift_plus;
	  while (signed_shift >= this->n_matches_byshift_plus) {
	    this->n_matches_byshift_plus *= 2;
	  }
	  newarray = (int *) CALLOC(this->n_matches_byshift_plus,sizeof(int));
	  memcpy(newarray,this->matches_byshift_plus,oldsize * sizeof(int));
	  FREE(this->matches_byshift_plus);
	  this->matches_byshift_plus = newarray;
	}
	this->matches_byshift_plus[signed_shift] += ncounts;
      } else {
	if (-signed_shift >= this->n_matches_byshift_minus) {
	  /* Resize array */
	  oldsize = this->n_matches_byshift_minus;
	  while (-signed_shift >= this->n_matches_byshift_minus) {
	    this->n_matches_byshift_minus *= 2;
	  }
	  newarray = (int *) CALLOC(this->n_matches_byshift_minus,sizeof(int));
	  memcpy(newarray,this->matches_byshift_minus,oldsize * sizeof(int));
	  FREE(this->matches_byshift_minus);
	  this->matches_byshift_minus = newarray;
	}
	this->matches_byshift_minus[-signed_shift] += ncounts;
      }

      this->matches_bynm[nm] += ncounts;
      this->matches_byxs[xs] += ncounts;
    }

  } else if (toupper(querynt) == 'N' && ignore_query_Ns_p == true) {
    /* Skip query N's */

  } else {

    if ((mismatch = find_mismatch_byshift(this->mismatches_byshift,toupper(querynt),signed_shift)) != NULL) {
      mismatch->count += ncounts;
    } else {
#ifdef USE_MISMATCHPOOL
      this->mismatches_byshift = Mismatchpool_push(this->mismatches_byshift,this->mismatchpool,
                                           	   toupper(querynt),signed_shift,nm,xs,ncounts);
#else
      this->mismatches_byshift = List_push(this->mismatches_byshift,(void *) Mismatch_new(toupper(querynt),signed_shift,nm,xs,ncounts));
#endif
    }

    if ((mismatch = find_mismatch_bynm(this->mismatches_bynm,toupper(querynt),nm)) != NULL) {
      mismatch->count += ncounts;
    } else {
#ifdef USE_MISMATCHPOOL
      this->mismatches_bynm = Mismatchpool_push(this->mismatches_bynm,this->mismatchpool,
						toupper(querynt),signed_shift,nm,xs,ncounts);
#else
      this->mismatches_bynm = List_push(this->mismatches_bynm,(void *) Mismatch_new(toupper(querynt),signed_shift,nm,xs,ncounts));
#endif
    }

    if ((mismatch = find_mismatch_byxs(this->mismatches_byxs,toupper(querynt),xs)) != NULL) {
      mismatch->count += ncounts;
    } else {
#ifdef USE_MISMATCHPOOL
      this->mismatches_byxs = Mismatchpool_push(this->mismatches_byxs,this->mismatchpool,
						toupper(querynt),signed_shift,nm,xs,ncounts);
#else
      this->mismatches_byxs = List_push(this->mismatches_byxs,(void *) Mismatch_new(toupper(querynt),signed_shift,nm,xs,ncounts));
#endif
    }
  }

  /* Only needed if map_iit != NULL */
  this->readevidence = List_push(this->readevidence,Readevid_new(linei,toupper(querynt),signed_shift,nm,xs));

  return;
}


static void
revise_insertion (Genomicpos_T chrpos, char *query_insert, int mlength, int signed_shift, int nm, int xs,
		  Tally_T this, int ncounts) {
  Insertion_T ins;

  if ((ins = find_insertion_byshift(this->insertions_byshift,query_insert,mlength,signed_shift)) != NULL) {
    ins->count += ncounts;
  } else {
    this->insertions_byshift = List_push(this->insertions_byshift,(void *) Insertion_new(chrpos,query_insert,mlength,signed_shift,nm,xs,ncounts));
  }

  if ((ins = find_insertion_bynm(this->insertions_bynm,query_insert,mlength,nm)) != NULL) {
    ins->count += ncounts;
  } else {
    this->insertions_bynm = List_push(this->insertions_bynm,(void *) Insertion_new(chrpos,query_insert,mlength,signed_shift,nm,xs,ncounts));
  }

  if ((ins = find_insertion_byxs(this->insertions_byxs,query_insert,mlength,xs)) != NULL) {
    ins->count += ncounts;
  } else {
    this->insertions_byxs = List_push(this->insertions_byxs,(void *) Insertion_new(chrpos,query_insert,mlength,signed_shift,nm,xs,ncounts));
  }

  return;
}

static void
revise_deletion (Genomicpos_T chrpos, char *deletion, int mlength, int signed_shift, int nm, int xs, Tally_T this,
		 int ncounts) {
  Deletion_T del;

  if ((del = find_deletion_byshift(this->deletions_byshift,deletion,mlength,signed_shift)) != NULL) {
    del->count += ncounts;
  } else {
    this->deletions_byshift = List_push(this->deletions_byshift,(void *) Deletion_new(chrpos,deletion,mlength,signed_shift,nm,xs,ncounts));
  }

  if ((del = find_deletion_bynm(this->deletions_bynm,deletion,mlength,nm)) != NULL) {
    del->count += ncounts;
  } else {
    this->deletions_bynm = List_push(this->deletions_bynm,(void *) Deletion_new(chrpos,deletion,mlength,signed_shift,nm,xs,ncounts));
  }

  if ((del = find_deletion_byxs(this->deletions_byxs,deletion,mlength,xs)) != NULL) {
    del->count += ncounts;
  } else {
    this->deletions_byxs = List_push(this->deletions_byxs,(void *) Deletion_new(chrpos,deletion,mlength,signed_shift,nm,xs,ncounts));
  }

  return;
}



static double
average (char *quality_scores, int n) {
  double total = 0.0;
  int i;

  for (i = 0; i < n; i++) {
    total += (double) quality_scores[i];
  }
  return total/(double) n;
}


static void
revise_read (Tally_T *alloc_tallies, Genomicpos_T chrstart, Genomicpos_T chrend, Genomicpos_T chrpos_low,
	     unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
	     char *shortread, int nm, char splice_strand, bool terminalp, Genomicpos_T alloc_low,
	     Genome_T genome, Genomicpos_T chroffset, bool ignore_query_Ns_p,
	     bool print_indels_p, bool readlevel_p, int max_softclip, unsigned int linei,
	     __m128i *counts, int nreps) {
  Tally_T this, right;
  int alloci;
  int shift, signed_shift;
  char strand;
  Genomicpos_T pos;
  int r, i;
  char *genomic = NULL, *p, *q;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int ncounts;
  int type;
  int quality_score;
  int xs;
  short block[8];
  __m128i foo;


  debug1(printf("Revising read at %u\n",chrpos_low));

  if (splice_strand == '+') {
    xs = 1;
  } else if (splice_strand == '-') {
    xs = 2;
  } else {
    xs = 0;
  }

  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
  }

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */

  if (terminalp == false) {
    /* Don't handle soft clip for terminal alignments */
    max_softclip = 0;
  }

  p = shortread;
  r = 0;
  if (counts) _mm_store_si128((__m128i *) block,counts[0]);
  if (max_softclip > 0 && types != NULL) {
    if (Intlist_head(types) == 'S') {
      /* Revise pos, so we handle the initial soft clip */
      mlength = Uintlist_head(npositions);
      if (pos < mlength) {
	/* Make sure initial soft clip does not extend past beginning of chromosome */
	pos = 0U;
	p += (mlength - pos);
	r += (mlength - pos);
	if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
	Uintlist_head_set(npositions,pos);
      } else {
	pos -= mlength;
      }

      mlength = Uintlist_head(npositions);
      if (mlength > max_softclip) {
	/* Make sure initial soft clip does not extend past max_softclip */
	fprintf(stderr,"Read has initial soft clip %d greater than max_softclip %d\n",mlength,max_softclip);
	pos += (mlength - max_softclip);
	p += (mlength - max_softclip);
	r += (mlength - max_softclip);
	if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
	Uintlist_head_set(npositions,max_softclip);
      }
    }
  }

  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S' && max_softclip == 0) {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;
      r += mlength;
      if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
      /* r += mlength; -- hard clip means quality string is absent */
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'N') {
      pos += mlength;

    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'I') {

      if (print_indels_p == true) {
	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing insertion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = alloc_tallies[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	/* quality_score = (int) count_average(counts,r,mlength); */
	revise_insertion(pos,/*query_insert*/p,mlength,signed_shift,nm,xs,this,nreps);
      }

      p += mlength;
      r += mlength;
      if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
      shift += ((strand == '+') ? +mlength : -mlength);

    } else if (type == 'D') {

      if (print_indels_p == true) {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  /* q = &(*p); */
	  fprintf(stderr,"Need genome to print deletions\n");
	  exit(9);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing deletion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = alloc_tallies[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	revise_deletion(pos,/*deletion*/q,mlength,signed_shift,nm,xs,this,nreps);
	/* Revise deletion counts for coverage */
	if (signed_shift > 0) {
	  for (i = 0; i < mlength; i++) {
	    this = alloc_tallies[alloci+i];
	    this->delcounts_plus += nreps;
	  }
	} else {
	  for (i = 0; i < mlength; i++) {
	    this = alloc_tallies[alloci+i];
	    this->delcounts_minus += nreps;
	  }
	}
      }

      pos += mlength;

    } else if (type == 'M' || (type == 'S' && max_softclip > 0)) {
      if (0 /* mlength < min_mlength */) {
	p += mlength;
	r += mlength;
	if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	if (type == 'S' && mlength > max_softclip) {
	  /* Must be final softclip, because we handled initial one already */
	  fprintf(stderr,"Read has final soft clip %d greater than max_softclip %d\n",mlength,max_softclip);
	  mlength = max_softclip;
	}

	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  q = &(*p);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	/* psave = p; qsave = q; */
	debug1(printf("Processing query %.*s and genomic %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));
	  this = alloc_tallies[alloci];
	  if (mlength == 0) {
	    right = (Tally_T) NULL;
	  } else {
	    right = alloc_tallies[alloci+1];
	  }

	  if (genome == NULL) {
	    this->refnt = ' ';

	  } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	    this->refnt = toupper(*q);
	    debug1(printf("Line %d assigning %c at %u to tally %p\n",linei,this->refnt,pos+1U,this));
	    
	  } else if (this->refnt != toupper(*q)) {
	    fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		    this->refnt,*q,pos+1U);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    this->nmatches,List_length(this->mismatches_byshift));
	    abort();
	  }

	  signed_shift = (strand == '+') ? shift : -shift;
	  if (counts) {
	    ncounts = block[r % 8];
	  } else {
	    ncounts = nreps;
	  }
	  debug1(printf("ncounts at r %d is %d\n",r,ncounts));
	  revise_position(/*querynt*/*p,/*genomicnt*/*q,nm,xs,signed_shift,
			  this,right,ignore_query_Ns_p,readlevel_p,linei,ncounts);
	  if (readlevel_p == true && pos >= chrstart && pos <= chrend) {
	    FREE(genomic);
	    return;
	  }

	  p++;
	  q++;
	  if (++r % 8 == 0) {
	    if (counts) _mm_store_si128((__m128i *) block,counts[r/8]);
	  }
	  pos++;
	  shift += ((strand == '+') ? +1 : -1);
	}
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }
  
  FREE(genomic);
  return;
}


static int
count_nsplices (Intlist_T types) {
  int nsplices = 0;
  Intlist_T p;

  for (p = types; p != NULL; p = Intlist_next(p)) {
    if (Intlist_head(p) == 'N') {
      nsplices++;
    }
  }
  return nsplices;
}


static void
revise_read_lh (Tally_T *alloc_tallies_low, Tally_T *alloc_tallies_high, Genomicpos_T chrstart, Genomicpos_T chrend,
		bool lowend_p, Genomicpos_T chrpos_low,
		unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
		char *shortread, char *quality_string, int nm, char splice_strand, bool terminalp,
		Genomicpos_T alloc_low, Genome_T genome, Genomicpos_T chroffset,
		bool ignore_query_Ns_p, bool print_indels_p, bool readlevel_p, int max_softclip,
		unsigned int linei) {
  Tally_T this, right;
  int alloci;
  int shift, signed_shift;
  char strand;
  Genomicpos_T pos;
  char *genomic = NULL, *p, *q;
  char *r;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;
  int quality_score;
  int xs;


  debug1(printf("Revising read at %u\n",chrpos_low));

  if (splice_strand == '+') {
    xs = 1;
  } else if (splice_strand == '-') {
    xs = 2;
  } else {
    xs = 0;
  }


  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
  }

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */

  if (terminalp == false) {
    /* Don't handle soft clip for terminal alignments */
    max_softclip = 0;
  }

  p = shortread;
  r = quality_string;
  if (max_softclip > 0 && types != NULL) {
    if (Intlist_head(types) == 'S') {
      /* Revise pos, so we handle the initial soft clip */
      mlength = Uintlist_head(npositions);
      if (pos < mlength) {
	/* Make sure initial soft clip does not extend past beginning of chromosome */
	pos = 0U;
	p += (mlength - pos);
	r += (mlength - pos);
	Uintlist_head_set(npositions,pos);
      } else {
	pos -= mlength;
      }

      mlength = Uintlist_head(npositions);
      if (mlength > max_softclip) {
	/* Make sure initial soft clip does not extend past max_softclip */
	fprintf(stderr,"Read has initial soft clip %d greater than max_softclip %d\n",mlength,max_softclip);
	pos += (mlength - max_softclip);
	p += (mlength - max_softclip);
	r += (mlength - max_softclip);
	Uintlist_head_set(npositions,max_softclip);
      }
    }
  }

  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S' && max_softclip == 0) {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
      /* r += mlength; -- hard clip means quality string is absent */
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'N') {
      pos += mlength;

    } else if (type == 'P') {
      /* Do nothing */

    } else if (type == 'I') {

      if (print_indels_p == true) {
	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing insertion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	if (quality_string == NULL) {
	  quality_score = DEFAULT_QUALITY;
	} else {
	  quality_score = (int) average(r,mlength);
	}

	revise_insertion(pos,/*query_insert*/p,mlength,signed_shift,nm,xs,this,/*nreps*/1);
      }

      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);

    } else if (type == 'D') {

      if (print_indels_p == true) {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  /* q = &(*p); */
	  fprintf(stderr,"Need genome to print deletions\n");
	  exit(9);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing deletion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	revise_deletion(pos,/*deletion*/q,mlength,signed_shift,nm,xs,this,/*nreps*/1);
      }

      pos += mlength;

    } else if (type == 'M' || (type == 'S' && max_softclip > 0)) {
      if (0 /* mlength < min_mlength */) {
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	if (type == 'S' && mlength > max_softclip) {
	  /* Must be final softclip, because we handled initial one already */
	  fprintf(stderr,"Read has final soft clip %d greater than max_softclip %d\n",mlength,max_softclip);
	  mlength = max_softclip;
	}

	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  q = &(*p);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	/* psave = p; qsave = q; */
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));

	  signed_shift = (strand == '+') ? shift : -shift;
	  if (quality_string == NULL) {
	    quality_score = DEFAULT_QUALITY;
	  } else {
	    quality_score = (int) *r;
	  }

	  this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];
	  if (mlength == 0) {
	    right = (Tally_T) NULL;
	  } else {
	    right = (lowend_p == true) ? alloc_tallies_low[alloci+1] : alloc_tallies_high[alloci+1];
	  }

	  if (genome == NULL) {
	    this->refnt = ' ';

	  } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	    this->refnt = toupper(*q);
	    
	  } else if (this->refnt != toupper(*q)) {
	    fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		    this->refnt,*q,pos+1U);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    this->nmatches,List_length(this->mismatches_byshift));
	    abort();
	  }

	  revise_position(/*querynt*/*p,/*genomicnt*/*q,nm,xs,signed_shift,
			  this,right,ignore_query_Ns_p,readlevel_p,linei,/*ncounts*/1);
	  if (readlevel_p == true && pos >= chrstart && pos <= chrend) {
	    FREE(genomic);
	    return;
	  }

	  p++;
	  q++;
	  r++;
	  pos++;
	  shift += ((strand == '+') ? +1 : -1);
	}
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }
  
  FREE(genomic);
  return;
}


/* Also appears in gstruct.c */
static bool
best_mapping_p (Tableuint_T resolve_low_table, Tableuint_T resolve_high_table, char *acc,
		Genomicpos_T chrpos_low, Genomicpos_T chroffset) {
  Genomicpos_T genomicpos, genomicpos_low, genomicpos_high;

  genomicpos_low = Tableuint_get(resolve_low_table,(void *) acc);
  genomicpos_high = Tableuint_get(resolve_high_table,(void *) acc);


  if (genomicpos_low == 0 && genomicpos_high == 0) {
    /* Was already unique.  Could also check nhits. */
    return true;
  } else {
    genomicpos = chroffset + chrpos_low;
    if (genomicpos == genomicpos_low || genomicpos == genomicpos_high) {
      return true;
    } else {
      return false;
    }
  }
}


static void
get_passing_counts (__m128i *counts, int ncounts, int readlength, Bamline_T *bamlines, int nreps, int minimum_quality_score) {
  int i, r, b;
  char *quality_string, buffer[16];
  int x;

  __m128i block, cmp16, cmp8, threshold, ones;
  __m128i foo;

  threshold = _mm_set1_epi8(minimum_quality_score);
  ones = _mm_set1_epi8(1);

  /* Initialize counts */
  b = 0; r = 0;
  while (r + 16 < readlength) {
    counts[b++] = _mm_set1_epi16(0);
    counts[b++] = _mm_set1_epi16(0);
    r += 16;
  }
  if (r < readlength) {
    counts[b++] = _mm_set1_epi16(0);
    r += 8;
    if (r < readlength) {
      counts[b] = _mm_set1_epi16(0);
    }
  }
  if (b > ncounts) {
    fprintf(stderr,"Allocated only %d vectors, but need %d\n",ncounts,b);
    abort();
  }

  for (i = 0; i < nreps; i++) {
    quality_string = Bamline_quality_string(bamlines[0]);
    /* printf("quality string: %s\n",quality_string); */
    b = 0; r = 0;
    while (r + 16 < readlength) {
      block = _mm_loadu_si128((__m128i *) &(quality_string[r])); /* Put 16 bytes of quality string into block */

      /* Count bytes where (threshold > quality).  "true": 0xFF => 0.
	 "false": 0x00 => 1.  Equivalently, cmp8 contains a 1 whenever
	 quality >= threshold. */
      cmp8 = _mm_add_epi8(_mm_cmpgt_epi8(threshold,block),ones);

      /* Tally the lower 8 bytes */
      cmp16 = _mm_cvtepi8_epi16(cmp8);
      counts[b] = _mm_add_epi16(counts[b],cmp16); b++;

      /* Tally the upper 8 bytes */
      cmp16 = _mm_cvtepi8_epi16(_mm_srli_si128(cmp8,8));
      counts[b] = _mm_add_epi16(counts[b],cmp16); b++;

      r += 16;
    }
    if (r < readlength) {
      /* Handle the last 8 bytes */
      strncpy(buffer,&(quality_string[r]),readlength-r);
      block = _mm_loadu_si128((__m128i *) buffer);
      cmp8 = _mm_add_epi8(_mm_cmpgt_epi8(threshold,block),ones);
      cmp16 = _mm_cvtepi8_epi16(cmp8);
      counts[b] = _mm_add_epi16(counts[b],cmp16); b++;
      r += 8;

      if (r < readlength) {
	cmp16 = _mm_cvtepi8_epi16(_mm_srli_si128(cmp8,8));
	counts[b] = _mm_add_epi16(counts[b],cmp16); /* b++; */
      }
    }

#if 0
    b = 0; r = 0;
    while (r + 16 < readlength) {
      for (x = 0; x < 16; x++) {
	printf("%d ",quality_string[r+x]);
      }
      printf("\n");

      foo = counts[b++];
      print_vector_16_dec(foo);
      foo = counts[b++];
      print_vector_16_dec(foo);
      printf("\n");

      r += 16;
    }

    if (r < readlength) {
      x = 0;
      while (r + x < readlength) {
	printf("%d ",quality_string[r+x]);
	x++;
      }
      printf("\n");

      foo = counts[b++];
      print_vector_16_dec(foo);
      r += 8;
      if (r < readlength) {
	foo = counts[b++];
	print_vector_16_dec(foo);
      }
      printf("\n");
    }

    exit(0);
#endif
  }

  return;

}



/* A modified version, Indelfix_run, is in indelfix.c */
/* alloc_ptr is the maximal position where information has been seen so far */
long int
Bamtally_run (long int **tally_matches, long int **tally_mismatches,
	      List_T *intervallist, List_T *labellist, List_T *datalist,
	      Bamreader_T bamreader, Genome_T genome, char *printchr,
	      Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, IIT_T map_iit,
	      int alloclength, Tableuint_T resolve_low_table, Tableuint_T resolve_high_table,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq,
	      int minimum_quality_score, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      bool ignore_lowend_p, bool ignore_highend_p,
	      Tally_outputtype_T output_type, bool blockp, int blocksize,
	      int quality_score_adj, int min_depth, int variant_strands,
	      bool genomic_diff_p, bool signed_counts_p, bool ignore_query_Ns_p,
	      bool print_indels_p, bool print_totals_p, bool print_cycles_p, bool print_nm_scores_p, bool print_xs_scores_p,
	      bool want_genotypes_p, bool verbosep, bool readlevel_p, int max_softclip, bool print_noncovered_p,
	      char *bamfile) {
  long int grand_total = 0;
  Tally_T this;
  Genomicpos_T alloc_ptr, alloc_low, alloc_high, chrpos_low, chrpos_high, chrpos;
  Genomicpos_T blockptr = 0U, blockstart = 0U, blockend = 0U;
  Genomicpos_T prev_blockptr = 0U, prev_blockstart = 0U, prev_blockend = 0U;
  int delta, alloci;
  bool goodp;
  Bamline_T *bamlines, bamline, bamline_rep, prev_bamline = NULL;
  char *rsequence;
  int readlength;
  unsigned int linei = 0, linei_start, linei_end;
  int nlines;
  int nreps;
  __m128i *counts, foo;
  int ncounts;
  int i;


  Genomicpos_T exonstart = 0U, lastpos = chrstart;

  Tally_T *alloc_tallies, *block_tallies, *alloc_tallies_alloc, *block_tallies_alloc;
  Tally_T *prev_block_tallies, *prev_block_tallies_alloc;
  
#if 0
  /* How gmapR calls Bamtally_iit */
  minimum_mapq = 13;
  variant_strands = 1;
  min_depth = 0;
  print_indels_p = true;
  ignore_duplicates_p = true;
#endif

  /*
  fprintf(stderr,"Called with alloclength %d, minimum_mapq %d, good_unique_mapq %d, minimum_quality_score %d, maximum_nhits %d, need_concordant_p %d, need_unique_p %d, need_primary_p %d, ignore_duplicates_p %d, min_depth %d, variant_strands %d, ignore_query_Ns_p %d, blocksize %d, verbosep %d, readlevel_p %d, max_softclip %d\n",
	  alloclength,minimum_mapq,good_unique_mapq,minimum_quality_score,maximum_nhits,need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,min_depth,variant_strands,ignore_query_Ns_p,blocksize,verbosep,readlevel_p,max_softclip);
  */


  /* Create tally at position N to store n_fromleft */
  alloc_tallies_alloc = (Tally_T *) CALLOC(alloclength + 2*max_softclip + 1,sizeof(Tally_T));
  for (i = 0; i < alloclength + 2*max_softclip + 1; i++) {
    alloc_tallies_alloc[i] = Tally_new();
  }
  alloc_tallies = &(alloc_tallies_alloc[0]);

  block_tallies_alloc = (Tally_T *) CALLOC(blocksize+1,sizeof(Tally_T));
  for (i = 0; i < blocksize+1; i++) {
    block_tallies_alloc[i] = Tally_new();
  }
  block_tallies = &(block_tallies_alloc[0]);

  prev_block_tallies_alloc = (Tally_T *) CALLOC(blocksize+1,sizeof(Tally_T));
  for (i = 0; i < blocksize+1; i++) {
    prev_block_tallies_alloc[i] = Tally_new();
  }
  prev_block_tallies = &(prev_block_tallies_alloc[0]);

  *tally_matches = (long int *) NULL;
  *tally_mismatches = (long int *) NULL;
  *intervallist = (List_T) NULL;
  *labellist = (List_T) NULL;
  *datalist = (List_T) NULL;

  if (output_type == OUTPUT_TALLY) {
    *tally_matches = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));
    *tally_mismatches = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));
  }
  
  alloc_ptr = 0U;
  alloc_low = 0U;
  alloc_high = alloc_low + alloclength + 2*max_softclip;
  goodp = false;

  blockstart = 0U;
  blockend = 0U;
  blockptr = 0U;

  while (goodp == false &&
	 (bamlines = Bamread_next_bamline_set(&nlines,&prev_bamline,bamreader,
					      desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					      need_unique_p,need_primary_p,ignore_duplicates_p,
					      need_concordant_p)) != NULL) {
    for (linei_start = 0; linei_start < nlines; linei_start = linei_end) {
      linei_end = linei_start + 1;
      while (linei_end < nlines && !strcmp(Bamline_read(bamlines[linei_start]),Bamline_read(bamlines[linei_end]))) {
	linei_end++;
      }

      bamline_rep = bamlines[linei_start];
      nreps = linei_end - linei_start;
      
      debug0(printf("\n"));
      debug0(printf(">%s:%u..%u ",printchr,Bamline_chrpos_low(bamline_rep),Bamline_chrpos_high(bamline_rep)));
      debug0(Bamread_print_cigar(stdout,bamline_rep));
      debug0(printf("\n"));

      if (resolve_low_table != NULL &&
	  best_mapping_p(resolve_low_table,resolve_high_table,
			 Bamline_acc(bamline_rep),Bamline_chrpos_low(bamline_rep),
			 chroffset) == false) {
	/* Skip */
      } else {
	chrpos_low = Bamline_chrpos_low(bamline_rep);
	chrpos_high = Bamline_chrpos_high(bamline_rep);

	alloc_low = (chrpos_low < max_softclip) ? 0U : (chrpos_low - max_softclip);
	alloc_high = alloc_low + alloclength + 2*max_softclip;
	alloc_ptr = alloc_low;

	debug0(printf("    initialize alloc_low %u, alloc_high = %u, blockstart = %u, blockend = %u\n",
		      alloc_low,alloc_high,blockstart,blockend));

	if (chrpos_high + max_softclip >= alloc_high) {
	  if (verbosep == true) {
	    fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		    Bamline_acc(bamline_rep),printchr,chrpos_low,chrpos_high,alloclength);
	  }
	} else {
	  if (chrpos_high + max_softclip + 1U > alloc_ptr) {
	    alloc_ptr = chrpos_high + max_softclip + 1U;
	    debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	  }

	  rsequence = Bamline_read(bamline_rep);
	  readlength = Bamline_readlength(bamline_rep);

	  if (Bamline_quality_string(bamline_rep) == NULL) {
	    revise_read(alloc_tallies,chrstart,chrend,chrpos_low,Bamline_flag(bamline_rep),
			Bamline_cigar_types(bamline_rep),Bamline_cigar_npositions(bamline_rep),
			Bamline_cigar_querylength(bamline_rep),rsequence,Bamline_nm(bamline_rep),
			Bamline_splice_strand(bamline_rep),Bamline_terminalp(bamline_rep),alloc_low,
			genome,chroffset,ignore_query_Ns_p,print_indels_p,readlevel_p,
			max_softclip,linei+linei_start,/*counts*/NULL,nreps);
	  } else {
	    ncounts = (readlength + 7)/8;
	    counts = (__m128i *) _mm_malloc(ncounts * sizeof(__m128i),16);
	    get_passing_counts(counts,ncounts,readlength,&(bamlines[linei_start]),nreps,minimum_quality_score);
	    revise_read(alloc_tallies,chrstart,chrend,chrpos_low,Bamline_flag(bamline_rep),
			Bamline_cigar_types(bamline_rep),Bamline_cigar_npositions(bamline_rep),
			Bamline_cigar_querylength(bamline_rep),rsequence,Bamline_nm(bamline_rep),
			Bamline_splice_strand(bamline_rep),Bamline_terminalp(bamline_rep),alloc_low,
			genome,chroffset,ignore_query_Ns_p,print_indels_p,readlevel_p,
			max_softclip,linei+linei_start,counts,nreps);
	    _mm_free(counts);
	  }
	  goodp = true;
	}
      }
    }

    for (linei_start = 0; linei_start < nlines; linei_start++) {
      bamline = bamlines[linei_start];
      Bamline_free(&bamline);
    }
    FREE(bamlines);

    linei += nlines;
  }


  while ((bamlines = Bamread_next_bamline_set(&nlines,&prev_bamline,bamreader,
					      desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					      need_unique_p,need_primary_p,ignore_duplicates_p,
					      need_concordant_p)) != NULL) {
    for (linei_start = 0; linei_start < nlines; linei_start = linei_end) {
      linei_end = linei_start + 1;
      while (linei_end < nlines && !strcmp(Bamline_read(bamlines[linei_start]),Bamline_read(bamlines[linei_end]))) {
	linei_end++;
      }

      bamline_rep = bamlines[linei_start];
      nreps = linei_end - linei_start;

      debug0(printf("*  alloc: %u..%u..%u  block: %u..%u..%u  lastpos: %u ",
		    alloc_low,alloc_ptr,alloc_high,blockstart,blockptr,blockend,lastpos));
      debug0(printf("\n"));
      debug0(printf(">%s:%u..%u ",printchr,Bamline_chrpos_low(bamline_rep),Bamline_chrpos_high(bamline_rep)));
      debug0(Bamread_print_cigar(stdout,bamline_rep));
      debug0(printf("\n"));
    
      if (resolve_low_table != NULL &&
	  best_mapping_p(resolve_low_table,resolve_high_table,
			 Bamline_acc(bamline_rep),Bamline_chrpos_low(bamline_rep),
			 chroffset) == false) {
	/* Skip */
      } else if (ignore_lowend_p == true && Bamline_lowend_p(bamline_rep) == true) {
	/* Skip */
      } else if (ignore_highend_p == true && Bamline_lowend_p(bamline_rep) == false) {
	/* Skip */
      } else {
	chrpos_low = Bamline_chrpos_low(bamline_rep);
	chrpos_high = Bamline_chrpos_high(bamline_rep);

	if (chrpos_low > alloc_ptr + max_softclip) {
	  debug0(printf("Case 1: chrpos_low %u > alloc_ptr %u + max_softclip %d\n",chrpos_low,alloc_ptr,max_softclip));
	  debug0(printf("    transfer from alloc_low %u to alloc_ptr %u\n",alloc_low,alloc_ptr));
	  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
	    alloci = chrpos - alloc_low;
	    this = alloc_tallies[alloci];
	    if (this->nmatches > 0 || this->delcounts_plus + this->delcounts_minus > 0 ||
		this->mismatches_byshift != NULL || this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
	      lastpos = transfer_position(&grand_total,alloc_tallies,block_tallies,&exonstart,lastpos,
					  &blockptr,&blockstart,&blockend,
					  prev_block_tallies,&prev_blockptr,&prev_blockstart,&prev_blockend,
					  *tally_matches,*tally_mismatches,
					  &(*intervallist),&(*labellist),&(*datalist),
					  chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
					  map_iit,output_type,blockp,blocksize,
					  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
					  signed_counts_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
					  want_genotypes_p,readlevel_p,print_noncovered_p);
	    }
	  }

	  debug0(printf("    reset alloc_low = chrpos_low - max_softclip\n"));
	  alloc_low = (chrpos_low < max_softclip) ? 0U : chrpos_low - max_softclip;
	  alloc_high = alloc_low + alloclength + 2*max_softclip;
	  alloc_tallies[alloclength + 2*max_softclip]->n_fromleft_plus = 0;
	  alloc_tallies[alloclength + 2*max_softclip]->n_fromleft_minus = 0;
	  
	} else if (chrpos_low > alloc_low + max_softclip) {
	  debug0(printf("Case 2: chrpos_low %u > alloc_low %u + max_softclip %d\n",chrpos_low,alloc_low,max_softclip));
	  debug0(printf("    transfer from alloc_low %u up to chrpos_low %u\n",alloc_low,chrpos_low));
	  for (chrpos = alloc_low; chrpos + max_softclip < chrpos_low; chrpos++) {
	    alloci = chrpos - alloc_low;
	    this = alloc_tallies[alloci];
	    if (this->nmatches > 0 || this->delcounts_plus + this->delcounts_minus > 0 ||
		this->mismatches_byshift != NULL || this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
	      lastpos = transfer_position(&grand_total,alloc_tallies,block_tallies,&exonstart,lastpos,
					  &blockptr,&blockstart,&blockend,
					  prev_block_tallies,&prev_blockptr,&prev_blockstart,&prev_blockend,
					  *tally_matches,*tally_mismatches,
					  &(*intervallist),&(*labellist),&(*datalist),
					  chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
					  map_iit,output_type,blockp,blocksize,
					  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
					  signed_counts_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
					  want_genotypes_p,readlevel_p,print_noncovered_p);
	    }
	  }

	  if (chrpos_low < max_softclip) {
	    delta = chrpos_low /* - alloc_low (should be 0) */;
	  } else {
	    delta = chrpos_low - (alloc_low + max_softclip);
	  }
	  debug0(printf("    shift alloc upward by %d so alloc_low = chrpos_low - max_softclip %d\n",delta,max_softclip));
	  for (alloci = 0; alloci < (int) (alloc_ptr - alloc_low - delta); alloci++) {
	    Tally_transfer(&(alloc_tallies[alloci]),&(alloc_tallies[alloci+delta]));
	  }
	  alloc_tallies[alloc_ptr - alloc_low - delta]->n_fromleft_plus = alloc_tallies[alloc_ptr - alloc_low]->n_fromleft_plus;
	  alloc_tallies[alloc_ptr - alloc_low - delta]->n_fromleft_minus = alloc_tallies[alloc_ptr - alloc_low]->n_fromleft_minus;
	  for ( ; alloci < (int) (alloc_ptr - alloc_low); alloci++) {
	    debug0(printf("Resetting alloci %d\n",alloci));
	    Tally_clear(alloc_tallies[alloci]);
	  }
	  alloc_low = (chrpos_low < max_softclip) ? 0U : (chrpos_low - max_softclip);
	  alloc_high = alloc_low + alloclength + 2*max_softclip;
	
	} else if (chrpos_low < max_softclip) {
	  debug0(printf("Case 3a: chrpos_low %u < max_softclip %d\n",chrpos_low,max_softclip));
	  
	} else if (chrpos_low == alloc_low + max_softclip) {
	  debug0(printf("Case 3b: chrpos_low %u == alloc_low %u + max_softclip %d\n",chrpos_low,alloc_low,max_softclip));

	} else {
	  fprintf(stderr,"Sequences are not in order.  Got chrpos_low %u < alloc_low %u + max_softclip %d\n",
		  chrpos_low,alloc_low,max_softclip);
	  abort();
	}

	if (chrpos_high + max_softclip >= alloc_high) {
	  if (verbosep == true) {
	    fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		    Bamline_acc(bamline_rep),printchr,chrpos_low,chrpos_high,alloc_high);
	  }
	} else {
	  if (chrpos_high + max_softclip + 1U > alloc_ptr) {
	    alloc_ptr = chrpos_high + max_softclip + 1U;
	    debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	  }
	  
	  if (Bamline_quality_string(bamline_rep) == NULL) {
	    revise_read(alloc_tallies,chrstart,chrend,chrpos_low,Bamline_flag(bamline_rep),
			Bamline_cigar_types(bamline_rep),Bamline_cigar_npositions(bamline_rep),
			Bamline_cigar_querylength(bamline_rep),Bamline_read(bamline_rep),
			Bamline_nm(bamline_rep),Bamline_splice_strand(bamline_rep),Bamline_terminalp(bamline_rep),alloc_low,
			genome,chroffset,ignore_query_Ns_p,print_indels_p,readlevel_p,
			max_softclip,linei + linei_start,/*counts*/NULL,nreps);
	  } else {
	    readlength = Bamline_readlength(bamline_rep);
	    ncounts = (readlength + 7)/8;
	    counts = (__m128i *) _mm_malloc(ncounts * sizeof(__m128i),16);
	    get_passing_counts(counts,ncounts,readlength,&(bamlines[linei_start]),nreps,minimum_quality_score);
	    revise_read(alloc_tallies,chrstart,chrend,chrpos_low,Bamline_flag(bamline_rep),
			Bamline_cigar_types(bamline_rep),Bamline_cigar_npositions(bamline_rep),
			Bamline_cigar_querylength(bamline_rep),Bamline_read(bamline_rep),
			Bamline_nm(bamline_rep),Bamline_splice_strand(bamline_rep),Bamline_terminalp(bamline_rep),alloc_low,
			genome,chroffset,ignore_query_Ns_p,print_indels_p,readlevel_p,
			max_softclip,linei + linei_start,counts,nreps);
	  }
	}
      }
    }

    for (linei_start = 0; linei_start < nlines; linei_start++) {
      bamline = bamlines[linei_start];
      Bamline_free(&bamline);
    }
    FREE(bamlines);

    linei += nlines;
  }

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));
  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
    alloci = chrpos - alloc_low;
    this = alloc_tallies[alloci];
    if (this->nmatches > 0 || this->delcounts_plus + this->delcounts_minus > 0 ||
	this->mismatches_byshift != NULL || this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
      lastpos = transfer_position(&grand_total,alloc_tallies,block_tallies,&exonstart,lastpos,
				  &blockptr,&blockstart,&blockend,
				  prev_block_tallies,&prev_blockptr,&prev_blockstart,&prev_blockend,
				  *tally_matches,*tally_mismatches,
				  &(*intervallist),&(*labellist),&(*datalist),
				  chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
				  map_iit,output_type,blockp,blocksize,
				  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
				  signed_counts_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				  want_genotypes_p,readlevel_p,print_noncovered_p);
    }
  }

  debug0(printf("print block from blockstart %u to blockptr %u\n",blockstart,blockptr));
  if (output_type == OUTPUT_RUNLENGTHS) {
    lastpos = print_runlength(block_tallies,&exonstart,lastpos,blockstart,blockptr,printchr);
  } else if (output_type == OUTPUT_TALLY) {
    tally_block(*tally_matches,*tally_mismatches,
		block_tallies,blockstart,blockptr,genome,printchr,chroffset,chrstart,
		quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		print_noncovered_p);
  } else if (output_type == OUTPUT_IIT) {
    iit_block(&(*intervallist),&(*labellist),&(*datalist),
	      block_tallies,blockstart,blockptr,
	      prev_block_tallies,prev_blockstart,prev_blockptr,
	      genome,printchr,chroffset,map_iit,
	      quality_score_adj,min_depth,variant_strands,
	      print_cycles_p,print_nm_scores_p,print_xs_scores_p,print_noncovered_p);
  } else if (output_type == OUTPUT_TOTAL) {
    grand_total += block_total(block_tallies,blockstart,blockptr);
  } else{
    if (print_noncovered_p == true) {
      debug0(printf("Printing zeroes from lastpos %u to blockstart %u\n",lastpos,blockstart));
      print_zeroes(lastpos,blockstart,printchr,blocksize,genome,chroffset,blockp);
    }
    print_block(block_tallies,blockstart,blockptr,
		prev_block_tallies,prev_blockstart,prev_blockptr,
		genome,printchr,chroffset,map_iit,
		blockp,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		signed_counts_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
		want_genotypes_p,readlevel_p,print_noncovered_p);
    if (print_noncovered_p == true) {
      debug0(printf("Printing zeroes from blockptr %u to chrend %u\n",blockptr,chrend));
      print_zeroes(blockptr,chrend,printchr,blocksize,genome,chroffset,blockp);
    }
  }


  for (i = 0; i < alloclength + 2*max_softclip + 1; i++) {
    Tally_free(&(alloc_tallies_alloc[i]));
  }
  FREE(alloc_tallies_alloc);

  for (i = 0; i < blocksize+1; i++) {
    Tally_free(&(block_tallies_alloc[i]));
  }
  FREE(block_tallies_alloc);

  for (i = 0; i < blocksize+1; i++) {
    Tally_free(&(prev_block_tallies_alloc[i]));
  }
  FREE(prev_block_tallies_alloc);

  return grand_total;
}


/* Assumes output_type is OUTPUT_TALLY */
void
Bamtally_run_lh (long int **tally_matches_low, long int **tally_mismatches_low,
		 long int **tally_matches_high, long int **tally_mismatches_high,
		 Bamreader_T bamreader, Genome_T genome, char *printchr,
		 Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, 
		 IIT_T map_iit, int alloclength,
		 char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		 bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		 int blocksize, int quality_score_adj, int min_depth, int variant_strands,
		 bool genomic_diff_p, bool ignore_query_Ns_p, bool verbosep, bool readlevel_p,
		 int max_softclip, bool print_noncovered_p) {
  Tally_T this_low, this_high;
  Genomicpos_T alloc_ptr, alloc_low, alloc_high, chrpos_low, chrpos_high, chrpos;
  Genomicpos_T blockptr = 0U, blockstart = 0U, blockend = 0U;
  int delta, alloci;
  bool goodp;
  Bamline_T bamline;
  unsigned int linei = 0;
  int i;

  Tally_T *alloc_tallies_low, *alloc_tallies_high, *block_tallies_low, *block_tallies_high;
  Tally_T *alloc_tallies_low_alloc, *alloc_tallies_high_alloc, *block_tallies_low_alloc, *block_tallies_high_alloc;


  alloc_tallies_low_alloc = (Tally_T *) CALLOC(alloclength + 2*max_softclip + 1,sizeof(Tally_T));
  for (i = 0; i < alloclength + 2*max_softclip + 1; i++) {
    alloc_tallies_low_alloc[i] = Tally_new();
  }
  alloc_tallies_low = &(alloc_tallies_low_alloc[0]);

  alloc_tallies_high_alloc = (Tally_T *) CALLOC(alloclength+1,sizeof(Tally_T));
  for (i = 0; i < alloclength+1; i++) {
    alloc_tallies_high_alloc[i] = Tally_new();
  }
  alloc_tallies_high = &(alloc_tallies_high_alloc[0]);

  block_tallies_low_alloc = (Tally_T *) CALLOC(blocksize+1,sizeof(Tally_T));
  for (i = 0; i < blocksize+1; i++) {
    block_tallies_low_alloc[i] = Tally_new();
  }
  block_tallies_low = &(block_tallies_low_alloc[0]);

  block_tallies_high_alloc = (Tally_T *) CALLOC(blocksize+1,sizeof(Tally_T));
  for (i = 0; i < blocksize+1; i++) {
    block_tallies_high_alloc[i] = Tally_new();
  }
  block_tallies_high = &(block_tallies_high_alloc[0]);

  *tally_matches_low = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));
  *tally_mismatches_low = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));
  *tally_matches_high = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));
  *tally_mismatches_high = (long int *) CALLOC(chrend-chrstart+1,sizeof(long int));

  alloc_ptr = 0U;
  alloc_low = 0U;
  alloc_high = alloc_low + alloclength + 2*max_softclip;
  goodp = false;

  blockstart = 0U;
  blockend = 0U;
  blockptr = 0U;

  while (goodp == false &&
	 (bamline = Bamread_next_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,ignore_duplicates_p,
					 need_concordant_p)) != NULL) {
    debug0(printf("\n"));
    debug0(printf(">%s:%u..%u ",printchr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline)));
    debug0(Bamread_print_cigar(stdout,bamline));
    debug0(printf("\n"));

    chrpos_low = Bamline_chrpos_low(bamline);
    chrpos_high = Bamline_chrpos_high(bamline);
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */

    alloc_low = (chrpos_low < max_softclip) ? 0U : (chrpos_low - max_softclip);
    alloc_high = alloc_low + alloclength + 2*max_softclip;
    alloc_ptr = alloc_low;

    debug0(printf("    initialize alloc_low %u, alloc_high = %u, blockstart = %u, blockend = %u\n",
		  alloc_low,alloc_high,blockstart,blockend));

    if (chrpos_high + max_softclip >= alloc_high) {
      if (verbosep == true) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		Bamline_acc(bamline),printchr,chrpos_low,chrpos_high,alloclength);
      }
    } else {
      if (chrpos_high + max_softclip + 1U > alloc_ptr) {
	alloc_ptr = chrpos_high + max_softclip + 1U;
	debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
      }

      revise_read_lh(alloc_tallies_low,alloc_tallies_high,chrstart,chrend,
		     Bamline_lowend_p(bamline),chrpos_low,Bamline_flag(bamline),
		     Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		     Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		     Bamline_quality_string(bamline),Bamline_nm(bamline),Bamline_splice_strand(bamline),
		     Bamline_terminalp(bamline),alloc_low,
		     genome,chroffset,ignore_query_Ns_p,/*print_indels_p*/false,
		     readlevel_p,max_softclip,linei);

      goodp = true;
    }

    Bamline_free(&bamline);
    linei++;
  }

  while ((bamline = Bamread_next_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,ignore_duplicates_p,
					 need_concordant_p)) != NULL) {
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */

    debug0(printf("*  alloc: %u..%u..%u  block: %u..%u..%u  ",
		  alloc_low,alloc_ptr,alloc_high,blockstart,blockptr,blockend));
    debug0(printf("\n"));
    debug0(printf(">%s:%u..%u ",printchr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline)));
    debug0(Bamread_print_cigar(stdout,bamline));
    debug0(printf("\n"));
    
    chrpos_low = Bamline_chrpos_low(bamline);
    chrpos_high = Bamline_chrpos_high(bamline);

    if (chrpos_low > alloc_ptr + max_softclip) {
      debug0(printf("Case 1: chrpos_low %u > alloc_ptr %u + max_softclip %d\n",chrpos_low,alloc_ptr,max_softclip));
      debug0(printf("    transfer from alloc_low %u to alloc_ptr %u\n",alloc_low,alloc_ptr));
      for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
	alloci = chrpos - alloc_low;
	this_low = alloc_tallies_low[alloci];
	this_high = alloc_tallies_high[alloci];
	if (this_low->nmatches > 0 || this_low->delcounts_plus + this_low->delcounts_minus > 0 ||
	    this_low->mismatches_byshift != NULL || this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	    this_high->nmatches > 0 || this_high->delcounts_plus + this_high->delcounts_minus > 0 ||
	    this_high->mismatches_byshift != NULL || this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
	  transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			       &blockptr,&blockstart,&blockend,
			       *tally_matches_low,*tally_mismatches_low,
			       *tally_matches_high,*tally_mismatches_high,
			       chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
			       map_iit,blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
			       print_noncovered_p);
	}
      }

      debug0(printf("    reset alloc_low = chrpos_low - max_softclip\n"));
      alloc_low = (chrpos_low < max_softclip) ? 0U : (chrpos_low - max_softclip);
      alloc_high = alloc_low + alloclength + 2*max_softclip;
	
    } else if (chrpos_low > alloc_low + max_softclip) {
      debug0(printf("Case 2: chrpos_low %u > alloc_low %u + max_softclip %d\n",chrpos_low,alloc_low,max_softclip));
      debug0(printf("    transfer from alloc_low %u up to chrpos_low %u\n",alloc_low,chrpos_low));
      for (chrpos = alloc_low; chrpos + max_softclip < chrpos_low; chrpos++) {
	alloci = chrpos - alloc_low;
	this_low = alloc_tallies_low[alloci];
	this_high = alloc_tallies_high[alloci];
	if (this_low->nmatches > 0 || this_low->delcounts_plus + this_low->delcounts_minus > 0 ||
	    this_low->mismatches_byshift != NULL || this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	    this_high->nmatches > 0 || this_high->delcounts_plus + this_high->delcounts_minus > 0 ||
	    this_high->mismatches_byshift != NULL || this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
	  transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			       &blockptr,&blockstart,&blockend,
			       *tally_matches_low,*tally_mismatches_low,
			       *tally_matches_high,*tally_mismatches_high,
			       chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
			       map_iit,blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
			       print_noncovered_p);
	}
      }

      if (chrpos_low < max_softclip) {
	delta = chrpos_low /* - alloc_low (should be 0) */;
      } else {
	delta = chrpos_low - (alloc_low + max_softclip);
      }
      debug0(printf("    shift alloc upward by %d so alloc_low = chrpos_low - max_softclip\n",delta));
      for (alloci = 0; alloci < (int) (alloc_ptr - alloc_low - delta); alloci++) {
	Tally_transfer(&(alloc_tallies_low[alloci]),&(alloc_tallies_low[alloci+delta]));
	Tally_transfer(&(alloc_tallies_high[alloci]),&(alloc_tallies_high[alloci+delta]));
      }
      for ( ; alloci < (int) (alloc_ptr - alloc_low); alloci++) {
	debug0(printf("Resetting alloci %d\n",alloci));
	Tally_clear(alloc_tallies_low[alloci]);
	Tally_clear(alloc_tallies_high[alloci]);
      }
      alloc_low = (chrpos_low < max_softclip) ? 0U : (chrpos_low - max_softclip);
      alloc_high = alloc_low + alloclength + 2*max_softclip;
	
    } else if (chrpos_low < max_softclip) {
      debug0(printf("Case 3a: chrpos_low %u < max_softclip %d\n",chrpos_low,max_softclip));

    } else if (chrpos_low == alloc_low + max_softclip) {
      debug0(printf("Case 3b: chrpos_low %u == alloc_low %u + max_softclip %d\n",chrpos_low,alloc_low,max_softclip));

    } else {
      fprintf(stderr,"Sequences are not in order.  Got chrpos_low %u < alloc_low %u + max_softclip %d\n",
	      chrpos_low,alloc_low,max_softclip);
      abort();
    }

    if (chrpos_high + max_softclip >= alloc_high) {
      if (verbosep == true) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		Bamline_acc(bamline),printchr,chrpos_low,chrpos_high,alloc_high);
      }
    } else {
      if (chrpos_high + max_softclip + 1U > alloc_ptr) {
	alloc_ptr = chrpos_high + max_softclip + 1U;
	debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
      }

      revise_read_lh(alloc_tallies_low,alloc_tallies_high,chrstart,chrend,
		     Bamline_lowend_p(bamline),chrpos_low,Bamline_flag(bamline),
		     Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		     Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		     Bamline_quality_string(bamline),Bamline_nm(bamline),Bamline_splice_strand(bamline),
		     Bamline_terminalp(bamline),alloc_low,
		     genome,chroffset,ignore_query_Ns_p,/*print_indels_p*/false,
		     readlevel_p,max_softclip,linei);
    }

    Bamline_free(&bamline);
    linei++;
  }

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));
  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
    alloci = chrpos - alloc_low;
    this_low = alloc_tallies_low[alloci];
    this_high = alloc_tallies_high[alloci];
    if (this_low->nmatches > 0 || this_low->delcounts_plus + this_low->delcounts_minus > 0 ||
	this_low->mismatches_byshift != NULL || this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	this_high->nmatches > 0 || this_high->delcounts_plus + this_high->delcounts_minus > 0 ||
	this_high->mismatches_byshift != NULL || this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
      transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			   &blockptr,&blockstart,&blockend,
			   *tally_matches_low,*tally_mismatches_low,
			   *tally_matches_high,*tally_mismatches_high,
			   chrpos,alloci,genome,printchr,chroffset,chrstart,chrend,
			   map_iit,blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
			   print_noncovered_p);
    }
  }

  debug0(printf("print block from blockstart %u to blockptr %u\n",blockstart,blockptr));
  tally_block(*tally_matches_low,*tally_mismatches_low,
	      block_tallies_low,blockstart,blockptr,genome,printchr,chroffset,chrstart,
	      quality_score_adj,min_depth,variant_strands,genomic_diff_p,
	      print_noncovered_p);
  tally_block(*tally_matches_high,*tally_mismatches_high,
	      block_tallies_high,blockstart,blockptr,genome,printchr,chroffset,chrstart,
	      quality_score_adj,min_depth,variant_strands,genomic_diff_p,
	      print_noncovered_p);


  for (i = 0; i < alloclength + 2*max_softclip + 1; i++) {
    Tally_free(&(alloc_tallies_low_alloc[i]));
    Tally_free(&(alloc_tallies_high_alloc[i]));
  }
  FREE(alloc_tallies_low_alloc);
  FREE(alloc_tallies_high_alloc);

  for (i = 0; i < blocksize+1; i++) {
    Tally_free(&(block_tallies_low_alloc[i]));
    Tally_free(&(block_tallies_high_alloc[i]));
  }
  FREE(block_tallies_low_alloc);
  FREE(block_tallies_high_alloc);


  return;
}



IIT_T
Bamtally_iit (Bamreader_T bamreader, char *desired_chr, char *bam_lacks_chr,
	      Genomicpos_T chrstart, Genomicpos_T chrend,
	      Genome_T genome, IIT_T chromosome_iit, IIT_T map_iit, int alloclength,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq,
	      int minimum_quality_score, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      int min_depth, int variant_strands, bool ignore_query_Ns_p,
	      bool print_indels_p, int blocksize, bool verbosep, bool readlevel_p,
	      int max_softclip, bool print_cycles_p, bool print_nm_scores_p, bool print_xs_scores_p,
	      bool print_noncovered_p) {
  IIT_T iit;

  long int *tally_matches, *tally_mismatches;
  List_T divlist = NULL, typelist = NULL;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL, p;
  Interval_T interval;
  char *divstring, *typestring, *label;
  char *chr, *chrptr;
  int bam_lacks_chr_length;

  int index;
  Genomicpos_T chroffset, chrlength;
  Table_T intervaltable, labeltable, datatable;
  bool allocp;

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  datatable = Table_new(65522,Table_string_compare,Table_string_hash);

  if (bam_lacks_chr == NULL) {
    bam_lacks_chr_length = 0;
  } else {
    bam_lacks_chr_length = strlen(bam_lacks_chr);
  }

  if (desired_chr == NULL) {
    /* Entire genome */
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      fprintf(stderr,"  processing chromosome %s...",chr);
      divlist = List_push(divlist,(void *) chr);

      chroffset = Interval_low(IIT_interval(chromosome_iit,index));
      chrlength = Interval_length(IIT_interval(chromosome_iit,index));

      if (bam_lacks_chr == NULL) {
	chrptr = chr;
      } else if (strncmp(chr,bam_lacks_chr,bam_lacks_chr_length) == 0) {
	chrptr = &(chr[bam_lacks_chr_length]);
      } else {
	chrptr = chr;
      }

      Bamread_limit_region(bamreader,chrptr,/*chrstart*/1,/*chrend*/chrlength);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   bamreader,genome,/*printchr*/chr,chroffset,
		   /*chrstart*/1U,/*chrend*/chrlength,map_iit,alloclength,
		   /*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   desired_read_group,minimum_mapq,good_unique_mapq,minimum_quality_score,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_IIT,/*blockp*/false,blocksize,
		   /*quality_score_adj*/0,min_depth,variant_strands,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,ignore_query_Ns_p,
		   print_indels_p,/*print_totals_p*/false,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
		   /*want_genotypes_p*/false,verbosep,readlevel_p,max_softclip,
		   print_noncovered_p,/*bamfile*/NULL);
      /* Reverse lists so we can specify presortedp == true */
      Table_put(intervaltable,(void *) chr,List_reverse(intervallist));
      Table_put(labeltable,(void *) chr,List_reverse(labellist));
      Table_put(datatable,(void *) chr,List_reverse(datalist));
      Bamread_unlimit_region(bamreader);
      /* Don't free chr yet, since divlist contains it */
      fprintf(stderr,"done\n");
    }

  } else if ((index = IIT_find_linear(chromosome_iit,desired_chr)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in genome\n",desired_chr);
    Table_free(&datatable);
    Table_free(&labeltable);
    Table_free(&intervaltable);

    return NULL;

  } else {
    /* Single chromosomal region */
    
    divlist = List_push(divlist,(void *) desired_chr);
    if (bam_lacks_chr == NULL) {
      chrptr = desired_chr;
    } else if (strncmp(desired_chr,bam_lacks_chr,bam_lacks_chr_length) == 0) {
      chrptr = &(desired_chr[bam_lacks_chr_length]);
    } else {
      chrptr = desired_chr;
    }

    if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == true) {
      /* Returns true only if some intervals are found */
      index = IIT_find_linear(chromosome_iit,desired_chr);
      chroffset = Interval_low(IIT_interval(chromosome_iit,index));
      
      Bamtally_run(&tally_matches,&tally_mismatches,
                   &intervallist,&labellist,&datalist,
                   bamreader,genome,/*printchr*/chrptr,chroffset,chrstart,chrend,map_iit,
		   alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
                   desired_read_group,minimum_mapq,good_unique_mapq,minimum_quality_score,maximum_nhits,
                   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
                   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
                   /*output_type*/OUTPUT_IIT,/*blockp*/false,blocksize,
                   /*quality_score_adj*/0,min_depth,variant_strands,
                   /*genomic_diff_p*/false,/*signed_counts_p*/false,
                   ignore_query_Ns_p, print_indels_p,/*print_totals_p*/false,
                   print_cycles_p,print_nm_scores_p,/*print_xs_scores_p*/false,
                   /*want_genotypes_p*/false,verbosep,readlevel_p,max_softclip,
		   print_noncovered_p,/*bamfile*/NULL);
    }
    /* Reverse lists so we can specify presortedp == true */
    Table_put(intervaltable,(void *) desired_chr,List_reverse(intervallist));
    Table_put(labeltable,(void *) desired_chr,List_reverse(labellist));
    Table_put(datatable,(void *) desired_chr,List_reverse(datalist));
    Bamread_unlimit_region(bamreader);
  }


  divlist = List_reverse(divlist);

  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(divlist,divstring);

  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(NULL,typestring);

  iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,intervaltable,
		   labeltable,datatable,/*divsort*/NO_SORT,/*version*/IIT_LATEST_VERSION,
		   /*presortedp*/true);

  FREE(typestring);
  List_free(&typelist);
  FREE(divstring);
  /* May need to free chr inside of divlist */
  List_free(&divlist);

  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      datalist = Table_get(datatable,(void *) chr);
      /* Do not free data yet, since tally_iit points to it */
      List_free(&datalist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    datalist = Table_get(datatable,(void *) desired_chr);
    /* Do not free data yet, since tally_iit points to it */
    List_free(&datalist);
  }
  Table_free(&datatable);


  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      labellist = Table_get(labeltable,(void *) chr);
      for (p = labellist; p != NULL; p = List_next(p)) {
	label = (char *) List_head(p);
	FREE(label);
      }
      List_free(&labellist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    labellist = Table_get(labeltable,(void *) desired_chr);
    for (p = labellist; p != NULL; p = List_next(p)) {
      label = (char *) List_head(p);
      FREE(label);
    }
    List_free(&labellist);
  }
  Table_free(&labeltable);
  

  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      intervallist = Table_get(intervaltable,(void *) chr);
      for (p = intervallist; p != NULL; p = List_next(p)) {
	interval = (Interval_T) List_head(p);
	Interval_free(&interval);
      }
      List_free(&intervallist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    intervallist = Table_get(intervaltable,(void *) desired_chr);
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&intervallist);
  }
  Table_free(&intervaltable);

  return iit;
}

