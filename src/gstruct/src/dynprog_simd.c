static char rcsid[] = "$Id: dynprog_simd.c 133254 2014-04-15 18:56:28Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_simd.h"
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

#include "mem.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* F loop */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



/************************************************************************
 *   Debugging procedures
 ************************************************************************/

#define NEG_INFINITY_DISPLAY -99



#ifdef DEBUG15
/* For debugging of SIMD procedures*/
static void
print_vector_8 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score8_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15]);
  return;
}

static void
print_vector_16 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score16_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d\n",r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif


#if defined(DEBUG2) || defined(DEBUG14)
static void
Matrix8_print (Score8_T **matrix, int rlength, int glength, char *rsequence,
	       char *gsequence, char *gsequencealt,
#ifdef DEBUG14
	       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
	       bool revp, int lband, int uband) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

#if 0
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("  . ");
      } else if (j > i + uband) {
	printf("  . ");
      } else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Matrix16_print (Score16_T **matrix, int rlength, int glength, char *rsequence,
		char *gsequence, char *gsequencealt,
#ifdef DEBUG14
		int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		bool revp, int lband, int uband) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

#if 0
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("  . ");
      } else if (j > i + uband) {
	printf("  . ");
      } else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

#endif

#if defined(DEBUG2) || defined(DEBUG14)
static void
Directions8_print (Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		   int rlength, int glength, char *rsequence, char *gsequence, char *gsequence_alt,
#ifdef DEBUG14
		   int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		   bool revp, int lband, int uband) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

#if 0
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("      ");
    } else {
      if (revp == false) {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("     ");
      } else if (j > i + uband) {
	printf("     ");
      } else {
	if (directions_Egap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be HORIZ */
	  printf("H");
	}
	printf("|");
	if (directions_nogap[j][i] == DIAG) {
	  printf("D");
	} else if (directions_nogap[j][i] == HORIZ) {
	  printf("H");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
	printf("|");
	if (directions_Fgap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}

static void
Directions16_print (Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		    int rlength, int glength, char *rsequence, char *gsequence, char *gsequence_alt,
#ifdef DEBUG14
		    int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		    bool revp, int lband, int uband) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

#if 0
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("      ");
    } else {
      if (revp == false) {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("     ");
      } else if (j > i + uband) {
	printf("     ");
      } else {
	if (directions_Egap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be HORIZ */
	  printf("H");
	}
	printf("|");
	if (directions_nogap[j][i] == DIAG) {
	  printf("D");
	} else if (directions_nogap[j][i] == HORIZ) {
	  printf("H");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
	printf("|");
	if (directions_Fgap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}

#endif



#ifdef DEBUG14
static void
banded_matrix8_compare (Score8_T **matrix1, Score32_T **matrix2, int rlength, int glength,
			int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_8 + 30 && matrix2[c][r] <= NEG_INFINITY_8 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix8_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			      goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
	exit(9);
      }
    }
  }

  return;
}

static void
banded_matrix16_compare (Score16_T **matrix1, Score32_T **matrix2, int rlength, int glength,
			 int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			 bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_16 + 30 && matrix2[c][r] <= NEG_INFINITY_16 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix16_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
	exit(9);
      }
    }
  }

  return;
}
#endif

#ifdef DEBUG14
static void
banded_directions8_compare_nogap (Score8_T **matrix, Direction8_T **directions1, Direction32_T **directions2,
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}


static void
banded_directions16_compare_nogap (Direction16_T **directions1, Direction32_T **directions2,
				   int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}
#endif

#ifdef DEBUG14
static void
banded_directions8_compare_Egap (Score8_T **matrix1, Direction8_T **directions1, Direction32_T **directions2,
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Egap (Direction16_T **directions1, Direction32_T **directions2,
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}
#endif


#ifdef DEBUG14
static void
banded_directions8_compare_Fgap (Score8_T **matrix1, Direction8_T **directions1, Direction32_T **directions2,
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  exit(9);
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Fgap (Direction16_T **directions1, Direction32_T **directions2,
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  exit(9);
	}
      }
    }
  }

  return;
}
#endif


/************************************************************************
 *   End of debugging procedures
 ************************************************************************/

/************************************************************************
 *   Note: These procedures expect SIMD registers to start at
 *   coordinate 1, not 0, and are different from the SIMD procedures
 *   for GMAP/GSNAP.  This is because we want row 0 to be a special
 *   case with no gap penalties.
 ************************************************************************/


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score8_T **
aligned_score8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));
#endif

  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score8_T **
aligned_directions8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score8_T **
aligned_directions8_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));

  return matrix;
}



/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score16_T **
aligned_score16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));
#endif
  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score16_T **
aligned_directions16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score16_T **
aligned_directions16_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));

  return matrix;
}
#endif




#define T Dynprog_T


#if defined(HAVE_SSE2)
/* Modified from Dynprog_simd_8_upper.  Operates by columns. */
Score8_T **
Dynprog_simd_8 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
		Direction8_T ***directions_Fgap,
		T this, char *rsequence, char *gsequence, char *gsequence_alt,
		int rlength, int glength,
#ifdef DEBUG14
		int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		Mismatchtype_T mismatchtype, int open, int extend,
		int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF;	/* Need to have the ability to go past NEG_INFINITY */
  Score8_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
#ifndef HAVE_SSE4_1
  __m128i pairscores_best, all_128;
#endif
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, all_one_bits;
  __m128i dir_horiz;
  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score8_T *pairscores_col0;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG14
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#endif


  rlength_ceil = (int) ((rlength + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("Dynprog_simd_8: "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned.one.matrix_ptrs,this->aligned.one.matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_0,this->aligned.one.directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned.one.directions_ptrs_1,this->aligned.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions8_calloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_2,this->aligned.one.directions_space_2);

  /* Row 0 initialization */
  /* penalty = 0; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NCHARS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);

  lband_ceil = (int) ((lband + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;
  pairscores_col0 = (Score8_T *) _mm_malloc(lband_ceil * sizeof(Score8_T),16);


  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
#ifdef HAVE_SSE4_1
  /* pairscores_col0[0] = (Score8_T) 0; */
  /* Initializion just to lband causes errors in dir_horiz for Egap */
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r-1] = (Score8_T) NEG_INFINITY_8;
  }
#else
  /* pairscores_col0[0] = (Score8_T) 0+128; */
  /* Initializion just to lband causes errors in dir_horiz for Egap */
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r-1] = (Score8_T) NEG_INFINITY_8+128;
  }
#endif

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r-1] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }


  all_one_bits = _mm_set1_epi8(-1);

  FF = (int *) MALLOC((glength + 1) * sizeof(int));

  gap_open = _mm_set1_epi8((Score8_T) open);
  gap_extend = _mm_set1_epi8((Score8_T) extend);

  if (jump_late_p) {
    for (rlo = 1; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H .  To fill in first column of each
	 row block with diags, make E < H. */
      E_r_gap = _mm_set1_epi8(NEG_INFINITY_8);
      H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8+1);

      if ((c = rlo - lband) < 0) {
	c = 0;
      }
      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;
	} else {
	  na2 = revp ? gsequence[1-c] : gsequence[c-1];
	  na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
	  switch (na2) {
	  case 'A': pairscores_std_ptr = pairscores[0]; break;
	  case 'C': pairscores_std_ptr = pairscores[1]; break;
	  case 'G': pairscores_std_ptr = pairscores[2]; break;
	  case 'T': pairscores_std_ptr = pairscores[3]; break;
	  default: pairscores_std_ptr = pairscores[4];
	  }
	  switch (na2_alt) {
	  case 'A': pairscores_alt_ptr = pairscores[0]; break;
	  case 'C': pairscores_alt_ptr = pairscores[1]; break;
	  case 'G': pairscores_alt_ptr = pairscores[2]; break;
	  case 'T': pairscores_alt_ptr = pairscores[3]; break;
	  default: pairscores_alt_ptr = pairscores[4];
	  }
	}

	if (rlo == 1) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (c == 0) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi8(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);

#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) > c + lband) {
	  rhigh_calc = c + lband;
	}
	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 1) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = score_column[rlo_calc-1];
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for (r = rlo_calc; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }

  } else {
    /* jump early */
    for (rlo = 1; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with diags, make E <= H. */
      E_r_gap = _mm_set1_epi8(NEG_INFINITY_8);
      H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8+0);

      if ((c = rlo - lband) < 0) {
	c = 0;
      }
      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;
	} else {
	  na2 = revp ? gsequence[1-c] : gsequence[c-1];
	  na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
	  switch (na2) {
	  case 'A': pairscores_std_ptr = pairscores[0]; break;
	  case 'C': pairscores_std_ptr = pairscores[1]; break;
	  case 'G': pairscores_std_ptr = pairscores[2]; break;
	  case 'T': pairscores_std_ptr = pairscores[3]; break;
	  default: pairscores_std_ptr = pairscores[4];
	  }
	  switch (na2_alt) {
	  case 'A': pairscores_alt_ptr = pairscores[0]; break;
	  case 'C': pairscores_alt_ptr = pairscores[1]; break;
	  case 'G': pairscores_alt_ptr = pairscores[2]; break;
	  case 'T': pairscores_alt_ptr = pairscores[3]; break;
	  default: pairscores_alt_ptr = pairscores[4];
	  }
	}

	if (rlo == 1) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (c == 0) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi8(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	/* Compare H + open with vert */
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi8(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);

#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) > c + lband) {
	  rhigh_calc = c + lband;
	}
	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 1) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = score_column[rlo_calc-1];
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for (r = rlo_calc; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }
  }

#ifdef DEBUG2
  printf("SIMD\n");
  Matrix8_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
#ifdef DEBUG14
		goffset,chroffset,chrhigh,watsonp,
#endif
		revp,lband,uband);
  Directions8_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			    rlength,glength,rsequence,gsequence,gsequence_alt,
#ifdef DEBUG14
			    goffset,chroffset,chrhigh,watsonp,
#endif
			    revp,lband,uband);
#endif
  
#ifdef DEBUG14
  matrix_std = compute_scores_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				       this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				       rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				       open,extend,lband,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_8);

#ifdef DEBUG2
  printf("Banded %s\n",revp ? "rev" : "fwd");
  Dynprog_Matrix32_print(matrix_std,rlength,glength,rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
  Dynprog_Directions32_print(directions_nogap_std,directions_Egap_std,directions_Fgap_std,
			     rlength,glength,rsequence,gsequence,gsequence_alt,
			     goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
#endif
  
  banded_matrix8_compare(matrix,matrix_std,rlength,glength,lband,uband,
			 rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions8_compare_nogap(matrix,*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);
  banded_directions8_compare_Egap(matrix,*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);
  banded_directions8_compare_Fgap(matrix,*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  FREE(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef HAVE_SSE2
/* Modified from Dynprog_simd_16_upper.  Operates by columns. */
Score16_T **
Dynprog_simd_16 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
		 Direction16_T ***directions_Fgap,
		 T this, char *rsequence, char *gsequence, char *gsequence_alt,
		 int rlength, int glength,
#ifdef DEBUG14
		 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		 Mismatchtype_T mismatchtype, int open, int extend,
		 int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF; /* Need to have the ability to go past NEG_INFINITY */
  Score16_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, all_one_bits;
  __m128i dir_horiz;
  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score16_T *pairscores_col0;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG14
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#endif

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned.one.matrix_ptrs,this->aligned.one.matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned.one.directions_ptrs_0,this->aligned.one.directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_1,this->aligned.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions16_calloc(rlength_ceil,glength,
						 this->aligned.one.directions_ptrs_2,this->aligned.one.directions_space_2);

  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NSHORTS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);

  lband_ceil = (int) ((lband + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairscores_col0 = (Score16_T *) _mm_malloc(lband_ceil * sizeof(Score16_T),16);


  /* pairscores_col0[0] = (Score16_T) 0; */
  /* Initialization just to lband causes errors in dir_horiz for Egap */
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r-1] = (Score16_T) NEG_INFINITY_16;
  }

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r-1] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }


  all_one_bits = _mm_set1_epi16(-1);
  
  FF = (int *) MALLOC((glength + 1) * sizeof(int));

  gap_open = _mm_set1_epi16((Score16_T) open);
  gap_extend = _mm_set1_epi16((Score16_T) extend);

  if (jump_late_p) {
    for (rlo = 1; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H .  To fill in first column of each
	 row block with diags, make E < H. */
      E_r_gap = _mm_set1_epi16(NEG_INFINITY_16);
      H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16+1);

      if ((c = rlo - lband) < 0) {
	c = 0;
      }
      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;
	} else {
	  na2 = revp ? gsequence[1-c] : gsequence[c-1];
	  na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
	  switch (na2) {
	  case 'A': pairscores_std_ptr = pairscores[0]; break;
	  case 'C': pairscores_std_ptr = pairscores[1]; break;
	  case 'G': pairscores_std_ptr = pairscores[2]; break;
	  case 'T': pairscores_std_ptr = pairscores[3]; break;
	  default: pairscores_std_ptr = pairscores[4];
	  }
	  switch (na2_alt) {
	  case 'A': pairscores_alt_ptr = pairscores[0]; break;
	  case 'C': pairscores_alt_ptr = pairscores[1]; break;
	  case 'G': pairscores_alt_ptr = pairscores[2]; break;
	  case 'T': pairscores_alt_ptr = pairscores[3]; break;
	  default: pairscores_alt_ptr = pairscores[4];
	  }
	}

	if (rlo == 1) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (c == 0) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi16(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) > c + lband) {
	  rhigh_calc = c + lband;
	}
	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 1) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = score_column[rlo_calc-1];
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for (r = rlo_calc; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }

  } else {
    /* jump early */
    for (rlo = 1; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with diags, make E <= H. */
      E_r_gap = _mm_set1_epi16(NEG_INFINITY_16);
      H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16+0);

      if ((c = rlo - lband) < 0) {
	c = 0;
      }
      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;
	} else {
	  na2 = revp ? gsequence[1-c] : gsequence[c-1];
	  na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
	  switch (na2) {
	  case 'A': pairscores_std_ptr = pairscores[0]; break;
	  case 'C': pairscores_std_ptr = pairscores[1]; break;
	  case 'G': pairscores_std_ptr = pairscores[2]; break;
	  case 'T': pairscores_std_ptr = pairscores[3]; break;
	  default: pairscores_std_ptr = pairscores[4];
	  }
	  switch (na2_alt) {
	  case 'A': pairscores_alt_ptr = pairscores[0]; break;
	  case 'C': pairscores_alt_ptr = pairscores[1]; break;
	  case 'G': pairscores_alt_ptr = pairscores[2]; break;
	  case 'T': pairscores_alt_ptr = pairscores[3]; break;
	  default: pairscores_alt_ptr = pairscores[4];
	  }
	}

	if (rlo == 1) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (c == 0) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi16(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo-1]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi16(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) > c + lband) {
	  rhigh_calc = c + lband;
	}
	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 1) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = score_column[rlo_calc-1];
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for (r = rlo_calc; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }
  }

#ifdef DEBUG2
  printf("SIMD\n");
  Matrix16_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
#ifdef DEBUG14
			 goffset,chroffset,chrhigh,watsonp,
#endif
			 revp,lband,uband);
  Directions16_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			     rlength,glength,rsequence,gsequence,gsequence_alt,
#ifdef DEBUG14
			     goffset,chroffset,chrhigh,watsonp,
#endif
			     revp,lband,uband);
#endif

#ifdef DEBUG14
  matrix_std = compute_scores_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				       this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				       rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				       open,extend,lband,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_16);

#ifdef DEBUG2
  printf("Banded\n");
  Dynprog_Matrix32_print(matrix_std,rlength,glength,rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
  Dynprog_Directions32_print(directions_nogap_std,directions_Egap_std,directions_Fgap_std,
			     rlength,glength,rsequence,gsequence,gsequence_alt,
			     goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
#endif
  
  banded_matrix16_compare(matrix,matrix_std,rlength,glength,lband,uband,
			  rsequence,gsequence,gsequence_alt,
			  goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions16_compare_nogap(*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);
  banded_directions16_compare_Egap(*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);
  banded_directions16_compare_Fgap(*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  FREE(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
  }
#endif



static List_T
push_token (bool *startp, List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);

  *startp = false;
  return List_push(tokens,(void *) copy);
}



#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
char *
Dynprog_cigar_8 (int *finalc, Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		 int r, int c, char *rsequence, char *gsequence, char *gsequence_alt,
		 char *nindels, int queryoffset, int genomeoffset, bool revp,
		 Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int insertion_width = 0;
  bool startp = true;

  int dist;
  Direction8_T dir;

  debug(printf("Starting Dynprog_cigar_8 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

#if 0
  /* Handle initial indel */
  if ((dir = directions_nogap[c][r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (c > 1 && directions_Egap[c][r] != DIAG) {
      dist++;
      c--;
    }
    c--;
    /* dir = directions_nogap[c][r]; */

    sprintf(token,"0M");
    tokens = push_token(&startp,tokens,token);
    sprintf(token,"%dD",dist);
    tokens = push_token(&startp,tokens,token);

  } else {
    /* Must be VERT */
    dist = 1;
    while (r > 1 && directions_Fgap[c][r] != DIAG) {
      dist++;
      r--;
    }
    r--;
    /* dir = directions_nogap[c][r]; */

    sprintf(token,"0M");
    tokens = push_token(&startp,tokens,token);
    sprintf(token,"%dI",dist);
    tokens = push_token(&startp,tokens,token);
  }
#endif

  *finalc = c;
  while (c > 0) {  /* dir != STOP */
    if (r == 0) {
      /* Ignore gap in row 0 */
      if (Ilength > 0 && Ilength != insertion_width) {
	if (Mlength_postins > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	}
	sprintf(token,"%dS",Ilength);
	tokens = push_token(&startp,tokens,token);
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      } else {
	if (Ilength > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(&startp,tokens,token);
	}
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      }

      *finalc = c;
      tokens = List_reverse(tokens);
      cigar = Dynprog_tokens_string(List_reverse(tokens));
      Dynprog_tokens_free(&tokens);

      return cigar;

    } else if ((dir = directions_nogap[c][r]) == DIAG) {
      /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */

      if (nindels[c-1] < 0) {
	/* Genome modified with a deletion */
	if (startp == true && Ilength > 0 && Ilength != insertion_width) {
	  /* Incomplete insertion at beginning.  Mlength_postins should be 0. */
	  sprintf(token,"%dS",Ilength);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	} else {
	  if (Ilength > 0) {
	    sprintf(token,"%dM",Mlength_postins);
	    tokens = push_token(&startp,tokens,token);
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(&startp,tokens,token);
	  }
	  if (1 || Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(&startp,tokens,token);
	  }
	}
	Mlength = 1;

	sprintf(token,"%dD",-nindels[c-1]);
	tokens = push_token(&startp,tokens,token);

      } else if (nindels[c-1] > 0) {
	insertion_width = nindels[c-1];
	if (Ilength == 0) {
	  Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	  Mlength = 0;
	}
	Ilength += 1;
	/* printf("Incrementing Ilength to be %d\n",Ilength); */

      } else {
	Mlength += 1;
	/* printf("Incrementing Mlength to be %d\n",Mlength); */
      }
      r--; c--;

    } else if (dir == HORIZ) {
      /* printf("At r %d, c %d, dir is HORIZ\n",r,c); */

      if (startp == true && Ilength > 0 && Ilength != insertion_width) {
	/* Incomplete insertion at beginning.  Mlength_postins should be 0. */
	sprintf(token,"%dS",Ilength);
	tokens = push_token(&startp,tokens,token);
	sprintf(token,"%dM",Mlength);
	tokens = push_token(&startp,tokens,token);
      } else {
	if (Ilength > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(&startp,tokens,token);
	}
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      }
      Mlength = 0;

      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dD",dist);
      tokens = push_token(&startp,tokens,token);

    } else {
      /* Must be VERT */
      /* printf("At r %d, c %d, dir is VERT\n",r,c); */
      if (Ilength == 0) {
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	  Mlength = 0;
	}
      }

      Ilength += 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	Ilength += 1;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dI",Ilength);
      tokens = push_token(&startp,tokens,token);
      Ilength = 0;
    }
  }

  if (Ilength > 0 && Ilength != insertion_width) {
    if (Mlength_postins > 0) {
      sprintf(token,"%dM",Mlength_postins);
      tokens = push_token(&startp,tokens,token);
    }
    sprintf(token,"%dS",Ilength);
    tokens = push_token(&startp,tokens,token);
    if (Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      tokens = push_token(&startp,tokens,token);
    }
  } else {
    if (Ilength > 0) {
      sprintf(token,"%dM",Mlength_postins);
      tokens = push_token(&startp,tokens,token);
      sprintf(token,"%dI",Ilength);
      tokens = push_token(&startp,tokens,token);
    }
    if (1 || Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      tokens = push_token(&startp,tokens,token);
    }
  }

  *finalc = c;
  tokens = List_reverse(tokens);
  cigar = Dynprog_tokens_string(List_reverse(tokens));
  Dynprog_tokens_free(&tokens);

  return cigar;
}
#endif



#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
char *
Dynprog_cigar_16 (int *finalc, Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		  int r, int c, char *rsequence, char *gsequence, char *gsequence_alt,
		  char *nindels, int queryoffset, int genomeoffset, bool revp,
		  Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int insertion_width = 0;
  bool startp = true;

  int dist;
  Direction16_T dir;

  debug(printf("Starting Dynprog_cigar_16 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

#if 0
  /* Handle initial indel */
  if ((dir = directions_nogap[c][r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (c > 1 && directions_Egap[c][r] != DIAG) {
      dist++;
      c--;
    }
    c--;
    /* dir = directions_nogap[c][r]; */

    sprintf(token,"0M");
    tokens = push_token(&startp,tokens,token);
    sprintf(token,"%dD",dist);
    tokens = push_token(&startp,tokens,token);

  } else {
    /* Must be VERT */
    dist = 1;
    while (r > 1 && directions_Fgap[c][r] != DIAG) {
      dist++;
      r--;
    }
    r--;
    /* dir = directions_nogap[c][r]; */

    sprintf(token,"0M");
    tokens = push_token(&startp,tokens,token);
    sprintf(token,"%dI",dist);
    tokens = push_token(&startp,tokens,token);
  }
#endif

  *finalc = c;
  while (c > 0) {  /* dir != STOP */
    if (r == 0) {
      /* Ignore gap in row 0 */
      if (Ilength > 0 && Ilength != insertion_width) {
	if (Mlength_postins > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	}
	sprintf(token,"%dS",Ilength);
	tokens = push_token(&startp,tokens,token);
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      } else {
	if (Ilength > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(&startp,tokens,token);
	}
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      }

      *finalc = c;
      tokens = List_reverse(tokens);
      cigar = Dynprog_tokens_string(List_reverse(tokens));
      Dynprog_tokens_free(&tokens);

      return cigar;

    } else if ((dir = directions_nogap[c][r]) == DIAG) {
      /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */

      if (nindels[c-1] < 0) {
	/* Genome modified with a deletion */
	if (startp == true && Ilength > 0 && Ilength != insertion_width) {
	  /* Incomplete insertion at beginning.  Mlength_postins should be 0. */
	  sprintf(token,"%dS",Ilength);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	} else {
	  if (Ilength > 0) {
	    sprintf(token,"%dM",Mlength_postins);
	    tokens = push_token(&startp,tokens,token);
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(&startp,tokens,token);
	  }
	  if (1 || Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(&startp,tokens,token);
	  }
	}
	Mlength = 1;

	sprintf(token,"%dD",-nindels[c-1]);
	tokens = push_token(&startp,tokens,token);

      } else if (nindels[c-1] > 0) {
	insertion_width = nindels[c-1];
	if (Ilength == 0) {
	  Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	  Mlength = 0;
	}
	Ilength += 1;
	/* printf("Incrementing Ilength to be %d\n",Ilength); */

      } else {
	Mlength += 1;
	/* printf("Incrementing Mlength to be %d\n",Mlength); */
      }
      r--; c--;

    } else if (dir == HORIZ) {
      /* printf("At r %d, c %d, dir is HORIZ\n",r,c); */

      if (startp == true && Ilength > 0 && Ilength != insertion_width) {
	/* Incomplete insertion at beginning.  Mlength_postins should be 0. */
	sprintf(token,"%dS",Ilength);
	tokens = push_token(&startp,tokens,token);
	sprintf(token,"%dM",Mlength);
	tokens = push_token(&startp,tokens,token);
      } else {
	if (Ilength > 0) {
	  sprintf(token,"%dM",Mlength_postins);
	  tokens = push_token(&startp,tokens,token);
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(&startp,tokens,token);
	}
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	}
      }
      Mlength = 0;

      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dD",dist);
      tokens = push_token(&startp,tokens,token);

    } else {
      /* Must be VERT */
      /* printf("At r %d, c %d, dir is VERT\n",r,c); */
      if (Ilength == 0) {
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(&startp,tokens,token);
	  Mlength = 0;
	}
      }

      Ilength += 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	Ilength += 1;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dI",Ilength);
      tokens = push_token(&startp,tokens,token);
      Ilength = 0;
    }
  }

  if (Ilength > 0 && Ilength != insertion_width) {
    if (Mlength_postins > 0) {
      sprintf(token,"%dM",Mlength_postins);
      tokens = push_token(&startp,tokens,token);
    }
    sprintf(token,"%dS",Ilength);
    tokens = push_token(&startp,tokens,token);
    if (Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      tokens = push_token(&startp,tokens,token);
    }
  } else {
    if (Ilength > 0) {
      sprintf(token,"%dM",Mlength_postins);
      tokens = push_token(&startp,tokens,token);
      sprintf(token,"%dI",Ilength);
      tokens = push_token(&startp,tokens,token);
    }
    if (1 || Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      tokens = push_token(&startp,tokens,token);
    }
  }

  *finalc = c;
  tokens = List_reverse(tokens);
  cigar = Dynprog_tokens_string(List_reverse(tokens));
  Dynprog_tokens_free(&tokens);

  return cigar;
}
#endif


