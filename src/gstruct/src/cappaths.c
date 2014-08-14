static char rcsid[] = "$Id: cappaths.c 136513 2014-05-16 17:58:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For isdigit */
#include <unistd.h>		/* For getopt */
#include <math.h>		/* For sqrt */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include "mem.h"
#include "fopen.h"

#include "datadir.h"
#include "parserange.h"
#include "genome.h"
#include "iit-read.h"
#include "list.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* trace rightward and leftward by loglik */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif



static double slope_threshold_utr = 0.20; /* Looser criterion */
static int flat_length_utr = 200;

/************************************************************************
 *   Signal processing functions
 ************************************************************************/


#define HALFDATA 10		/* ramp size */
static double xcoef[HALFDATA+1+HALFDATA];
static double linear_meanx;
static double linear_den;

void
Cappaths_setup () {
  double sumx = 0.0, sumxx = 0.0;
  int i, j;

  /* Want x_up to be 0:(NDATA-1), and x_down to be -(NDATA-1):0 */
  for (i = 0, j = -HALFDATA; j <= HALFDATA; i++, j++) {
    xcoef[i] = (double) j;
    sumx += xcoef[i];
    sumxx += xcoef[i]*xcoef[i];
  }

  linear_meanx = sumx/(double) (HALFDATA+1+HALFDATA);

  linear_den = sumxx - sumx*sumx/(double) (HALFDATA+1+HALFDATA);

  return;
}


#if 0
static double
linear_slope (long int *heights, double *xcoef, double meanx, double den) {
  double slope;
  long int sumy;
  double sumxy;
  int i;

  sumy = 0;
  sumxy = 0.0;
  for (i = 0; i < HALFDATA+1+HALFDATA; i++) {
    sumy += heights[i];
    sumxy += xcoef[i] * (double) heights[i];
  }
  if (sumy == 0) {
    return 0.0;	
  } else {
    slope = sumxy - sumy * meanx;
    slope /= den;
    return slope;
  }
}
#endif


static double
linear_xintercept (long int *heights1, long int *heights2, double *xcoef, double meanx) {
  double xintercept,  yintercept, slope;
  long int sumy;
  double meany, sumxy;
  int i;

  sumy = 0;
  for (i = 0; i < HALFDATA+1+HALFDATA; i++) {
    sumy += heights1[i] + heights2[i];
  }
  if (sumy == 0) {
    return 0.0;			/* Terminates search for a low x-intercept */
  } else {
    meany = (double) sumy/(double) (HALFDATA+1+HALFDATA);
  }

  sumxy = 0.0;
  for (i = 0; i < HALFDATA+1+HALFDATA; i++) {
    sumxy += xcoef[i] * (double) (heights1[i] + heights2[i]);
  }

  slope = sumxy - sumy * meanx;
  slope /= linear_den;

  yintercept = meany - slope * meanx;

  if (slope == 0.0) {
    return 1000000.0;
  } else {
    xintercept = -yintercept/slope;
    return xintercept;
  }
}



#if 0
static double *
compute_slopes (long int *counts_pairing, Genomicpos_T chrlength) {
  double *slopes;
  Genomicpos_T pos;

  slopes = (double *) CALLOC(chrlength+1,sizeof(double));
  for (pos = HALFDATA; pos + HALFDATA < chrlength; pos++) {
    slopes[pos] = linear_slope(&(counts_pairing[pos-HALFDATA]),xcoef,linear_meanx,linear_den);
  }
  return slopes;
}
#endif


#if 0
static double *
compute_xintercepts (long int *counts_tally, Genomicpos_T chrlength) {
  double *xintercepts;
  Genomicpos_T pos;

  xintercepts = (double *) CALLOC(chrlength+1,sizeof(double));
  for (pos = HALFDATA; pos + HALFDATA < chrlength; pos++) {
    xintercepts[pos] = linear_xintercept(&(counts_tally[pos-HALFDATA]),xcoef,linear_meanx);
  }
  return xintercepts;
}
#endif

static double
compute_xintercept (Genomicpos_T pos, long int *tally_matches, long int *tally_mismatches,
		    Genomicpos_T chrlength) {
  if (pos < HALFDATA) {
    return 0.0;
  } else if (pos + HALFDATA >= chrlength) {
    return 0.0;
  } else {
    return linear_xintercept(&(tally_matches[pos-HALFDATA]),&(tally_mismatches[pos-HALFDATA]),
			     xcoef,linear_meanx);
  }
}

/************************************************************************/


#if 0
static Genomicpos_T
utr_zone_right_precomputed (double *xintercepts, Genomicpos_T splicepos, Genomicpos_T chrlength) {
  Genomicpos_T bestpos = splicepos, pos;
  double min_xintercept = 1000000.0;

  for (pos = splicepos + 20; pos < chrlength && pos < splicepos + 5000; pos++) {
    /* printf("%u %.2f vs %.2f",pos,xintercepts[pos],min_xintercept); */
    if (fabs(xintercepts[pos]) < min_xintercept) {
      /* printf(" **"); */
      min_xintercept = fabs(xintercepts[pos]);
      bestpos = pos;
    }
    /* printf("\n"); */
  }

  return bestpos;
}


static Genomicpos_T
utr_zone_left_precomputed (double *xintercepts, Genomicpos_T splicepos) {
  Genomicpos_T bestpos = splicepos, low, start, pos;
  double min_xintercept = 1000000.0;

  if (splicepos <= 20) {
    return splicepos;
  }

  if (splicepos < 5000) {
    low = 0;
  } else {
    low = splicepos - 5000;
  }
  for (pos = splicepos - 20; pos > low; pos--) {
    if (fabs(xintercepts[pos]) < min_xintercept) {
      min_xintercept = fabs(xintercepts[pos]);
      bestpos = pos;
    }
  }

  return bestpos;
}
#endif


static Genomicpos_T
utr_zone_right (long int *tally_matches, long int *tally_mismatches,
		Genomicpos_T splicepos, Genomicpos_T chrlength) {
  Genomicpos_T bestpos = splicepos, pos;
  double min_xintercept = 1000000.0, intercept;

  for (pos = splicepos + 20; pos < chrlength && pos < splicepos + 5000; pos++) {
    /* printf("%u %.2f vs %.2f",pos,xintercepts[pos],min_xintercept); */
    if ((intercept = fabs(compute_xintercept(pos,tally_matches,tally_mismatches,chrlength))) < min_xintercept) {
      /* printf(" **"); */
      min_xintercept = intercept;
      bestpos = pos;
    }
    /* printf("\n"); */
  }

  if (bestpos >= chrlength) {
    return chrlength - 1;
  } else {
    return bestpos;
  }
}


static Genomicpos_T
utr_zone_left (long int *tally_matches, long int *tally_mismatches,
	       Genomicpos_T splicepos, Genomicpos_T chrlength) {
  Genomicpos_T bestpos = splicepos, low, pos;
  double min_xintercept = 1000000.0, intercept;

  if (splicepos <= 20) {
    return splicepos;
  }

  if (splicepos < 5000) {
    low = 0;
  } else {
    low = splicepos - 5000;
  }
  for (pos = splicepos - 20; pos > low; pos--) {
    if ((intercept = fabs(compute_xintercept(pos,tally_matches,tally_mismatches,chrlength))) < min_xintercept) {
      min_xintercept = intercept;
      bestpos = pos;
    }
  }

  if (bestpos >= chrlength) {
    return chrlength - 1;
  } else {
    return bestpos;
  }
}


/************************************************************************
 *   Input
 ************************************************************************/

/* Eventually want to add extents_null */
Genomicpos_T
Cappaths_solve_genestart (Genomicpos_T pathstart,
			  long int *tally_matches, long int *tally_mismatches,
			  IIT_T end_exons_iit, char *chr, Genomicpos_T chrlength,
			  bool forwardp) {
  int divno;
  int *matches, nmatches, i;
  Interval_T interval;
  Genomicpos_T farthest;

  if (end_exons_iit != NULL) {
    farthest = pathstart;
    divno = IIT_divint(end_exons_iit,chr);
    if (forwardp == true) {
      matches = IIT_get_signed_with_divno(&nmatches,end_exons_iit,divno,
					  /*x*/pathstart,/*y*/pathstart,/*sortp*/false,/*sign*/+1);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(end_exons_iit,matches[i]);
	if (Interval_high(interval) == pathstart) {
	  if (Interval_low(interval) < farthest) {
	    farthest = Interval_low(interval);
	  }
	}
      }

    } else {
      matches = IIT_get_signed_with_divno(&nmatches,end_exons_iit,divno,
					  /*x*/pathstart,/*y*/pathstart,/*sortp*/false,/*sign*/-1);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(end_exons_iit,matches[i]);
	if (Interval_low(interval) == pathstart) {
	  if (Interval_high(interval) > farthest) {
	    farthest = Interval_high(interval);
	  }
	}
      }
    }

    FREE(matches);
    return farthest;

  } else if (forwardp == true) {
    return utr_zone_left(tally_matches,tally_mismatches,pathstart,chrlength);
  } else {
    return utr_zone_right(tally_matches,tally_mismatches,pathstart,chrlength);
  }
}

Genomicpos_T
Cappaths_solve_geneend (Genomicpos_T pathend,
			long int *tally_matches, long int *tally_mismatches,
			IIT_T end_exons_iit, char *chr, Genomicpos_T chrlength,
			bool forwardp) {
  int divno;
  int *matches, nmatches, i;
  Interval_T interval;
  Genomicpos_T farthest;

  if (end_exons_iit != NULL) {
    farthest = pathend;
    divno = IIT_divint(end_exons_iit,chr);
    if (forwardp == true) {
      matches = IIT_get_signed_with_divno(&nmatches,end_exons_iit,divno,
					  /*x*/pathend,/*y*/pathend,/*sortp*/false,/*sign*/+1);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(end_exons_iit,matches[i]);
	if (Interval_low(interval) == pathend) {
	  if (Interval_high(interval) > farthest) {
	    farthest = Interval_high(interval);
	  }
	}
      }

    } else {
      matches = IIT_get_signed_with_divno(&nmatches,end_exons_iit,divno,
					  /*x*/pathend,/*y*/pathend,/*sortp*/false,/*sign*/-1);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(end_exons_iit,matches[i]);
	if (Interval_high(interval) == pathend) {
	  if (Interval_low(interval) < farthest) {
	    farthest = Interval_low(interval);
	  }
	}
      }
    }

    FREE(matches);
    return farthest;

  } else if (forwardp == true) {
    return utr_zone_right(tally_matches,tally_mismatches,pathend,chrlength);
  } else {
    return utr_zone_left(tally_matches,tally_mismatches,pathend,chrlength);
  }
}

