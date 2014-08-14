static char rcsid[] = "$Id: tally.old.c 136513 2014-05-16 17:58:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tally.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For index */
#include <ctype.h>		/* For isdigit */
#include <math.h>		/* For log, sqrt */

#include "mem.h"
#include "interval.h"
#include "orderstat.h"
#include "lgamma.h"


/* store_counts */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Median filtering */
#ifdef DEBUGM
#define debugM(x) x
#else
#define debugM(x)
#endif

/* tracing */
#ifdef DEBUGT
#define debugT(x) x
#else
#define debugT(x)
#endif

/* loglikelihood */
#ifdef DEBUGL
#define debugL(x) x
#else
#define debugL(x)
#endif



static long int
get_total_tally (const char *ptr) {
  long int tally = 0, n;
  char *end;

  if ((end = index(ptr,'\n')) == NULL) {
    fprintf(stderr,"Premature end of line %s\n",ptr);
    return 0;
  }
  /* fprintf(stderr,"Getting tally for %.*s\n",end-ptr,ptr); */

  while (ptr < end) {
    while (ptr < end && !isdigit((int) *ptr)) {
      ptr++;
    }
    if (ptr < end) {
      sscanf(ptr,"%ld",&n);
      tally += n;
      while (ptr < end && !isspace(*ptr)) {
	ptr++;
      }
      while (ptr < end && isspace(*ptr)) {
	ptr++;
      }
    }
  }

  return tally;
}


static void
store_counts_aux (long int *maxcount, long int *counts_tally,
#if 0
		  long int *counts_thresholded,
#endif
		  unsigned int coordstart, unsigned int coordend, 
		  int indexi, IIT_T iit) {
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
  bool allocp;
  unsigned int chrpos, intervalend;
  long int count;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      /* Can happen because we retrieved intervals for exon and intron */
      if (allocp == true) {
	FREE(restofheader);
      }
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }
    
  while (chrpos <= intervalend && chrpos <= coordend) {
    count = get_total_tally(ptr);
    if (count > *maxcount) {
      *maxcount = count;
    }
    counts_tally[chrpos] = count;
#if 0
    if (count > XMAX) {
      counts_thresholded[chrpos] = XMAX;
    } else {
      counts_thresholded[chrpos] = count;
    }
#endif

    debug1(printf("Putting %ld into chrpos %u\n",counts_tally[chrpos],chrpos));

    if ((ptr = index(ptr,'\n')) == NULL) {
      if (allocp == true) {
	FREE(restofheader);
      }
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }

  if (allocp == true) {
    FREE(restofheader);
  }

  return;
}

long int
Tally_store_counts (long int *counts_tally, IIT_T tally_iit,
		    char *chr, unsigned int coordstart, unsigned int coordend) {
  long int maxcount = 0;
  int *matches;
  int nmatches, i;

  matches = IIT_get(&nmatches,tally_iit,/*divstring*/chr,coordstart,coordend,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    store_counts_aux(&maxcount,counts_tally,
#if 0
		     counts_thresholded,
#endif
		     coordstart,coordend,matches[i],tally_iit);
  }
  
  FREE(matches);

  return maxcount;
}


static void
store_runlength_aux (long int *maxcount, long int *counts,
		     unsigned int coordstart, unsigned int coordend, 
		     int indexi, IIT_T iit) {
  Interval_T interval;
  char *label;
  bool allocp;
  unsigned int chrpos, intervalend;
  long int count;

  label = IIT_label(iit,indexi,&allocp);
  count = strtoul(label,NULL,10);

  if (count > 0) {
    if (count > *maxcount) {
      *maxcount = count;
    }

    interval = IIT_interval(iit,indexi);
    chrpos = Interval_low(interval);
    intervalend = Interval_high(interval);

    /* Advance if necessary */
    while (chrpos < coordstart) {
      chrpos++;
    }
    
    while (chrpos <= intervalend && chrpos <= coordend) {
      /* printf("Storing runlength %d at chrpos %u\n",count,chrpos); */
      counts[chrpos] = count;
      chrpos++;
    }
  }

  if (allocp == true) {
    FREE(label);
  }

  return;
}



long int
Tally_store_runlength (long int *counts, IIT_T runlength_iit,
		       char *chr, unsigned int coordstart, unsigned int coordend) {
  long int maxcount = 0;
  int *matches;
  int nmatches, i;

  matches = IIT_get(&nmatches,runlength_iit,/*divstring*/chr,coordstart,coordend,/*sortp*/false);

  /* printf("Got %d matches for %u..%u\n",nmatches,coordstart,coordend); */
  for (i = 0; i < nmatches; i++) {
    /* printf("Index %d\n",matches[i]); */
    store_runlength_aux(&maxcount,counts,
			coordstart,coordend,matches[i],runlength_iit);
  }
  
  FREE(matches);

  return maxcount;
}



static void
add_runlength_aux (long int *counts,
		   unsigned int coordstart, unsigned int coordend, 
		   int indexi, IIT_T iit) {
  Interval_T interval;
  char *label;
  bool allocp;
  unsigned int chrpos, intervalend;
  long int count;

  label = IIT_label(iit,indexi,&allocp);
  count = strtoul(label,NULL,10);

  if (count > 0) {
    interval = IIT_interval(iit,indexi);
    chrpos = Interval_low(interval);
    intervalend = Interval_high(interval);

    /* Advance if necessary */
    while (chrpos < coordstart) {
      chrpos++;
    }
    
    while (chrpos <= intervalend && chrpos <= coordend) {
      /* printf("Adding runlength %d at chrpos %u\n",count,chrpos); */
      counts[chrpos] += count;
      chrpos++;
    }
  }

  if (allocp == true) {
    FREE(label);
  }

  return;
}



void
Tally_add_runlength (long int *counts, IIT_T runlength_iit,
		     char *chr, unsigned int coordstart, unsigned int coordend) {
  int *matches;
  int nmatches, i;

  matches = IIT_get(&nmatches,runlength_iit,/*divstring*/chr,coordstart,coordend,/*sortp*/false);

  /* printf("Got %d matches for %u..%u\n",nmatches,coordstart,coordend); */
  for (i = 0; i < nmatches; i++) {
    /* printf("Index %d\n",matches[i]); */
    add_runlength_aux(counts,coordstart,coordend,matches[i],runlength_iit);
  }
  
  FREE(matches);

  return;
}



void
Tally_compute_log (double *log_tally, int *counts_tally, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;

  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    log_tally[chrpos] = log((double) (counts_tally[chrpos]+1));
  }
  
  return;
}

long int
Tally_sum (long int *maxcount, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  long int sum = 0;
  Genomicpos_T chrpos;

  *maxcount = 0;
  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    /* printf("%u %d %d\n",chrpos,x[chrpos],sum); */
    sum = x[chrpos] + sum;
    if (x[chrpos] > *maxcount) {
      *maxcount = x[chrpos];
    }
  }
  return sum;
}

double
Tally_mean (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  long int sum = 0;
  Genomicpos_T chrpos;

  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    /* printf("%u %d %d\n",chrpos,x[chrpos],sum); */
    sum = x[chrpos] + sum;
  }

  return (double) sum/(double) (chrend - chrstart + 1);
}

double
Tally_mean_double (double *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  double sum = 0.0;
  Genomicpos_T chrpos;

  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    /* printf("%u %d %d\n",chrpos,x[chrpos],sum); */
    sum = x[chrpos] + sum;
  }

  return sum/(double) (chrend - chrstart + 1);
}


long int
Tally_median (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  return Orderstat_long_int_pct(&(x[chrstart]),chrend-chrstart+1,/*percentile*/0.50);
}

double
Tally_median_double (double *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  return Orderstat_double_pct(&(x[chrstart]),chrend-chrstart+1,/*percentile*/0.50);
}


long int
Tally_quantile (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend, double percentile) {
  return Orderstat_long_int_pct(&(x[chrstart]),/*n*/chrend-chrstart+1,percentile);
}


void
Tally_stats (long int *minx, long int *maxx, double *mean, double *sdev,
	     long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;
  int sumx = 0, sumxx = 0, n;
  double var;

  *minx = *maxx = x[chrstart];

  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    if (x[chrpos] < *minx) {
      *minx = x[chrpos];
    }
    if (x[chrpos] > *maxx) {
      *maxx = x[chrpos];
    }

    sumx += x[chrpos];
    sumxx += x[chrpos] * x[chrpos];
  }

  n = (int) (chrend - chrstart + 1);
  *mean = (double) sumx/(double) n;
  var = ((double) sumxx - (*mean) * (double) sumx)/(double) (n - 1);
  *sdev = sqrt(var);

  return;
}


void
Tally_range (long int *mincount, long int *maxcount, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;

  *mincount = *maxcount = x[chrstart];
  for (chrpos = chrstart+1; chrpos <= chrend; chrpos++) {
    if (x[chrpos] < *mincount) {
      *mincount = x[chrpos];
    }
    if (x[chrpos] > *maxcount) {
      *maxcount = x[chrpos];
    }
  }
  return;
}

long int
Tally_maxcount (long int *x, Genomicpos_T chrlength, Genomicpos_T chrstart, Genomicpos_T chrend) {
  long int maxcount = 0;
  Genomicpos_T chrpos;

  if (chrend > chrlength) {
    chrend = chrlength;
  }

  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    if (x[chrpos] > maxcount) {
      maxcount = x[chrpos];
    }
  }
  return maxcount;
}





void
Tally_cumulate_int (long int *cum, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;

  cum[(chrstart - 1U)] = 0;
  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    cum[chrpos] = x[chrpos] + cum[(chrpos - 1U)];
  }
  return;
}

void
Tally_cumulate_int_to_double (double *cum, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;

  cum[chrstart - 1U] = 0.0;
  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    cum[chrpos] = (double) x[chrpos] + cum[chrpos - 1U];
  }
  return;
}

void
Tally_cumulate_double (double *cum, double *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
  Genomicpos_T chrpos;

  cum[(chrstart - 1U)] = 0.0;
  for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
    cum[chrpos] = x[chrpos] + cum[(chrpos - 1U)];
  }
  return;
}


/************************************************************************
 *   Intron test.  Sees if edges (exon) are greater than some low percentile of middle (intron)
 ************************************************************************/

bool
Tally_introntest (long int *counts, Genomicpos_T low, Genomicpos_T high,
		  double percentile, int testregion) {
  long int pct;
  Genomicpos_T pos;

  pct = Orderstat_long_int_pct(&(counts[low]),high-low+1,percentile);

  /* Test below intron */
  for (pos = low-1; pos >= low - testregion; pos--) {
    if (counts[pos] <= pct) {
      return false;
    }
  }

  /* Test above intron */
  for (pos = high+1; pos <= high + testregion; pos++) {
    if (counts[pos] <= pct) {
      return false;
    }
  }
  
  return true;
}


/************************************************************************
 *   Exon test.  Sees if edges (intron) are less than some low percentile of middle (exon)
 ************************************************************************/

bool
Tally_exontest_leftward (long int *counts, Genomicpos_T exonlow, Genomicpos_T exonhigh,
			 double percentile, int testregion) {
  long int pct;
  Genomicpos_T pos;

  pct = Orderstat_long_int_pct(&(counts[exonlow]),exonhigh-exonlow+1,percentile);

  /* Test below exon */
  for (pos = exonlow-1; pos >= exonlow - testregion; pos--) {
    if (counts[pos] >= pct) {
      return false;
    }
  }

#if 0
  /* Test above exon */
  for (pos = exonhigh+1; pos <= exonhigh + testregion; pos++) {
    if (counts[pos] >= pct) {
      return false;
    }
  }
#endif
  
  return true;
}


bool
Tally_exontest_rightward (long int *counts, Genomicpos_T exonlow, Genomicpos_T exonhigh,
			  double percentile, int testregion) {
  long int pct;
  Genomicpos_T pos;

  pct = Orderstat_long_int_pct(&(counts[exonlow]),exonhigh-exonlow+1,percentile);

#if 0
  /* Test below exon */
  for (pos = exonlow-1; pos >= exonlow - testregion; pos--) {
    if (counts[pos] >= pct) {
      return false;
    }
  }
#endif

  /* Test above exon */
  for (pos = exonhigh+1; pos <= exonhigh + testregion; pos++) {
    if (counts[pos] >= pct) {
      return false;
    }
  }
  
  return true;
}


static bool
exontest_left (long int *counts, Genomicpos_T exonlow, int pct, int testregion) {
  Genomicpos_T pos;

  /* Test below exon */
  for (pos = exonlow-1; pos >= exonlow - testregion; pos--) {
    if (counts[pos] >= pct) {
      return false;
    }
  }

  return true;
}


static bool
exontest_right (long int *counts, Genomicpos_T exonhigh, long int pct, int testregion) {
  Genomicpos_T pos;

  /* Test above exon */
  for (pos = exonhigh+1; pos <= exonhigh + testregion; pos++) {
    /* printf("  At position %d, got count %d >= %d\n",pos,counts[pos],pct); */
    if (counts[pos] >= pct) {
      return false;
    }
  }

  return true;
}



Genomicpos_T
Tally_trace_leftward_exontest (long int *counts, Genomicpos_T exonhigh, int min_exonlength, int max_exonlength,
			       double percentile, int testregion) {
  Genomicpos_T exonlow;
  long int pct;

  if (exonhigh < min_exonlength) {
    return 0;
  } else {
    exonlow = exonhigh - min_exonlength;
  }

  while (exonlow + max_exonlength >= exonhigh) {
    pct = Orderstat_long_int_pct(&(counts[exonlow]),exonhigh-exonlow+1,percentile);
    /* printf("Percentile %.2f for %u..%u is %d\n",percentile,exonlow,exonhigh,pct); */
    if (exontest_left(counts,exonlow,pct,testregion) == true) {
      return exonlow;
    }
    exonlow--;
  }

  return exonlow;
}


Genomicpos_T
Tally_trace_rightward_exontest (long int *counts, Genomicpos_T exonlow, int min_exonlength, int max_exonlength,
				double percentile, int testregion) {
  Genomicpos_T exonhigh;
  long int pct;

  for (exonhigh = exonlow + min_exonlength; exonhigh <= exonlow + max_exonlength; exonhigh++) {
    pct = Orderstat_long_int_pct(&(counts[exonlow]),exonhigh-exonlow+1,percentile);
    /* printf("Percentile %.2f for %u..%u is %d\n",percentile,exonlow,exonhigh,pct); */
    if (exontest_right(counts,exonhigh,pct,testregion) == true) {
      return exonhigh;
    }
  }

  return exonhigh;
}



Genomicpos_T
Tally_trace_leftward_pairing (long int *counts, Genomicpos_T low, Genomicpos_T high,
			      int min_exonlength, int max_exonlength,
			      double foldchange_down, double foldchange_up, int testregion) {
  Genomicpos_T bound;
  double mean;
  long int sum, mincount, maxcount;

  sum = Tally_sum(&maxcount,counts,low,high);
  mean = sum/(double) (high - low);

  if (low < min_exonlength) {
    return 0;
  } else {
    bound = low - min_exonlength;
  }

  while (bound + max_exonlength >= low) {
    Tally_range(&mincount,&maxcount,counts,bound-testregion,bound);
    if (maxcount < mean/foldchange_down || mincount > mean*foldchange_up) {
      return bound;
    }
    bound--;
  }

  return bound;
}


Genomicpos_T
Tally_trace_rightward_pairing (long int *counts, Genomicpos_T low, Genomicpos_T high,
			       int min_exonlength, int max_exonlength,
			       double foldchange_down, double foldchange_up, int testregion) {
  Genomicpos_T bound;
  double mean;
  long int sum, mincount, maxcount;

  sum = Tally_sum(&maxcount,counts,low,high);
  mean = sum/(double) (high - low);

  bound = high + min_exonlength;
  while (bound <= high + max_exonlength) {
    Tally_range(&mincount,&maxcount,counts,bound,bound+testregion);
    if (maxcount < mean/foldchange_down || mincount > mean*foldchange_up) {
      return bound;
    }
    bound++;
  }

  return bound;
}


/************************************************************************
 *   Gap test.  Sees if a region exists that is significantly below the endpoints
 ************************************************************************/

bool
Tally_gaptest (long int *counts, Genomicpos_T low, Genomicpos_T high,
	       int testregion) {
  long int ceiling, maxcount, minlevel;
  Genomicpos_T pos;
  
  /* printf("Starting gaptest from %u to %u\n",low,high); */
  if (counts[low] < counts[high]) {
    minlevel = counts[low];
  } else {
    minlevel = counts[high];
  }
  if (minlevel < 10) {
    ceiling = 0;
  } else {
    ceiling = minlevel / 10;
  }
  /* printf("  minlevel = %d, ceiling = %d\n",minlevel,ceiling); */

  for (pos = low; pos < high - testregion; pos++) {
    maxcount = Tally_maxcount(counts,/*chrlength*/-1U,pos,pos+testregion);
    if (maxcount <= ceiling) {
      /* printf("  maxcount at position %u is %d, so returning true\n",pos,maxcount); */
      return true;
    }
  }

  /* printf("  returning false, no gap\n"); */
  return false;
}



/************************************************************************
 *   Median filtering
 ************************************************************************/

static long int
shift_up (long int *histogram, int *k1, int *k2, long int last_pct) {

  if (*k2 > 0) {
    (*k2)--;
    *k1 = histogram[last_pct] - (*k2) - 1;
    debugM(printf("shift up by adjusting k1 and k2\n"));
    return last_pct;

  } else {
    last_pct++;
    while (histogram[last_pct] == 0) {
      last_pct++;
    }
    
    *k1 = 0;
    *k2 = histogram[last_pct] - 1;
    
    debugM(printf("shift up to %d\n",last_pct));
    return last_pct;
  }
}

static long int
shift_down (long int *histogram, int *k1, int *k2, long int last_pct) {
  
  if (*k1 > 0) {
    (*k1)--;
    (*k2) = histogram[last_pct] - (*k1) - 1;
    debugM(printf("shift down by adjusting k1 and k2\n"));
    return last_pct;

  } else {
    last_pct--;
    while (histogram[last_pct] == 0) {
      last_pct--;
    }

    *k1 = histogram[last_pct] - 1;
    *k2 = 0;

    debugM(printf("shift down to %d\n",last_pct));
    return last_pct;
  }
}


#ifdef DEBUGM
static void
print_histogram (long int *histogram, int k1, int k2, long int last_pct, long int maxcount) {
  long int count;

  for (count = 0; count < maxcount; count++) {
    if (histogram[count] > 0) {
      printf(" %d:%d",count,histogram[count]);
      if (count == last_pct) {
	printf(" (%d ^ %d)",k1,k2);
      }
    }
  }
  printf("\n");
  return;
}
#endif


/* For a sliding window */
static long int
get_pct (long int last_pct, long int *histogram, int *k1, int *k2, long int oldcount, long int newcount) {
  
  histogram[oldcount] -= 1;
  histogram[newcount] += 1;

  debugM(printf("last_pct:%d, oldcount:%d, newcount:%d, k1:%d, k2:%d ",
		last_pct,oldcount,newcount,*k1,*k2));
  if (oldcount < last_pct) {
    if (newcount < last_pct) {
      /* No change */
      debugM(printf("no change\n"));
      return last_pct;

    } else if (newcount == last_pct) {
      /* Record one at last_pct below current pointer */
      (*k1)++;
      debugM(printf("increment k1 to %d\n",*k1));
      return last_pct;

    } else if (newcount > last_pct) {
      return shift_up(histogram,&(*k1),&(*k2),last_pct);
    }

  } else if (oldcount == last_pct) {
    if (newcount < last_pct) {
      return shift_down(histogram,&(*k1),&(*k2),last_pct);

    } else if (newcount == last_pct) {
      /* No change */
      debugM(printf("no change\n"));
      return last_pct;

    } else if (newcount > last_pct) {
      return shift_up(histogram,&(*k1),&(*k2),last_pct);
    }
    
  } else if (oldcount > last_pct) {
    if (newcount < last_pct) {
      return shift_down(histogram,&(*k1),&(*k2),last_pct);

    } else if (newcount == last_pct) {
      (*k2)++;
      debugM(printf("increment k2 to %d\n",*k2));
      return last_pct;
      
    } else if (newcount > last_pct) {
      /* No change */
      debugM(printf("no change\n"));
      return last_pct;
    }
  }

  abort();
  return 0;
}

      
static long int
initialize_histogram (long int *histogram, int *k1, int *k2, long int *data, int width,
		      double percentile) {
  int pct, desired_nbelow, nbelow, count;
  int i;

  for (i = 0; i < width; i++) {
    histogram[data[i]] += 1;
  }

  pct = Orderstat_long_int_pct(data,width,percentile);
  if (histogram[pct] == 1) {
    *k1 = *k2 = 0;
    return pct;

  } else {
    desired_nbelow = (int) (percentile * (double) width);

    nbelow = 0;
    for (count = pct-1; count >= 0; count--) {
      nbelow += histogram[count];
    }

    *k1 = desired_nbelow - nbelow;
    *k2 = histogram[pct] - (*k1) - 1;
    
    debugM(printf("At pct %d, k1 %d = desired %d - nbelow %d, k2 %d = hist %d - k1\n",
		  pct,*k1,desired_nbelow,nbelow,*k2,histogram[pct]));
    return pct;
  }
}


void
Tally_median_filter (long int *new_tally, long int *raw_tally, long int maxcount, Genomicpos_T median_halfwidth,
		     Genomicpos_T chrstart, Genomicpos_T chrend, Genomicpos_T chrlength) {
  long int *histogram, *data;
  int count, width;
  Genomicpos_T chrpos;
  long int last_pct;
  int k1, k2;

  histogram = (long int *) CALLOC(maxcount+1,sizeof(long int));
  for (count = 0; count < maxcount; count++) {
    histogram[count] = 0;
  }

  fprintf(stderr,"   Performing median filtering...");
  if (chrend - chrstart < median_halfwidth + median_halfwidth) {
    median_halfwidth = (chrend - chrstart)/2;
  }

  /* Filtering near left end of chromosome */
  for (chrpos = chrstart, width = median_halfwidth; chrpos < chrstart + median_halfwidth;
       chrpos++, width++) {
    data = &(raw_tally[chrstart]);
    new_tally[chrpos] = Orderstat_long_int_pct(data,width,/*percentile*/0.50);
  }

  /* Percentile in middle */
  last_pct = initialize_histogram(histogram,&k1,&k2,&(raw_tally[chrpos-median_halfwidth]),
				  /*width*/median_halfwidth+median_halfwidth+1,/*percentile*/0.50);
  new_tally[chrpos] = last_pct;

  debugM(printf("Initial histogram:\n"));
  debugM(print_histogram(histogram,k1,k2,last_pct,maxcount));

  for (chrpos++; chrpos < chrend - median_halfwidth; chrpos++) {
    debugM(printf("pct should be %d  ",Orderstat_long_int_pct(&(raw_tally[chrpos-median_halfwidth]),
							 median_halfwidth+1+median_halfwidth,/*percentile*/0.50)));
    debugM(print_histogram(histogram,k1,k2,last_pct,maxcount));

    last_pct = get_pct(last_pct,histogram,&k1,&k2,raw_tally[chrpos-median_halfwidth-1],
		       raw_tally[chrpos+median_halfwidth]);
    new_tally[chrpos] = last_pct;
  }

  /* Percentile near right end of chromosome */
  for (width--; chrpos <= chrend; chrpos++, width--) {
    data = &(raw_tally[chrpos - median_halfwidth]);
    new_tally[chrpos] = Orderstat_long_int_pct(data,width,/*percentile*/0.50);
  }

  FREE(histogram);

  return;
}


/************************************************************************
 *   Linear fit:  Good for finding gradual drops
 ************************************************************************/


#define NDATA 30		/* ramp size */
#define NBASELINE 40
static double xcoef_up[NDATA];
static double xcoef_down[NDATA];
static double linear_xmean_up;
static double linear_xmean_down;
static double linear_den;

void
Tally_setup_xcoef () {
  double xi_up, xi_down;
  int i, j;

  linear_xmean_up = (double) (NDATA+1)/2;
  linear_xmean_down = -linear_xmean_up;

  linear_den = 0.0;
  for (i = 0, j = NDATA-1; i < NDATA; i++, j--) {
    xi_up = (double) (i+1);
    xi_down = -xi_up;
    xcoef_up[i] = xi_up - linear_xmean_up;
    xcoef_down[j] = xi_down - linear_xmean_down;
    linear_den += xcoef_up[i] * xcoef_up[i];
  }


  return;
}


static double
linear_xintercept (long int *counts_tally, long int *cum_tally, Genomicpos_T chrpos, double *xcoef,
		   double xmean, double baseline) {
  double xintercept,  yintercept, slope;
  long int sumy;
  double meany;
  int i;

  sumy = cum_tally[chrpos + NDATA - 1] - cum_tally[chrpos-1];
  if (sumy == 0) {
    return 0.0;			/* Terminates search for a low x-intercept */
  } else {
    meany = sumy/(double) NDATA;
  }

  slope = 0.0;
  for (i = 0; i < NDATA; i++, chrpos++) {
    slope += xcoef[i] * (double) counts_tally[chrpos];
  }
  slope /= linear_den;

  yintercept = meany - slope*xmean;

  if (slope == 0.0) {
    return 1000000.0;
  } else {
    xintercept = (baseline-yintercept)/slope;
    return xintercept;
  }
}



static double
linear_xintercept_double (double *log_tally, double *cumlog_tally, Genomicpos_T chrpos, double *xcoef,
			  double xmean, double baseline) {
  double xintercept,  yintercept, slope;
  double sumy;
  double meany;
  int i;

  sumy = log_tally[chrpos + NDATA - 1] - cumlog_tally[chrpos-1];
  if (sumy == 0.0) {
    return 0.0;			/* Terminates search for a low x-intercept */
  } else {
    meany = sumy/(double) NDATA;
  }

  slope = 0.0;
  for (i = 0; i < NDATA; i++, chrpos++) {
    debugT(printf("%f ",log_tally[chrpos]));
    slope += xcoef[i] * log_tally[chrpos];
  }
  debugT(printf("\n"));
  slope /= linear_den;

  yintercept = meany - slope*xmean;

  if (slope == 0.0) {
    return 1000000.0;
  } else {
    xintercept = (baseline-yintercept)/slope;
    return xintercept;
  }
}




/* Will not go to left of chrstart + NBASELINE + 10 */
/* For tolerance, lower numbers (like 5.0) extend exons further.  Higher numbers (like 10.0) truncate them */
Genomicpos_T
Tally_trace_leftward_linearfit (Genomicpos_T startpos, long int *counts_median, long int *cum_median,
				Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance) {
  Genomicpos_T chrpos;
  double xintercept;
#if 0
  int sumy;
  double baseline;
#endif

  chrpos = startpos;
  debugT(printf("trace_leftward at %u\n",startpos));
  while (chrpos >= chrstart + NBASELINE + 10U) {
#if 0
    sumy = cum_median[chrpos-1] - cum_median[chrpos - NBASELINE - 1];
    if (sumy == 0) {
      return 0.0;			/* Terminates search for a low x-intercept */
    } else {
      baseline = sumy/(double) NBASELINE;
    }
    xintercept = linear_xintercept(counts_median,cum_median,chrpos,xcoef_up,
				   linear_xmean_up,baseline);
#else
    xintercept = linear_xintercept(counts_median,cum_median,chrpos,xcoef_up,
				   linear_xmean_up,/*baseline*/0.0);
#endif
    debugT(printf("%u %.2f\n",chrpos,xintercept));
    if (fabs(xintercept) < tolerance) {
      debugT(printf("Done\n"));
      return chrpos;
    }

    chrpos--;
  }

  return chrpos;
}


/* For tolerance, lower numbers (like 5.0) extend exons further.  Higher numbers (like 10.0) truncate them */
Genomicpos_T
Tally_trace_leftward_loglinearfit (Genomicpos_T startpos, double *log_median, double *cumlog_median,
				   Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance) {
  Genomicpos_T chrpos;
  double xintercept;

  chrpos = startpos;
  debugT(printf("trace_leftward at %u\n",startpos));
  while (chrpos >= chrstart + NBASELINE + 10U) {
    xintercept = linear_xintercept_double(log_median,cumlog_median,chrpos,xcoef_up,
					  linear_xmean_up,/*baseline*/0.0);
    debugT(printf("%u %.2f\n",chrpos,xintercept));
    if (fabs(xintercept) < tolerance) {
      debugT(printf("Done\n"));
      return chrpos;
    }

    chrpos--;
  }

  return chrpos;
}




/* Will not go to right of chrend - NBASELINE - 10 */
Genomicpos_T
Tally_trace_rightward_linearfit (Genomicpos_T startpos, long int *counts_median, long int *cum_median,
				 Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance) {
  Genomicpos_T chrpos;
  double xintercept;
#if 0
  int sumy;
  double baseline;
#endif

  chrpos = startpos;
  debugT(printf("trace_rightward at %u\n",startpos));
  while (chrpos <= chrend - NBASELINE - 10U) {
#if 0
    sumy = cum_median[chrpos+NBASELINE] - cum_median[chrpos];
    if (sumy == 0) {
      return 0.0;			/* Terminates search for a low x-intercept */
    } else {
      baseline = sumy/(double) NBASELINE;
    }
    xintercept = linear_xintercept(counts_median,cum_median,chrpos-NDATA+1,xcoef_down,
				   linear_xmean_down,baseline);
#else
    xintercept = linear_xintercept(counts_median,cum_median,chrpos-NDATA+1,xcoef_down,
				   linear_xmean_down,/*baseline*/0.0);
#endif
    debugT(printf("%u %.2f\n",chrpos,xintercept));
    if (fabs(xintercept) < tolerance) {
      debugT(printf("Done\n"));
      return chrpos;
    }

    chrpos++;
  }

  return chrpos;
}


/* Will not go to right of chrend - NBASELINE - 10 */
Genomicpos_T
Tally_trace_rightward_loglinearfit (Genomicpos_T startpos, double *log_median, double *cumlog_median,
				    Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance) {
  Genomicpos_T chrpos;
  double xintercept;

  chrpos = startpos;
  debugT(printf("trace_rightward at %u\n",startpos));
  while (chrpos <= chrend - NBASELINE - 10U) {
    xintercept = linear_xintercept_double(log_median,cumlog_median,chrpos-NDATA+1,xcoef_down,
					  linear_xmean_down,/*baseline*/0.0);
    debugT(printf("%u %.2f\n",chrpos,xintercept));
    if (fabs(xintercept) < tolerance) {
      debugT(printf("Done\n"));
      return chrpos;
    }

    chrpos++;
  }

  return chrpos;
}




/************************************************************************/

#define TESTREGION 100

#if 0
static double
loglik (int *data, int n, double lambda) {
  double loglik = 0.0;
  int sumx = 0;
  int i;

  if (lambda < 0.1) {
    /* Prevent taking log of zero */
    lambda = 0.1;
  }

  for (i = 0; i < n; i++) {
    /* printf("Testing data %d against lambda %f\n",data[i],lambda); */
    sumx += data[i];
  }

  loglik = log(lambda)*(double) sumx - (double) n*lambda;
  for (i = 0; i < n; i++) {
    loglik -= gammln(data[i]+1.0);
  }
  /* printf("loglik = %g\n",loglik); */

  return loglik;
}
#endif



Genomicpos_T
Tally_trace_leftward_loglik (Genomicpos_T downpos, long int *counts_tally, long int *cum_tally,
			     Genomicpos_T chrstart, Genomicpos_T chrend, 
			     int min_exonlength, int max_exonlength) {
  Genomicpos_T changepoint, bestchangepoint = downpos, endpos;
  double best_ll = 0.0, ll;
  long int exon_sumx, intron_sumx;
  double exon_mean, intron_mean;

  if (downpos < max_exonlength) {
    endpos = 0U;
  } else {
    endpos = downpos - max_exonlength;
  }
  if (endpos < chrstart + TESTREGION) {
    endpos = chrstart + TESTREGION;
  }

  for (changepoint = downpos - min_exonlength; changepoint > endpos; changepoint--) {
    intron_sumx = cum_tally[changepoint] - cum_tally[changepoint-TESTREGION];
    exon_sumx = cum_tally[downpos] - cum_tally[changepoint];

    intron_mean = (double) intron_sumx/(double) TESTREGION;
    exon_mean = (double) exon_sumx/(double) (downpos - changepoint);

    debugL(printf("testing up at %u, exon_mean:%.1f intron_mean:%.1f",
		  changepoint,exon_mean,intron_mean));
    if (exon_mean > intron_mean) {
      ll = Lgamma_dpois_log_n(&(counts_tally[changepoint]),downpos-changepoint,intron_mean) +
	Lgamma_dpois_log_n(&(counts_tally[changepoint-TESTREGION]),TESTREGION,exon_mean);
      ll /= (double) (downpos - changepoint + TESTREGION);
      debugL(printf(" loglik:%.1f",ll));
      if (ll < best_ll) {
	debugL(printf(" **"));
	best_ll = ll;
	bestchangepoint = changepoint;
      }
    }
    debugL(printf("\n"));
  }

  return bestchangepoint;
}


Genomicpos_T
Tally_trace_rightward_loglik (Genomicpos_T uppos, long int *counts_tally, long int *cum_tally,
			      Genomicpos_T chrstart, Genomicpos_T chrend,
			      int min_exonlength, int max_exonlength) {
  Genomicpos_T changepoint, bestchangepoint = uppos, endpos;
  double best_ll = 0.0, ll;
  long int exon_sumx, intron_sumx;
  double exon_mean, intron_mean;

  if ((endpos = uppos + max_exonlength) > chrend - TESTREGION) {
    endpos = chrend - TESTREGION;
  }

  for (changepoint = uppos + min_exonlength; changepoint < endpos; changepoint++) {
    exon_sumx = cum_tally[changepoint] - cum_tally[uppos];
    intron_sumx = cum_tally[changepoint+TESTREGION] - cum_tally[changepoint];
	      
    exon_mean = (double) exon_sumx/(double) (changepoint - uppos);
    intron_mean = (double) intron_sumx/(double) TESTREGION;

    debugL(printf("testing down at %u, exon_mean: %.1f intron_mean:%.1f",
		  changepoint,exon_mean,intron_mean));
    if (exon_mean > intron_mean) {
      ll = Lgamma_dpois_log_n(&(counts_tally[uppos]),changepoint-uppos,intron_mean) + 
	Lgamma_dpois_log_n(&(counts_tally[changepoint]),TESTREGION,exon_mean);
      ll /= (double) (changepoint - uppos + TESTREGION);
      debugL(printf(" loglik:%.1f",ll));
      if (ll < best_ll) {
	debugL(printf(" **"));
	best_ll = ll;
	bestchangepoint = changepoint;
      }
    }
    debugL(printf("\n"));
  }

  return bestchangepoint;
}



/************************************************************************/


#if 0
static double *prior;

void
Tally_prior_init (int tau, int max_exonlength) {
  int i;

  prior = (double *) CALLOC(max_exonlength+1,sizeof(double));

  for (i = 0; i <= max_exonlength; i++) {
    prior[i] = exp(-(double) i/(double) tau);
  }
  return;
}



static int
region_max (int *counts, int n) {
  int maxcount = 0;
  int i;

  for (i = 0; i < n; i++) {
    if (counts[i] > maxcount) {
      maxcount = counts[i];
    }
  }
  return maxcount;
}


static void
region_stats (int *sumx, int *sumxx, int *counts, int n) {
  int count;
  int i;

  *sumx = *sumxx = 0;
  for (i = 0; i < n; i++) {
    count = counts[i];
    *sumx += count;
    *sumxx += count*count;
  }
  return;
}


Genomicpos_T
Tally_trace_leftward_varmean (Genomicpos_T downpos, int *counts_tally,
			      Genomicpos_T chrstart, Genomicpos_T chrend, 
			      int auto_exonlength, int max_exonlength, bool debugp) {
  Genomicpos_T bestchangepoint = downpos, changepoint, endpos;
  int maxcount, oldcount, newcount, sumx, sumxx;
  double maxvalue = 0.0, value, variance, meanx;

  if (downpos < max_exonlength) {
    endpos = 0U;
  } else {
    endpos = downpos - max_exonlength;
  }
  if (endpos < chrstart + TESTREGION) {
    endpos = chrstart + TESTREGION;
  }


  changepoint = downpos - auto_exonlength + 1;
  region_stats(&sumx,&sumxx,&(counts_tally[changepoint - TESTREGION/2]),TESTREGION);

  for (changepoint--; changepoint > endpos; changepoint--) {
    /* Compute variance centered at changepoint */
    oldcount = counts_tally[changepoint - TESTREGION/2 + TESTREGION];
    sumx -= oldcount;
    sumxx -= oldcount*oldcount;
    
    newcount = counts_tally[changepoint - TESTREGION/2];
    sumx += newcount;
    sumxx += newcount*newcount;

    meanx = (double) sumx / (double) TESTREGION;
    variance = ((double) sumxx - 2.0*meanx * (double) sumx + ((double) TESTREGION)*meanx*meanx) / (double) (TESTREGION-1);

    /* Compute max past changepoint */
    maxcount = region_max(&(counts_tally[changepoint - TESTREGION]),TESTREGION);

    /* Compute value */
    if ((value = prior[downpos - changepoint] * variance / (double) (maxcount+1)) > maxvalue) {
      bestchangepoint = changepoint;
      maxvalue = value;
    }

    if (debugp) {
      printf("%u %d %f %f %f %d\n",
	     changepoint,counts_tally[changepoint],value,prior[downpos - changepoint],variance,maxcount);
    }

  }

  return bestchangepoint;
}


Genomicpos_T
Tally_trace_rightward_varmean (Genomicpos_T uppos, int *counts_tally,
			       Genomicpos_T chrstart, Genomicpos_T chrend,
			       int auto_exonlength, int max_exonlength, bool debugp) {
  Genomicpos_T bestchangepoint = uppos, changepoint, endpos;
  int maxcount, oldcount, newcount, sumx, sumxx;
  double maxvalue = 0.0, value, variance, meanx;

  if ((endpos = uppos + max_exonlength) > chrend - TESTREGION) {
    endpos = chrend - TESTREGION;
  }

  changepoint = uppos + auto_exonlength - 1;
  region_stats(&sumx,&sumxx,&(counts_tally[changepoint - TESTREGION/2]),TESTREGION);

  for (changepoint++; changepoint < endpos; changepoint++) {
    /* Compute variance centered at changepoint */
    oldcount = counts_tally[changepoint - TESTREGION/2 - 1];
    sumx -= oldcount;
    sumxx -= oldcount*oldcount;

    newcount = counts_tally[changepoint - TESTREGION/2 - 1 + TESTREGION];
    sumx += newcount;
    sumxx += newcount*newcount;

    meanx = (double) sumx / (double) TESTREGION;
    variance = ((double) sumxx - 2.0*meanx * (double) sumx + ((double) TESTREGION)*meanx*meanx) / (double) (TESTREGION-1);

    /* Compute max past changepoint */
    maxcount = region_max(&(counts_tally[changepoint]),TESTREGION);

    /* Compute value */
    if ((value = prior[changepoint - uppos] * variance / (double) (maxcount+1)) > maxvalue) {
      bestchangepoint = changepoint;
      maxvalue = value;
    }

    if (debugp) {
      printf("%u %d %f %f %f %d\n",
	     changepoint,counts_tally[changepoint],value,prior[changepoint - uppos],variance,maxcount);
    }

  }

  return bestchangepoint;
}


#endif
