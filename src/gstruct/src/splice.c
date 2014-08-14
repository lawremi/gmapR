static char rcsid[] = "$Id: splice.c 136513 2014-05-16 17:58:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splice.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For qsort */
#include <ctype.h>		/* For isspace */
#include <string.h>
#include "mem.h"
#include "assert.h"
#include "interval.h"


char *
Splicestatus_string (Splicestatus_T status) {
  switch (status) {
  case VALID: return "";
  case MINCOUNT: return "mincount";
  case MINSUPPORT: return "minsupport";
  case BADPROB: return "badprob";
  case NONCANONICAL: return "noncanonical";
  case WRONGDIR_CMP: return "wrongdir-cmp";
  case WRONGDIR_EXTENT: return "wrongdir-extent";
  }
  return (char *) NULL;
}


#define T Splice_T
struct T {
  bool intron_usedp;
  int known_count;
  int primary_count;
  int crosshyb_count;
  int nunique;

  int maxminsupport;
  int max_overhang1_lowend;
  int max_overhang2_lowend;
  int max_overhang1_highend;
  int max_overhang2_highend;
  int max_querylength_lowend;
  int max_querylength_highend;

  int sign;
  Genomicpos_T donorpos;
  Genomicpos_T acceptorpos;
  Genomicpos_T low;
  Genomicpos_T high;

  bool canonicalp;
  double donorprob;
  double acceptorprob;

  char donor1;
  char donor2;
  char acceptor1;
  char acceptor2;

  int alpha;
  int omega;

  Splicestatus_T status;
  Spliceknownp_T knownp;
  bool valid_end_p;
  int level;
};


void
Splice_free (T *old) {
  FREE(*old);
  return;
}

void
Splice_gc (List_T *splices) {
  List_T p;
  T splice;

  for (p = *splices; p != NULL; p = List_next(p)) {
    splice = (T) List_head(p);
    Splice_free(&splice);
  }
  List_free(&(*splices));
  return;
}

bool
Splice_intron_usedp (T this) {
  return this->intron_usedp;
}

void
Splice_set_intron_usedp (T this) {
  this->intron_usedp = true;
}


int
Splice_count (T this) {
  return this->primary_count + this->crosshyb_count;
}

int
Splice_known_count (T this) {
  return this->known_count;
}

int
Splice_primary_count (T this) {
  return this->primary_count;
}

int
Splice_crosshyb_count (T this) {
  return this->crosshyb_count;
}

int
Splice_maxminsupport (T this) {
  return this->maxminsupport;
}

#if 0
int
Splice_npositions_obs (T this) {
  int npositions, npositions_lowend, npositions_highend;

  if ((npositions_lowend = this->max_overhang1_lowend + this->max_overhang2_lowend - this->max_querylength_lowend) < 0) {
    npositions_lowend = 0;
  }
  if ((npositions_highend = this->max_overhang1_highend + this->max_overhang2_highend - this->max_querylength_highend) < 0) {
    npositions_highend = 0;
  }
  
  if ((npositions = npositions_lowend + npositions_highend) == 0) {
    return 1;
  } else {
    return npositions;
  }
}
#endif

int
Splice_npositions (T this, Genomicpos_T genestart, Genomicpos_T geneend, int insertlength, int readlength,
		   int min_overhang) {
  Genomicpos_T genelow, genehigh, pos;
  int npositions, npositions_lowend, npositions_highend;
  int max_overhang1_lowend, max_overhang2_lowend, max_overhang1_highend, max_overhang2_highend;
  int adj_readlength;

  if (genestart < geneend) {
    genelow = genestart;
    genehigh = geneend;
  } else {
    genelow = geneend;
    genehigh = genestart;
  }

  adj_readlength = readlength - min_overhang;

  if (genelow > this->low - adj_readlength) {
    max_overhang1_lowend = this->low - genelow;
  } else {
    max_overhang1_lowend = adj_readlength;
  }
#ifdef USE_OBSERVED
  if (this->max_overhang1_lowend > max_overhang1_lowend) {
    max_overhang1_lowend = this->max_overhang1_lowend;
  }
#endif


  if ((pos = genelow + insertlength - adj_readlength) > this->low) {
    max_overhang1_highend = 0;
  } else if (pos > this->low - adj_readlength) {
    max_overhang1_highend = this->low - pos;
  } else {
    max_overhang1_highend = adj_readlength;
  }
#ifdef USE_OBSERVED
  if (this->max_overhang1_highend > max_overhang1_highend) {
    max_overhang1_highend = this->max_overhang1_highend;
  }
#endif


  if (genehigh < this->high + adj_readlength) {
    max_overhang2_highend = genehigh - this->high;
  } else {
    max_overhang2_highend = adj_readlength;
  }
#ifdef USE_OBSERVED
  if (this->max_overhang2_highend > max_overhang2_highend) {
    max_overhang2_highend = this->max_overhang2_highend;
  }
#endif


  if ((pos = genehigh - insertlength + adj_readlength) < this->high) {
    max_overhang2_lowend = 0;
  } else if (pos < this->high + adj_readlength) {
    max_overhang2_lowend = pos - this->high;
  } else {
    max_overhang2_lowend = adj_readlength;
  }
#ifdef USE_OBSERVED
  if (this->max_overhang2_lowend > max_overhang2_lowend) {
    max_overhang2_lowend = this->max_overhang2_lowend;
  }
#endif
  

  if ((npositions_lowend = max_overhang1_lowend + max_overhang2_lowend - readlength) < 0) {
    npositions_lowend = 0;
  }
  if ((npositions_highend = max_overhang1_highend + max_overhang2_highend - readlength) < 0) {
    npositions_highend = 0;
  }

#if 0
  if ((npositions = npositions_lowend + npositions_highend) == 0) {
    return 1;
  } else {
    return npositions;
  }
#else
  return npositions_lowend + npositions_highend;
#endif
}



double
Splice_density (T this, Genomicpos_T genestart, Genomicpos_T geneend, int insertlength, int readlength,
		int min_overhang) {
  int npositions = Splice_npositions(this,genestart,geneend,insertlength,readlength,min_overhang);

  return (double) (this->primary_count + this->crosshyb_count)/(double) (npositions + 1);
}



int
Splice_sign (T this) {
  return this->sign;
}

Genomicpos_T
Splice_donorpos (T this) {
  return this->donorpos;
}

Genomicpos_T
Splice_acceptorpos (T this) {
  return this->acceptorpos;
}

Genomicpos_T
Splice_low (T this) {
  return this->low;
}

Genomicpos_T
Splice_high (T this) {
  return this->high;
}

bool
Splice_canonicalp (T this) {
  return this->canonicalp;
}

double
Splice_donorprob (T this) {
  return this->donorprob;
}

double
Splice_acceptorprob (T this) {
  return this->acceptorprob;
}

Splicestatus_T
Splice_status (T this) {
  return this->status;

}

Spliceknownp_T
Splice_knownp (T this) {
  return this->knownp;
}


bool
Splice_validp (T this) {
  if (this->status == VALID) {
    return true;
  } else {
    return false;
  }
}

bool
Splice_valid_end_p (T this) {
  return this->valid_end_p;
}

int
Splice_level (T this) {
  return this->level;
}

T *
Splice_array_copy (T *array, int nsplices) {
  T *copy;
  int i;

  if (nsplices == 0) {
    return (T *) NULL;
  } else {
    copy = (T *) CALLOC(nsplices,sizeof(T));
    for (i = 0; i < nsplices; i++) {
      copy[i] = array[i];
    }
    return copy;
  }
}

List_T
Splice_valid_list (List_T splices) {
  List_T valid = NULL, p;
  T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (splice->status == VALID) {
      valid = List_push(valid,(void *) splice);
    }
  }

  return valid;
}



int
Splice_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->donorpos < y->donorpos) {
    return -1;
  } else if (y->donorpos < x->donorpos) {
    return +1;
  } else if (x->acceptorpos < y->acceptorpos) {
    return -1;
  } else if (y->acceptorpos < x->acceptorpos) {
    return +1;
  } else {
    return 0;
  }
}

int
Splice_low_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1; 
  } else {
    return 0;
  }
}

int
Splice_high_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->high > y->high) {
    return -1;
  } else if (y->high > x->high) {
    return +1;
  } else {
    return 0;
  }
}

int
Splice_bycount_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->primary_count + x->crosshyb_count < y->primary_count + y->crosshyb_count) {
    return -1;
  } else if (y->primary_count + y->crosshyb_count < x->primary_count + x->crosshyb_count) {
    return +1;
  } else {
    return 0;
  }
}


static int
level_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->low < b->low) {
    return -1;
  } else if (a->low > b->low) {
    return +1;
  } else if (a->high < b->high) {
    return -1;
  } else if (a->high > b->high) {
    return +1;
  } else if (a->sign > b->sign) {
    return -1;
  } else if (a->sign < b->sign) {
    return +1;
  } else {
    return 0;
  }
}



int
Splice_compute_levels (List_T splices, Genomicpos_T mincoord,
		       Genomicpos_T maxcoord, int max_allowed_levels,
		       double xfactor, int mincount) {
  int nsplices, i;
  int maxlevel = -1, level;
  bool donep;
  T *array, splice;
  double *rightmost, xlow;

  if ((nsplices = List_length(splices)) > 0) {
    array = (T *) List_to_array(splices,NULL);
    qsort(array,nsplices,sizeof(T),level_cmp);

    rightmost = (double *) CALLOC(max_allowed_levels,sizeof(double));
    for (i = 0; i < max_allowed_levels; i++) {
      rightmost[i] = 0.0;
    }

    for (i = 0; i < nsplices; i++) {
      splice = array[i];
      if (splice->primary_count + splice->crosshyb_count < mincount) {
	/* Skip */
	splice->level = -1;
      } else {
	/* Find appropriate level */
	level = 0;
	donep = false;
	while (level < max_allowed_levels && !donep) {
	  xlow = xfactor * splice->low;
	  if (level > maxlevel) {
	    donep = true;
	    maxlevel = level;
	  } else if (rightmost[level] < xlow) {
	    donep = true;
	  } else {
	    level++;
	  }
	}

	if (level < max_allowed_levels) {
	  rightmost[level] = xfactor * splice->high + 3.0;

	  if (splice->high < mincoord || splice->low > maxcoord) {
	    /* Skip printing if both ends outside of region */
	  } else if (splice->low < mincoord || splice->high > maxcoord) {
	    /* Skip printing if either end outside of region */
	  } else {
	    splice->level = level;
	  }
	}
      }
    }

    FREE(rightmost);
    FREE(array);
  }

  return maxlevel + 1;
}


T
Splice_new_known (Genomicpos_T donorpos, Genomicpos_T acceptorpos, int sign,
		  bool mainpath_p) {
  T new = (T) MALLOC(sizeof(*new));

  new->intron_usedp = false;
  new->known_count = 1;
  new->primary_count = 0;	/* Do not add to observed counts */
  new->crosshyb_count = 0;	/* Do not add to observed counts */
  new->nunique = 0;		/* Do not add to observed counts */
  new->maxminsupport = 100;

  new->sign = sign;
  new->donorpos = donorpos;
  new->acceptorpos = acceptorpos;
  if (sign > 0) {
    new->low = donorpos;
    new->high = acceptorpos;
  } else if (sign < 0) {
    new->low = acceptorpos;
    new->high = donorpos;
  } else {
    fprintf(stderr,"sign is 0\n");
    abort();
  }

  new->canonicalp = true;
  new->donorprob = 2.0;
  new->acceptorprob = 2.0;
  new->donor1 = ' ';
  new->donor2 = ' ';
  new->acceptor1 = ' ';
  new->acceptor2 = ' ';
  new->status = VALID;
  if (mainpath_p == true) {
    new->knownp = KNOWN_MAIN;
  } else {
    new->knownp = KNOWN_ALT;
  }
  new->level = -1;

  new->valid_end_p = true;

  return new;
}


void
Splice_transfer_info (T knownsplice, T obssplice) {
  knownsplice->primary_count = obssplice->primary_count;
  knownsplice->crosshyb_count = obssplice->crosshyb_count;
  knownsplice->nunique = obssplice->nunique;

  knownsplice->max_overhang1_lowend = obssplice->max_overhang1_lowend;
  knownsplice->max_overhang2_lowend = obssplice->max_overhang2_lowend;
  knownsplice->max_overhang1_highend = obssplice->max_overhang1_highend;
  knownsplice->max_overhang2_highend = obssplice->max_overhang2_highend;
  knownsplice->max_querylength_lowend = obssplice->max_querylength_lowend;
  knownsplice->max_querylength_highend = obssplice->max_querylength_highend;

  return;
}


List_T
Splice_add_at_donor (List_T list, int sign, Genomicpos_T donorpos, Genomicpos_T acceptorpos,
		     bool crosshybp, int nhits, int support1, int support2,
		     int overhang1, int overhang2, int querylength, bool lowend_p, bool canonicalp,
		     double donorprob, double acceptorprob,
		     char donor1, char donor2, char acceptor1, char acceptor2, bool knownp) {
  List_T p;
  T splice;
  int minsupport;

  if (support1 < support2) {
    minsupport = support1;
  } else {
    minsupport = support2;
  }

  for (p = list; p != NULL; p = List_next(p)) {
    splice = List_head(p);
    if (splice->acceptorpos == acceptorpos) {
      if (knownp == true) {
	splice->known_count += 1;
      } else if (crosshybp == true) {
	splice->crosshyb_count += 1;
      } else {
	splice->primary_count += 1;
      }
      if (nhits == 1) {
	splice->nunique += 1;
      }
      if (minsupport > splice->maxminsupport) {
	splice->maxminsupport = minsupport;
      }
      if (lowend_p == true) {
	if (overhang1 > splice->max_overhang1_lowend) {
	  splice->max_overhang1_lowend = overhang1;
	}
	if (overhang2 > splice->max_overhang2_lowend) {
	  splice->max_overhang2_lowend = overhang2;
	}
	if (querylength > splice->max_querylength_lowend) {
	  splice->max_querylength_lowend = querylength;
	}

      } else {
	if (overhang1 > splice->max_overhang1_highend) {
	  splice->max_overhang1_highend = overhang1;
	}
	if (overhang2 > splice->max_overhang2_lowend) {
	  splice->max_overhang2_highend = overhang2;
	}
	if (querylength > splice->max_querylength_highend) {
	  splice->max_querylength_highend = querylength;
	}
      }

      return list;
    }
  }

  /* Not found, so add to list */
  splice = (T) MALLOC(sizeof(*splice));
  splice->intron_usedp = false;
  if (knownp == true) {
    splice->known_count = 1;
    splice->primary_count = 0;
    splice->crosshyb_count = 0;
  } else if (crosshybp == true) {
    splice->known_count = 0;
    splice->primary_count = 0;
    splice->crosshyb_count = 1;
  } else {
    splice->known_count = 0;
    splice->primary_count = 1;
    splice->crosshyb_count = 0;
  }
  splice->nunique = 0;
  if (nhits == 1) {
    splice->nunique += 1;
  }
  splice->maxminsupport = minsupport;
  if (lowend_p == true) {
    splice->max_overhang1_lowend = overhang1;
    splice->max_overhang2_lowend = overhang2;
    splice->max_overhang1_highend = splice->max_overhang2_highend = 0;
    splice->max_querylength_lowend = querylength;
    splice->max_querylength_highend = 0;
  } else {
    splice->max_overhang1_lowend = splice->max_overhang2_lowend = 0;
    splice->max_overhang1_highend = overhang1;
    splice->max_overhang2_highend = overhang2;
    splice->max_querylength_lowend = 0;
    splice->max_querylength_highend = querylength;
  }

  splice->sign = sign;
  splice->donorpos = donorpos;
  splice->acceptorpos = acceptorpos;
  if (sign > 0) {
    splice->low = donorpos;
    splice->high = acceptorpos;
  } else if (sign < 0) {
    splice->low = acceptorpos;
    splice->high = donorpos;
  } else {
    fprintf(stderr,"sign is 0\n");
    abort();
  }

  splice->canonicalp = canonicalp;
  splice->donorprob = donorprob;
  splice->acceptorprob = acceptorprob;
  splice->donor1 = donor1;
  splice->donor2 = donor2;
  splice->acceptor1 = acceptor1;
  splice->acceptor2 = acceptor2;
  splice->status = VALID;
  splice->knownp = NOVEL;
  splice->level = -1;

  splice->valid_end_p = true;

  return List_push(list,(void *) splice);
}


#if 0
static int
strongest (T fwd_splice, T rev_splice, 
	   long int *fwd_runlengths, long int *rev_runlengths,
	   T *array, int a, int b, int nsplices, char *chr,
	   Genomicpos_T chrlength, Genomicpos_T max_exonlength) {
  long int fwd_max, rev_max;
  int i;

  fwd_max = Tally_maxcount(fwd_runlengths,chrlength,fwd_splice->low,fwd_splice->high);
  rev_max = Tally_maxcount(rev_runlengths,chrlength,rev_splice->low,rev_splice->high);
  if (fwd_max + 1 > 5*(rev_max + 1)) {
    fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld)\n",
	    chr,fwd_splice->low,fwd_splice->high,fwd_max,
	    chr,rev_splice->high,rev_splice->low,rev_max);
    return +1;
  } else if (rev_max + 1 > 5*(fwd_max + 1)) {
    fprintf(stderr,"Choose -%s:%u..%u (%ld) over +%s:%u..%u (%ld)\n",
	    chr,rev_splice->high,rev_splice->low,rev_max,
	    chr,fwd_splice->low,fwd_splice->high,fwd_max);
    return -1;
  } else {
    i = a-1;
    while (i >= 0 && array[i]->high + max_exonlength > fwd_splice->low && array[i]->high > rev_splice->low) {
      fwd_max = Tally_maxcount(fwd_runlengths,chrlength,array[i]->low,array[i]->high);
      rev_max = Tally_maxcount(rev_runlengths,chrlength,array[i]->low,array[i]->high);
      if (fwd_max + 1 > 5*(rev_max + 1)) {
	fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld), secondary left\n",
		chr,fwd_splice->low,fwd_splice->high,fwd_max,
		chr,rev_splice->high,rev_splice->low,rev_max);
	return +1;
      } else if (rev_max + 1 > 5*(fwd_max + 1)) {
	fprintf(stderr,"Choose -%s:%u..%u (%ld) over +%s:%u..%u (%ld), secondary left\n",
		chr,rev_splice->high,rev_splice->low,rev_max,
		chr,fwd_splice->low,fwd_splice->high,fwd_max);
	return -1;
      }
      i--;
    }

    i = b+1;
    while (i < nsplices && array[i]->low < fwd_splice->high + max_exonlength && array[i]->low < rev_splice->high + max_exonlength) {
      fwd_max = Tally_maxcount(fwd_runlengths,chrlength,array[i]->low,array[i]->high);
      rev_max = Tally_maxcount(rev_runlengths,chrlength,array[i]->low,array[i]->high);
      if (fwd_max + 1 > 5*(rev_max + 1)) {
	fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld), secondary right\n",
		chr,fwd_splice->low,fwd_splice->high,fwd_max,
		chr,rev_splice->high,rev_splice->low,rev_max);
	return +1;
      } else if (rev_max + 1 > 5*(fwd_max + 1)) {
	fprintf(stderr,"Choose -%s:%u..%u (%ld) over +%s:%u..%u (%ld), secondary right\n",
		chr,rev_splice->high,rev_splice->low,rev_max,
		chr,fwd_splice->low,fwd_splice->high,fwd_max);
	return -1;
      }
      i++;
    }
  }

  fprintf(stderr,"Cannot choose between +%s:%u..%u (%ld) and -%s:%u..%u (%ld)\n",
	  chr,fwd_splice->low,fwd_splice->high,fwd_max,chr,rev_splice->high,rev_splice->low,rev_max);

  return 0;
}
#endif


static long int
compute_maxcount (long int *x, Genomicpos_T chrlength, Genomicpos_T chrstart, Genomicpos_T chrend) {
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



static bool
consistentp (T splice, long int *fwd_runlengths, long int *rev_runlengths,
	     char *chr, Genomicpos_T chrlength) {
  long int fwd_max, rev_max;

  fwd_max = compute_maxcount(fwd_runlengths,chrlength,splice->low,splice->high);
  rev_max = compute_maxcount(rev_runlengths,chrlength,splice->low,splice->high);
  if (splice->sign > 0) {
    if (rev_max + 1 > 5*(fwd_max + 1)) {
      fprintf(stderr,"Inconsistent: +%s:%u..%u (%ld fwd, %ld rev)\n",
	      chr,splice->low,splice->high,fwd_max,rev_max);
      return false;
    }
  } else if (splice->sign < 0) {
    if (fwd_max + 1 > 5*(rev_max + 1)) {
      fprintf(stderr,"Inconsistent: -%s:%u..%u (%ld fwd, %ld rev)\n",
	      chr,splice->high,splice->low,fwd_max,rev_max);
      return false;
    }
  }

  return true;
}



static int
sufficient_splice_prob_local (int support, double spliceprob) {
  if (support < 14) {
    return (spliceprob > 0.95);
  } else if (support < 20) {
    return (spliceprob > 0.90);
  } else if (support < 25) {
    return (spliceprob > 0.85);
  } else if (support < 30) {
    return (spliceprob > 0.70);
  } else {
    return 1;
  }
}


void
Splice_filter_list (List_T splices, int mincount_end_alt, int minsupport, bool need_canonical_p) {
  List_T p;
  T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (T) List_head(p);
    if (need_canonical_p == true && splice->canonicalp == false) {
      splice->status = NONCANONICAL;
    }

    if (sufficient_splice_prob_local(splice->maxminsupport,splice->donorprob) == 0) {
      splice->valid_end_p = false;
    } else if (sufficient_splice_prob_local(splice->maxminsupport,splice->acceptorprob) == 0) {
      splice->valid_end_p = false;
    } else if (splice->primary_count + splice->crosshyb_count < mincount_end_alt) {
      splice->valid_end_p = false;
    } else if (splice->maxminsupport < minsupport) {
      splice->valid_end_p = false;
    } else {
      splice->valid_end_p = true;
    }
  }

  return;
}


void
Splice_filter_array (T *array, int nsplices, int mincount_end_alt, int minsupport, bool need_canonical_p) {
  T splice;
  int a;

  for (a = 0; a < nsplices; a++) {
    splice = array[a];
    if (need_canonical_p == true && splice->canonicalp == false) {
      splice->status = NONCANONICAL;
    }

    if (sufficient_splice_prob_local(splice->maxminsupport,splice->donorprob) == 0) {
      splice->valid_end_p = false;
    } else if (sufficient_splice_prob_local(splice->maxminsupport,splice->acceptorprob) == 0) {
      splice->valid_end_p = false;
    } else if (splice->primary_count + splice->crosshyb_count < mincount_end_alt) {
      splice->valid_end_p = false;
    } else if (splice->maxminsupport < minsupport) {
      splice->valid_end_p = false;
    } else {
      splice->valid_end_p = true;
    }
  }

  return;
}


void
Splice_resolve (T *array, int nsplices, long int *fwd_runlengths, long int *rev_runlengths,
		Genomicpos_T chrlength, Genomicpos_T max_exonlength, char *chr) {
  T splice1;
  int a;

#if 0
  qsort(array,nsplices,sizeof(T),Splice_low_cmp);
  for (a = 0; a < nsplices; a++) {
    splice1 = array[a];
    for (b = a+1; b < nsplices; b++) {
      splice2 = array[b];
      if (splice2->low > splice1->low + 20) {
	b = nsplices;	/* End loop */

      } else if (splice1->sign == splice2->sign) {
	/* Don't resolve between splices on the same strand */

      } else if (splice1->sign > 0) {
	if ((result = strongest(splice1,splice2,fwd_runlengths,rev_runlengths,array,a,b,nsplices,chr,
				chrlength,max_exonlength)) > 0) {
	  splice2->status = WRONGDIR_CMP;
	} else if (result < 0) {
	  splice1->status = WRONGDIR_CMP;
	}
	      
      } else {
	if ((result = strongest(splice2,splice1,fwd_runlengths,rev_runlengths,array,a,b,nsplices,chr,
				chrlength,max_exonlength)) > 0) {
	  splice1->status = WRONGDIR_CMP;
	} else if (result < 0) {
	  splice2->status = WRONGDIR_CMP;
	}
      }
    }
  }

  qsort(array,nsplices,sizeof(T),Splice_high_cmp);
  for (a = 0; a < nsplices; a++) {
    splice1 = array[a];
    for (b = a+1; b < nsplices; b++) {
      splice2 = array[b];
      if (splice2->high < splice1->high - 20) {
	b = nsplices;	/* End loop */

      } else if (splice1->sign == splice2->sign) {
	/* Don't resolve between splices on the same strand */

      } else if (splice1->sign > 0) {
	if ((result = strongest(splice1,splice2,fwd_runlengths,rev_runlengths,array,a,b,nsplices,chr,
				chrlength,max_exonlength)) > 0) {
	  splice2->status = WRONGDIR_CMP;
	} else if (result < 0) {
	  splice1->status = WRONGDIR_CMP;
	}

      } else {
	if ((result = strongest(splice2,splice1,fwd_runlengths,rev_runlengths,array,a,b,nsplices,chr,
				chrlength,max_exonlength)) > 0) {
	  splice1->status = WRONGDIR_CMP;
	} else if (result < 0) {
	  splice2->status = WRONGDIR_CMP;
	}
      }
    }
  }
#endif

  for (a = 0; a < nsplices; a++) {
    splice1 = array[a];
    if (splice1->status == VALID) {
      if (consistentp(splice1,fwd_runlengths,rev_runlengths,chr,chrlength) == false) {
	splice1->status = WRONGDIR_EXTENT;
      }
    }
  }

  return;
}



#if 0
/* Not used anymore */
void
Splice_turn (T *array, int nsplices, Genomicpos_T max_exonlength, char *chr) {
  int i, j, k;
  T splice, left_splice, right_splice;
  T *alpha_array, *omega_array;

  alpha_array = Splice_array_copy(array,nsplices);
  qsort(alpha_array,nsplices,sizeof(T),Splice_low_cmp);
  for (i = 0; i < nsplices; i++) {
    alpha_array[i]->alpha = i;
  }

  omega_array = Splice_array_copy(array,nsplices);
  qsort(omega_array,nsplices,sizeof(T),Splice_high_cmp);
  for (i = 0; i < nsplices; i++) {
    omega_array[i]->omega = i;
  }

  qsort(array,nsplices,sizeof(T),Splice_bycount_cmp);

  for (i = 0; i < nsplices; i++) {
    splice = array[i];
    if (splice->status == VALID) {
      /* fprintf(stderr,"Splice %c%s:%u..%u (%d)\n",
	 splice->sign > 0 ? '+' : '-',chr,splice->donorpos,splice->acceptorpos,splice->count); */

      /* Get left splices */
      j = splice->omega;
      while (j < nsplices && omega_array[j]->high >= splice->low) {
	/* fprintf(stderr,"  advancing: omega %d, high %u\n",j,omega_array[j]->high); */
	j++;
      }
      for ( ; j < nsplices && omega_array[j]->high + max_exonlength >= splice->low && splice->status == VALID; j++) {
	/* fprintf(stderr,"  evaluating: omega %d, high %u\n",j,omega_array[j]->high); */
	left_splice = omega_array[j];
	if (left_splice->status == VALID) {

	  /* Get right splices */
	  k = splice->alpha + 1;
	  while (k < nsplices && alpha_array[k]->low <= splice->high) {
	    k++;
	  }
	  for ( ; k < nsplices && alpha_array[k]->low <= splice->high + max_exonlength && splice->status == VALID; k++) {
	    right_splice = alpha_array[k];
	    if (right_splice->status == VALID) {

#if 0
	      fprintf(stderr,"Splice %c%s:%u..%u (%d) compared with %c%s:%u..%u (%d) and  %c%s:%u..%u (%d)\n",
		      splice->sign > 0 ? '+' : '-',chr,splice->donorpos,splice->acceptorpos,splice->count,
		      left_splice->sign > 0 ? '+' : '-',chr,left_splice->donorpos,left_splice->acceptorpos,left_splice->count,
		      right_splice->sign > 0 ? '+' : '-',chr,right_splice->donorpos,right_splice->acceptorpos,right_splice->count);
#endif


	      if (left_splice->sign != right_splice->sign) {
#if 0
		fprintf(stderr,"Splices %c%s:%u..%u (%d) and %c%s:%u..%u (%d) are inconsistent\n",
			left_splice->sign > 0 ? '+' : '-',chr,left_splice->donorpos,left_splice->acceptorpos,left_splice->count,
			right_splice->sign > 0 ? '+' : '-',chr,right_splice->donorpos,right_splice->acceptorpos,right_splice->count);
#endif
	      } else if (splice->sign != left_splice->sign) {
		fprintf(stderr,"Splice %c%s:%u..%u (%d) eliminated by %c%s:%u..%u (%d) and  %c%s:%u..%u (%d)\n",
			splice->sign > 0 ? '+' : '-',chr,splice->donorpos,splice->acceptorpos,splice->count,
			left_splice->sign > 0 ? '+' : '-',chr,left_splice->donorpos,left_splice->acceptorpos,left_splice->count,
			right_splice->sign > 0 ? '+' : '-',chr,right_splice->donorpos,right_splice->acceptorpos,right_splice->count);
		splice->boundedp = false;
	      }
	    }
	  }
	}
      }
    }
  }

  FREE(omega_array);
  FREE(alpha_array);

  return;
}
#endif



static char *
get_donor_info (int *exoni, int *nexons, char *annot, Genomicpos_T chrpos) {
  char *genename;
  Genomicpos_T exonstart, exonend;
  char *p, *q;
  int gene_namelength, k;
  bool foundp = false, firstp = true, lastp = false;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }

  gene_namelength = p - annot;
  genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
  strncpy(genename,annot,gene_namelength);
  for (q = genename, k = 0; k < gene_namelength; q++, k++) {
    if (*q == '-') {
      *q = '.';
    }
  }

  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;


  *nexons = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      *nexons += 1;

      /* Advance to beginning of next exon to see if this is the last one */
      while (*p != '\0' && *p != '\n') p++;
      if (*p == '\n') p++;
      if (*p == '\0') {
	lastp = true;
      }

      if (exonend == chrpos && lastp == false) {
	*exoni = *nexons;
	foundp = true;
      }

      firstp = false;
    }
  }

  if (foundp == false) {
    FREE(genename);
    return (char *) NULL;
  } else {
    while (*p != '\0') {
      if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
	fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
	abort();
      } else {
	*nexons += 1;

	/* Advance to beginning of next exon to see if this is the last one */
	while (*p != '\0' && *p != '\n') p++;
	if (*p == '\n') p++;
#if 0
	if (*p == '\0') {
	  lastp = true;
	}
#endif
      }
    }
    return genename;
  }
}


static char *
get_acceptor_info (int *exoni, int *nexons, char *annot, Genomicpos_T chrpos) {
  char *genename;
  Genomicpos_T exonstart, exonend;
  char *p, *q;
  int gene_namelength, k;
  bool foundp = false, firstp = true, lastp = false;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }

  gene_namelength = p - annot;
  genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
  strncpy(genename,annot,gene_namelength);
  for (q = genename, k = 0; k < gene_namelength; q++, k++) {
    if (*q == '-') {
      *q = '.';
    }
  }

  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;


  *nexons = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      *nexons += 1;

      /* Advance to beginning of next exon to see if this is the last one */
      while (*p != '\0' && *p != '\n') p++;
      if (*p == '\n') p++;
#if 0
      if (*p == '\0') {
	lastp = true;
      }
#endif

      if (exonstart == chrpos && firstp == false) {
	*exoni = *nexons;
	foundp = true;
      }

      firstp = false;
    }
  }

  if (foundp == false) {
    FREE(genename);
    return (char *) NULL;
  } else {
    while (*p != '\0') {
      if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
	fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
	abort();
      } else {
	*nexons += 1;

	/* Advance to beginning of next exon to see if this is the last one */
	while (*p != '\0' && *p != '\n') p++;
	if (*p == '\n') p++;
#if 0
	if (*p == '\0') {
	  lastp = true;
	}
#endif

      }
    }
    return genename;
  }
}


static int
string_cmp (const void *a, const void *b) {
  char *x = * (char **) a;
  char *y = * (char **) b;

  return strcmp(x,y);
}


static void
list_gc (List_T *names) {
  List_T p;
  char *name;

  for (p = *names; p != NULL; p = List_next(p)) {
    name = (char *) List_head(p);
    FREE(name);
  }
  List_free(&(*names));
  return;
}


static void
print_list_unique (List_T names) {
  char **array;
  int n, i;

  if (names != NULL) {
    array = (char **) List_to_array_n(&n,names);
    qsort(array,n,sizeof(char *),string_cmp);
    printf("%s",array[0]);
    for (i = 1; i < n; i++) {
      if (!strcmp(array[i],array[i-1])) {
	/* Skip */
      } else {
	printf("|%s",array[i]);
      }
    }
    FREE(array);
  }

  return;
}

static void
print_list_limited (List_T names, int nmax) {
  char **array;
  int n, i;

  if (names != NULL) {
    array = (char **) List_to_array_n(&n,names);
    qsort(array,n,sizeof(char *),string_cmp);
    printf("%s",array[0]);
    for (i = 1; i < n && i < nmax; i++) {
      printf("|%s",array[i]);
    }
    FREE(array);
  }

  return;
}


void
Splice_print (T splice, char *chr, IIT_T genes_iit, bool show_invalid_p) {
  List_T donor_genenames = NULL, acceptor_genenames = NULL, donor_exonnames = NULL, acceptor_exonnames = NULL, p;
  char *donor_genename, *acceptor_genename, *exonname, Buffer[1024];
  int donor_exoni, donor_nexons, acceptor_exoni, acceptor_nexons;
  char *acc, *annot, *restofheader;
  int *matches, nmatches, i;
  bool allocp, alloc2p;


  if (show_invalid_p == false && splice->status != VALID) {
    /* Don't print */
  } else {
    printf(">");
    if (genes_iit == NULL) {
    
    } else {
      matches = IIT_get(&nmatches,genes_iit,chr,splice->donorpos,splice->donorpos,/*sortp*/false);
      for (i = 0; i < nmatches; i++) {
	annot = IIT_annotation(&restofheader,genes_iit,/*index*/matches[i],&allocp);
	if ((donor_genename = get_donor_info(&donor_exoni,&donor_nexons,annot,splice->donorpos)) != NULL) {
	  donor_genenames = List_push(donor_genenames,(void *) donor_genename);
	  acc = IIT_label(genes_iit,matches[i],&alloc2p);
	  sprintf(Buffer,"%s_exon%d/%d",acc,donor_exoni,donor_nexons);
	  exonname = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
	  strcpy(exonname,Buffer);
	  donor_exonnames = List_push(donor_exonnames,exonname);
	  if (alloc2p) {
	    FREE(acc);
	  }
	}
	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(matches);

      matches = IIT_get(&nmatches,genes_iit,chr,splice->acceptorpos,splice->acceptorpos,/*sortp*/false);
      for (i = 0; i < nmatches; i++) {
	annot = IIT_annotation(&restofheader,genes_iit,/*index*/matches[i],&allocp);
	if ((acceptor_genename = get_acceptor_info(&acceptor_exoni,&acceptor_nexons,annot,splice->acceptorpos)) != NULL) {
	  acceptor_genenames = List_push(acceptor_genenames,(void *) acceptor_genename);
	  acc = IIT_label(genes_iit,matches[i],&alloc2p);
	  sprintf(Buffer,"%s_exon%d/%d",acc,acceptor_exoni,acceptor_nexons);
	  exonname = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
	  strcpy(exonname,Buffer);
	  acceptor_exonnames = List_push(acceptor_exonnames,exonname);
	  if (alloc2p) {
	    FREE(acc);
	  }
	}
	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(matches);
    }

    if (donor_genenames == NULL) {
      printf("Novel");
    } else {
      print_list_unique(donor_genenames);
      list_gc(&donor_genenames);
    }
    printf("-");
    if (acceptor_genenames == NULL) {
      printf("Novel");
    } else {
      print_list_unique(acceptor_genenames);
      list_gc(&acceptor_genenames);
    }
    printf("-");
    if (donor_exonnames == NULL) {
      printf("NA");
    } else {
      print_list_limited(donor_exonnames,5);
      list_gc(&donor_exonnames);
    }
    printf("-");
    if (acceptor_exonnames == NULL) {
      printf("NA");
    } else {
      print_list_limited(acceptor_exonnames,5);
      list_gc(&acceptor_exonnames);
    }

    printf(" %s:%u..%u nprimary:%d ncrosshyb:%d nunique:%d",
	   chr,splice->donorpos,splice->acceptorpos,
	   splice->primary_count,splice->crosshyb_count,splice->nunique);
    printf(" maxminsupport:%d",splice->maxminsupport);
    if (splice->status != VALID) {
      printf(" invalid:%s",Splicestatus_string(splice->status));
    }
    if (splice->canonicalp == false) {
      printf(" %c%c-%c%c",splice->donor1,splice->donor2,splice->acceptor1,splice->acceptor2);
    }
    printf("\n");
  }

  return;
}

