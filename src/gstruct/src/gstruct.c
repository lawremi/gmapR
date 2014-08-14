static char rcsid[] = "$Id: gstruct.c 138418 2014-06-06 21:12:37Z twu $";
/* Note: Handles only paired-end data */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gstruct.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For qsort */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "chrom.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "iit-write.h"
#include "interval.h"
#include "table.h"
#include "tableint.h"
#include "tableuint.h"
#include "uinttable.h"
#include "maxent_hr.h"

#ifdef BAM_INPUT
#include "bamread.h"
#endif

#include "samflags.h"
#include "samread.h"
#include "genome.h"
#include "complement.h"
#include "splice.h"


#define MAX_INSERTLENGTH 1000
#define TEN_MILLION 10000000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Parsing SAM */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Adding splices and sites to graph */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Monitoring progress */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif



/************************************************************************
 *   Extents
 ************************************************************************/

typedef struct Extent_T *Extent_T;
struct Extent_T {
  int count;
  int nunique;

  Genomicpos_T low;
  Genomicpos_T high;
};

static void
Extent_free (Extent_T *old) {
  FREE(*old);
  return;
}


#if 0
static int
Extent_cmp (const void *a, const void *b) {
  Extent_T x = * (Extent_T *) a;
  Extent_T y = * (Extent_T *) b;

  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;
  } else {
    return 0;
  }
}
#endif

#if 0
static int
Extent_high_cmp (const void *a, const void *b) {
  Extent_T x = * (Extent_T *) a;
  Extent_T y = * (Extent_T *) b;

  if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;
  } else {
    return 0;
  }
}
#endif

#if 0
static int
Extent_count_cmp (const void *a, const void *b) {
  Extent_T x = * (Extent_T *) a;
  Extent_T y = * (Extent_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (y->count >  x->count) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
static List_T
Extent_cumulate_counts (List_T extents) {
  long int total_up, total_down;
  List_T p;
  Extent_T extent;

  total_up = 0;
  for (p = extents; p != NULL; p = List_next(p)) {
    extent = (Extent_T) List_head(p);
    total_up += extent->count;
    extent->total_up = total_up;
  }

  extents = List_reverse(extents);

  total_down = 0;
  for (p = extents; p != NULL; p = List_next(p)) {
    extent = (Extent_T) List_head(p);
    total_down += extent->count;
    extent->total_down = total_down;
  }

  return List_reverse(extents);
}
#endif


static List_T
Extent_add_at_low (List_T list, Genomicpos_T low, Genomicpos_T high, int nhits,
		   Genomicpos_T max_pairlength) {
  List_T p;
  Extent_T extent;

  debug1(printf("List has %d elements currently\n",List_length(list)));

  if (high - low > max_pairlength) {
    debug1(printf("pair is too long, so skipping\n"));
    return list;

  } else {
    for (p = list; p != NULL; p = List_next(p)) {
      extent = (Extent_T) List_head(p);
      debug1(printf("  Comparing high %u with list high %u",high,extent->high));
      if (extent->high == high) {
	extent->count += 1;
	if (nhits == 1) {
	  extent->nunique += 1;
	}
	debug1(printf("  equal, so count is now %d\n",extent->count));
	return list;
      }
      debug1(printf("\n"));
    }

    /* Not found, so add to list */
    debug1(printf("High not found, so creating new extent\n"));
    extent = (Extent_T) MALLOC(sizeof(*extent));
    extent->count = 1;
    extent->nunique = 0;
    if (nhits == 1) {
      extent->nunique += 1;
    }
    extent->low = low;
    extent->high = high;

    return List_push(list,extent);
  }
}


/************************************************************************
 *   Sites
 ************************************************************************/

typedef struct Site_T *Site_T;
struct Site_T {
  Genomicpos_T chrpos;
  List_T intervals;
};


static void
Site_free_extent (Site_T *old) {
  List_T p;
  Extent_T extent;

  for (p = (*old)->intervals; p != NULL; p = List_next(p)) {
    extent = (Extent_T) List_head(p);
    Extent_free(&extent);
  }
  List_free(&(*old)->intervals);

  FREE(*old);
  return;
}

static void
Site_free_splice (Site_T *old) {
  List_T p;
  Splice_T splice;

  /* Freed by procedure creating splices */
  for (p = (*old)->intervals; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    Splice_free(&splice);
  }
  List_free(&(*old)->intervals);

  FREE(*old);
  return;
}


static Site_T
Site_new (Genomicpos_T chrpos) {
  Site_T new = (Site_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->intervals = (List_T) NULL;

  return new;
}


#if 0
static void
Site_print_extent (Site_T this, char *chr) {
  Extent_T extent, *array;
  int n, i;
  
  n = List_length(this->intervals);
  if (n > 0) {
    array = (Extent_T *) List_to_array(this->intervals,NULL);
    qsort(array,n,sizeof(Extent_T),Extent_cmp);
    for (i = 0; i < n; i++) {
      extent = array[i];
      printf(">%d %s:%u..%u nunique:%d",
	     extent->count,chr,extent->low,extent->high,
	     extent->nunique);
      printf("\n");

    }
    FREE(array);
  }

  return;
}
#endif


#if 0
static void
Site_print_splice (Site_T this, char *chr, int mincount, int minsupport) {
  Splice_T *array;
  int n, i;
  
  n = List_length(this->intervals);
  if (n > 0) {
    array = (Splice_T *) List_to_array(this->intervals,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_cmp);
    for (i = 0; i < n; i++) {
      Splice_print(array[i],chr,mincount,minsupport,/*show_invalid_p*/true);
    }
    FREE(array);
  }

  return;
}
#endif



static void
Site_add_at_low (Table_T chrtable, char *chr, Genomicpos_T low, Genomicpos_T high,
		 Genomicpos_T max_pairlength, int nhits) {
  Site_T site;
  Uinttable_T sitetable;
  Chrom_T chrom;
  
  debug1(printf("Called Site_add_at_low with chr %s, low %u, high %u\n",
		chr,low,high));

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((sitetable = (Uinttable_T) Table_get(chrtable,(void *) chrom)) == NULL) {
    debug1(printf("Made new sitetable for chr %s\n",chr));
    sitetable = Uinttable_new(65522); /* estimate 65522 splice sites per chromosome */
    Table_put(chrtable,(void *) chrom,(void *) sitetable);
  } else {
    Chrom_free(&chrom);
  }

  if ((site = (Site_T) Uinttable_get(sitetable,low)) == NULL) {
    debug1(printf("Made new site\n"));
    site = Site_new(low);
    Uinttable_put(sitetable,low,(void *) site);
  }

  site->intervals = Extent_add_at_low(site->intervals,low,high,nhits,max_pairlength);

  return;
}


static List_T
accumulate_extents (List_T extents, Uinttable_T sitetable) {
  Genomicpos_T *keys;
  Site_T site;
  List_T p;
  Extent_T extent;
  int n, i;

  if (sitetable == NULL) {
    return extents;
  } else {
    n = Uinttable_length(sitetable);
    keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      for (p = site->intervals; p != NULL; p = List_next(p)) {
	extent = (Extent_T) List_head(p);
	extents = List_push(extents,(void *) extent);
      }
    }
    FREE(keys);
  
    return extents;
  }
}



static void
Site_add_at_donor (Table_T chrtable, char *chr, int sign,
		   Genomicpos_T donorpos, Genomicpos_T acceptorpos,
		   bool crosshybp, int nhits, int support1, int support2,
		   int overhang1, int overhang2, int querylength, bool lowend_p,
		   bool canonicalp, double donorprob, double acceptorprob,
		   char donor1, char donor2, char acceptor1, char acceptor2, bool knownp) {
  Site_T site;
  Uinttable_T sitetable;
  Chrom_T chrom;
  
  debug1(printf("Called Site_add_at_donor with chr %s, sign %d, donorpos %u, acceptorpos %u\n",
		chr,sign,donorpos,acceptorpos));

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((sitetable = (Uinttable_T) Table_get(chrtable,(void *) chrom)) == NULL) {
    debug1(printf("Made new sitetable for chr %s\n",chr));
    sitetable = Uinttable_new(65522); /* estimate 65522 splice sites per chromosome */
    Table_put(chrtable,(void *) chrom,(void *) sitetable);
  } else {
    Chrom_free(&chrom);
  }

  if ((site = (Site_T) Uinttable_get(sitetable,donorpos)) == NULL) {
    debug1(printf("Made new site\n"));
    site = Site_new(donorpos);
    Uinttable_put(sitetable,donorpos,(void *) site);
  }

  if (donorpos == 0 || acceptorpos == 0) {
    fprintf(stderr,"Unexpected donorpos %u or acceptorpos %u in sam_splices or gsnap_splices\n",
	    donorpos,acceptorpos);
    abort();
  }

  site->intervals = Splice_add_at_donor(site->intervals,sign,donorpos,acceptorpos,crosshybp,nhits,
					support1,support2,overhang1,overhang2,querylength,lowend_p,canonicalp,
					donorprob,acceptorprob,donor1,donor2,acceptor1,acceptor2,
					knownp);

  return;
}


static List_T
accumulate_splices (List_T splices, Uinttable_T sitetable) {
  Genomicpos_T *keys;
  Site_T site;
  List_T p;
  Splice_T splice;
  int n, i;

  if (sitetable == NULL) {
    return splices;
  } else {
    n = Uinttable_length(sitetable);
    keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      for (p = site->intervals; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (Splice_validp(splice) == true) {
	  splices = List_push(splices,(void *) splice);
	}
      }
    }
    FREE(keys);
  
    return splices;
  }
}



/************************************************************************
 *   GSTRUCT
 ************************************************************************/

#define T Gstruct_T
struct T {
  int npairs;

  Table_T fwd_extents_chrtable;
  Table_T rev_extents_chrtable;
  Table_T null_extents_chrtable;

  Table_T primary_extents_chrtable;
  Table_T crosshyb_extents_chrtable;

  Table_T fwd_splices_chrtable;
  Table_T rev_splices_chrtable;

#ifdef BAM_INPUT
  Tableuint_T chrlength_table;
  Table_T bamstore_chrtable;
#endif

  Chrom_T *fwd_chroms;
  Chrom_T *rev_chroms;
  
  int fwd_nchroms;
  int rev_nchroms;

  int fwd_chromi;
  int rev_chromi;
};


static void
sitetable_extent_gc (Uinttable_T *sitetable) {
  int n, i;
  Genomicpos_T *keys;
  Site_T site;

  n = Uinttable_length(*sitetable);
  keys = Uinttable_keys(*sitetable,/*sortp*/false);
  for (i = 0; i < n; i++) {
    site = (Site_T) Uinttable_get(*sitetable,keys[i]);
    Site_free_extent(&site);
  }

  FREE(keys);
  Uinttable_free(&(*sitetable));
  return;
}

static void
sitetable_splice_gc (Uinttable_T *sitetable) {
  int n, i;
  Genomicpos_T *keys;
  Site_T site;

  n = Uinttable_length(*sitetable);
  keys = Uinttable_keys(*sitetable,/*sortp*/false);
  for (i = 0; i < n; i++) {
    site = (Site_T) Uinttable_get(*sitetable,keys[i]);
    Site_free_splice(&site);
  }

  FREE(keys);
  Uinttable_free(&(*sitetable));
  return;
}


void
Gstruct_free (T *old) {
  int n, i;
  Uinttable_T sitetable;
  Chrom_T *chroms, chrom;
#ifdef BAM_INPUT
  Uinttable_T bamstore_table;
#endif

  /* Free fwd splices */
  for (i = 0; i < (*old)->fwd_nchroms; i++) {
    chrom = (*old)->fwd_chroms[i];
    sitetable = Table_get((*old)->fwd_splices_chrtable,(void *) chrom);
    sitetable_splice_gc(&sitetable);
  }
  for (i = 0; i < (*old)->fwd_nchroms; i++) {
    Chrom_free(&((*old)->fwd_chroms[i]));
  }
  FREE((*old)->fwd_chroms);
  Table_free(&(*old)->fwd_splices_chrtable);

  /* Free rev splices */
  for (i = 0; i < (*old)->rev_nchroms; i++) {
    chrom = (*old)->rev_chroms[i];
    sitetable = Table_get((*old)->rev_splices_chrtable,(void *) chrom);
    sitetable_splice_gc(&sitetable);
  }
  for (i = 0; i < (*old)->rev_nchroms; i++) {
    Chrom_free(&((*old)->rev_chroms[i]));
  }
  FREE((*old)->rev_chroms);
  Table_free(&(*old)->rev_splices_chrtable);


  /* Free primary extents */
  if ((*old)->primary_extents_chrtable != NULL) {
    if ((n = Table_length((*old)->primary_extents_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->primary_extents_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	sitetable = Table_get((*old)->primary_extents_chrtable,(void *) chrom);
	sitetable_extent_gc(&sitetable);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->primary_extents_chrtable);
  }

  /* Free crosshyb extents */
  if ((*old)->crosshyb_extents_chrtable != NULL) {
    if ((n = Table_length((*old)->crosshyb_extents_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->crosshyb_extents_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	sitetable = Table_get((*old)->crosshyb_extents_chrtable,(void *) chrom);
	sitetable_extent_gc(&sitetable);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->crosshyb_extents_chrtable);
  }

  /* Free fwd extents */
  if ((*old)->fwd_extents_chrtable != NULL) {
    if ((n = Table_length((*old)->fwd_extents_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->fwd_extents_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	sitetable = Table_get((*old)->fwd_extents_chrtable,(void *) chrom);
	sitetable_extent_gc(&sitetable);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->fwd_extents_chrtable);
  }

  /* Free rev extents */
  if ((*old)->rev_extents_chrtable != NULL) {
    if ((n = Table_length((*old)->rev_extents_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->rev_extents_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	sitetable = Table_get((*old)->rev_extents_chrtable,(void *) chrom);
	sitetable_extent_gc(&sitetable);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->rev_extents_chrtable);
  }

  /* Free null extents */
  if ((*old)->null_extents_chrtable != NULL) {
    if ((n = Table_length((*old)->null_extents_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->null_extents_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	sitetable = Table_get((*old)->null_extents_chrtable,(void *) chrom);
	sitetable_extent_gc(&sitetable);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->null_extents_chrtable);
  }


#ifdef BAM_INPUT
  if ((*old)->chrlength_table != NULL) {
    if ((n = Tableuint_length((*old)->chrlength_table)) > 0) {
      chroms = (Chrom_T *) Tableuint_keys((*old)->chrlength_table,NULL);
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Tableuint_free(&(*old)->chrlength_table);
  }

  if ((*old)->bamstore_chrtable != NULL) {
    if ((n = Table_length((*old)->bamstore_chrtable)) > 0) {
      chroms = (Chrom_T *) Table_keys((*old)->bamstore_chrtable,NULL);
      for (i = 0; i < n; i++) {
	chrom = chroms[i];
	bamstore_table = (Uinttable_T) Table_get((*old)->bamstore_chrtable,(void *) chrom);
	Bamstore_table_free(&bamstore_table);
	Uinttable_free(&bamstore_table);
      }
      for (i = 0; i < n; i++) {
	Chrom_free(&(chroms[i]));
      }
      FREE(chroms);
    }
    Table_free(&(*old)->bamstore_chrtable);
  }
#endif

  FREE(*old);
  return;
}


int
Gstruct_npairs (T this) {
  return this->npairs;
}



/************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static void
find_dinucleotides (bool *canonicalp, double *donorprob, double *acceptorprob,
		    char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		    Genomicpos_T firstpos, Genomicpos_T secondpos, char *chr,
		    Genome_T genome, IIT_T chromosome_iit, char truestrand, bool trust_sam_p) {
  Chrnum_T chrnum;
  Genomicpos_T chroffset;
  char nt1, nt2, nt3, nt4;

  chrnum = IIT_find_one(chromosome_iit,chr);
  chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum)) - 1U;

  if (truestrand != ' ') {
    if (trust_sam_p == true) {
      *canonicalp = true;
      *donor1 = *donor2 = *acceptor1 = *acceptor2 = ' ';

    } else {
      /* Look at genome inside of firstpos and secondpos to get dinucleotides */
      nt1 = Genome_get_char(genome,chroffset+firstpos+1);
      nt2 = Genome_get_char(genome,chroffset+firstpos+2);
      nt3 = Genome_get_char(genome,chroffset+secondpos-2);
      nt4 = Genome_get_char(genome,chroffset+secondpos-1);

      debug(printf("Got splice from %u to %u\n",firstpos,secondpos));
      debug(printf("Dinucleotides are %c%c to %c%c\n",nt1,nt2,nt3,nt4));

      if (truestrand == '+') {
	if (nt1 == 'G' && (nt2 == 'T' || nt2 == 'C') && nt3 == 'A' && nt4 == 'G') {
	  *canonicalp = true;
	} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A' && nt4 == 'C') {
	  *canonicalp = true;
	} else {
	  *canonicalp = false;
	}
	*donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;

      } else if (truestrand == '-') {
	if (nt1 == 'C' && nt2 == 'T' && (nt3 == 'A' || nt3 == 'G') && nt4 == 'C') {
	  *canonicalp = true;
	} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A' && nt4 == 'T') {
	  *canonicalp = true;
	} else {
	  *canonicalp = false;
	}
	*donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];

      } else {
	fprintf(stderr,"Unrecognized truestrand %c\n",truestrand);
	abort();
      }

    }

  } else {
    /* Look at genome inside of firstpos and secondpos to determine truestrand */
    nt1 = Genome_get_char(genome,chroffset+firstpos+1);
    nt2 = Genome_get_char(genome,chroffset+firstpos+2);
    nt3 = Genome_get_char(genome,chroffset+secondpos-2);
    nt4 = Genome_get_char(genome,chroffset+secondpos-1);

    debug(printf("Got splice from %u to %u\n",firstpos,secondpos));
    debug(printf("Dinucleotides are %c%c to %c%c\n",nt1,nt2,nt3,nt4));

    if (nt1 == 'G' && (nt2 == 'T' || nt2 == 'C') && nt3 == 'A' && nt4 == 'G') {
      truestrand = '+';
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
    } else if (nt1 == 'C' && nt2 == 'T' && (nt3 == 'A' || nt3 == 'G') && nt4 == 'C') {
      truestrand = '-';
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
    } else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A' && nt4 == 'C') {
      truestrand = '+';
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
    } else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A' && nt4 == 'T') {
      truestrand = '-';
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
    } else {
      /* In GSNAP, will want to output sense information in SAM output. */
#if 0
      fprintf(stderr,"Splice %s:%u..%u is not (semi-)canonical: %c%c...%c%c.  Cannot determine sense.\n",
	      chr,firstpos,secondpos,nt1,nt2,nt3,nt4);
#endif
      truestrand = ' ';
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = false;
    }
  }

  if (truestrand == '+') {
    *donorprob = Maxent_hr_donor_prob(chroffset+firstpos+1U);
    *acceptorprob = Maxent_hr_acceptor_prob(chroffset+secondpos);
    /* printf("fwd splice probs at %s:%u..%u are %f and %f\n",
       chr,firstpos,secondpos,*donorprob,*acceptorprob); */
  } else if (truestrand == '-') {
    *donorprob = Maxent_hr_antidonor_prob(chroffset+secondpos);
    *acceptorprob = Maxent_hr_antiacceptor_prob(chroffset+firstpos+1U);
    /* printf("rev splice probs at %s:%u..%u are %f and %f\n",
       chr,secondpos,firstpos,*donorprob,*acceptorprob); */
  } else {
    *donorprob = 0.0;
    *acceptorprob = 0.0;
  }

  return;
}


static void
add_splice (Genomicpos_T firstpos, Genomicpos_T secondpos,
	    char *chr, bool crosshybp, int nhits, int support1, int support2,
	    int overhang1, int overhang2, int querylength, bool lowend_p,
	    Table_T fwd_splices_chrtable, Table_T rev_splices_chrtable,
	    Genome_T genome, IIT_T chromosome_iit, char truestrand,
	    Genomicpos_T shortsplicedist, bool trust_sam_p, bool need_canonical_p,
	    bool knownp) {
  bool canonicalp;
  double donorprob, acceptorprob;
  char donor1, donor2, acceptor1, acceptor2;
  Genomicpos_T donorpos, acceptorpos;

  if (truestrand == '+') {
    donorpos = firstpos;
    acceptorpos = secondpos;
  } else if (truestrand == '-') {
    donorpos = secondpos;
    acceptorpos = firstpos;
  } else {
    return;
  }

  find_dinucleotides(&canonicalp,&donorprob,&acceptorprob,
		     &donor1,&donor2,&acceptor1,&acceptor2,
		     firstpos,secondpos,chr,genome,chromosome_iit,truestrand,
		     trust_sam_p);

  debug0(printf("Got splice from %u to %u\n",donorpos,acceptorpos));
  if (truestrand == '+' && acceptorpos > donorpos + shortsplicedist) {
    /* Distance too far */
  } else if (truestrand == '-' && donorpos > acceptorpos + shortsplicedist) {
    /* Distance too far */
  } else if (truestrand == ' ' && secondpos > firstpos + shortsplicedist) {
    /* Distance too far */
  } else {
    debug(printf("%c%s:%u..%u\n",truestrand,chr,donorpos,acceptorpos));
    
    if (truestrand == '+') {
      Site_add_at_donor(fwd_splices_chrtable,chr,/*sign*/+1,donorpos,acceptorpos,
			crosshybp,nhits,support1,support2,overhang1,overhang2,querylength,
			lowend_p,canonicalp,donorprob,acceptorprob,
			donor1,donor2,acceptor1,acceptor2,knownp);
    } else if (truestrand == '-') {
      Site_add_at_donor(rev_splices_chrtable,chr,/*sign*/-1,donorpos,acceptorpos,
			crosshybp,nhits,support1,support2,overhang1,overhang2,querylength,
			lowend_p,canonicalp,donorprob,acceptorprob,
			donor1,donor2,acceptor1,acceptor2,knownp);
    } else if (truestrand == ' ') {
      if (need_canonical_p == false) {
	fprintf(stderr,"Adding splice at %s:%u..%u to both strands\n",chr,firstpos,secondpos);
	/* Add to both sides */
	Site_add_at_donor(fwd_splices_chrtable,chr,/*sign*/+1,/*donorpos*/firstpos,/*acceptorpos*/secondpos,
			  crosshybp,nhits,support1,support2,overhang1,overhang2,querylength,
			  lowend_p,canonicalp,donorprob,acceptorprob,
			  donor1,donor2,acceptor1,acceptor2,knownp);
	Site_add_at_donor(rev_splices_chrtable,chr,/*sign*/-1,/*donorpos*/secondpos,/*acceptorpos*/firstpos,
			  crosshybp,nhits,support1,support2,overhang1,overhang2,querylength,
			  lowend_p,canonicalp,donorprob,acceptorprob,
			  donor1,donor2,acceptor1,acceptor2,knownp);
      }
    }
  }

  return;
}


static int
get_support2 (Intlist_T types, Uintlist_T npositions) {
  int support2 = 0;
  Intlist_T p = types;
  Uintlist_T q = npositions;
  int type;
  
  while (p != NULL) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      support2 += Uintlist_head(q);

    } else if (type == 'N') {
      return support2;

    } else if (type == 'I') {
      support2 += Uintlist_head(q);

    } else if (type == 'D') {
      if (Uintlist_head(q) > 50) {
	return support2;
      } else {
	/* Do nothing */
      }

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }

    p = Intlist_next(p);
    q = Uintlist_next(q);
  }

  return support2;
}


#if 0
/* Also appears in bamtally.c */
static bool
best_mapping_p (Tableuint_T resolve_low_table, Tableuint_T resolve_high_table, char *acc,
		char *chr, Genomicpos_T chrpos_low, IIT_T chromosome_iit) {
  Genomicpos_T genomicpos, genomicpos_low, genomicpos_high, chroffset;
  Interval_T interval;
  int index;

  genomicpos_low = Tableuint_get(resolve_low_table,(void *) acc);
  genomicpos_high = Tableuint_get(resolve_high_table,(void *) acc);

  if (genomicpos_low == 0 && genomicpos_high == 0) {
    /* Was already unique.  Could also check nhits. */
    return true;
  } else if ((index = IIT_find_linear(chromosome_iit,chr)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in genome\n",chr);
    return false;
  } else {
    interval = IIT_interval(chromosome_iit,index);
    chroffset = Interval_low(interval);
    genomicpos = chroffset + chrpos_low;
    if (genomicpos == genomicpos_low || genomicpos == genomicpos_high) {
      /* printf("acc %s is the best mapping at %s:%u\n",acc,chr,chrpos_low); */
      return true;
    } else {
      /* printf("acc %s is not the best mapping at %s:%u\n",acc,chr,chrpos_low); */
      return false;
    }
  }
}
#endif


static bool
determine_crosshybp (int nhits, char *chr, Genomicpos_T chrpos_low, Genomicpos_T chrpos_high,
		     char *acc, Tableint_T ngoodhits_table, IIT_T genes_iit) {
  int *matches, nmatches;
  int ngoodhits;

  if (nhits == 1) {
    return false;

  } else if (ngoodhits_table == NULL) {
    return true;		/* since nhits > 1 */

  } else {
    ngoodhits = Tableint_get(ngoodhits_table,(void *) acc);
    if (ngoodhits == 0) {
      return true;		/* since nhits > 1; all hits are bad */
    } else if (ngoodhits == 1) {
      matches = IIT_get(&nmatches,genes_iit,chr,/*coordstart*/chrpos_low,
			/*coordend*/chrpos_high,/*sortp*/false);
      if (nmatches > 0) {
	FREE(matches);
	return false;	/* This is the primary hit */
      } else {
	return true;	/* Another location is the primary hit */
      }
    } else {
      return true;
    }
  }
}


static void
parse_splices (int nhits, char *chr, Genomicpos_T chrpos_low, bool crosshybp,
	       Intlist_T types, Uintlist_T npositions, int querylength, bool lowend_p,
	       char truestrand, Table_T fwd_splices_chrtable, Table_T rev_splices_chrtable,
	       Genome_T genome, IIT_T chromosome_iit,
	       Genomicpos_T shortsplicedist, bool trust_sam_p, bool need_canonical_p) {
  Genomicpos_T firstpos, secondpos, chrpos;
  int type;
  int support1, support2, overhang1, overhang2;
  Intlist_T p;
  Uintlist_T q;

  /* Get splice coordinates */
  chrpos = chrpos_low;
  support1 = 0;
  overhang1 = 0;
  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */
      overhang1 += Uintlist_head(q);

    } else if (type == 'H') {
      /* Ignore */
      overhang1 += Uintlist_head(q);

    } else if (type == 'M') {
      chrpos += Uintlist_head(q);
      support1 += Uintlist_head(q);
      overhang1 += Uintlist_head(q);

    } else if (type == 'N') {
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;
      support2 = get_support2(/*types*/Intlist_next(p),/*npositions*/Uintlist_next(q));
      overhang2 = querylength - overhang1;

      add_splice(firstpos,secondpos,chr,crosshybp,nhits,support1,support2,
		 overhang1,overhang2,querylength,lowend_p,
		 fwd_splices_chrtable,rev_splices_chrtable,
		 genome,chromosome_iit,truestrand,shortsplicedist,
		 trust_sam_p,need_canonical_p,/*knownp*/false);

      support1 = 0;

    } else if (type == 'I') {
      /* Do nothing */
      support1 += Uintlist_head(q);
      overhang1 += Uintlist_head(q);

    } else if (type == 'D') {
      /* CHECK */
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;
      if (0 && secondpos - firstpos > 50) {
	/* Note: This requires the user to provide the genome */
	support2 = get_support2(/*types*/Intlist_next(p),/*npositions*/Uintlist_next(q));

	add_splice(firstpos,secondpos,chr,crosshybp,nhits,support1,support2,
		   overhang1,overhang2,querylength,lowend_p,
		   fwd_splices_chrtable,rev_splices_chrtable,
		   genome,chromosome_iit,truestrand,
		   shortsplicedist,trust_sam_p,need_canonical_p,/*knownp*/false);

	support1 = 0;
      }

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
    debug(printf("  type = %c, chrpos = %u\n",type,chrpos));
  }


  return;
}


/************************************************************************
 *   Output
 ************************************************************************/

#if 0
static void
print_extents (Uinttable_T sitetable, char *chr) {
  Genomicpos_T *keys;
  int n, i;
  Site_T site;

  if ((n = Uinttable_length(sitetable)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      site = Uinttable_get(sitetable,keys[i]);
      if (site == NULL) {
	fprintf(stderr,"key is %u, value is NULL\n",keys[i]);
	abort();
      } else {
	Site_print_extent(site,chr);
	Site_free_extent(&site);
      }
    }
    FREE(keys);
  }

  return;
}
#endif


#if 0
void
Gstruct_print_splices_sitetable (Uinttable_T sitetable, char *chr, int mincount, int minsupport) {
  Genomicpos_T *keys;
  int n, i;
  Site_T site;

  n = Uinttable_length(sitetable);
  keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    site = (Site_T) Uinttable_get(sitetable,keys[i]);
    if (site == NULL) {
      fprintf(stderr,"key is %u, value is NULL\n",keys[i]);
      abort();
    } else {
      Site_print_splice(site,chr,mincount,minsupport);
      /* Site_free_splice(&site); */
    }
  }
  FREE(keys);

  return;
}
#endif

#if 0
void
Gstruct_print_splices_list (List_T splices, char *chr) {
  List_T p;

  for (p = splices; p != NULL; p = List_next(p)) {
    Splice_print((Splice_T) List_head(p),chr,/*show_invalid_p*/true);
  }

  return;
}
#endif


#if 0
List_T
Gstruct_splices_sitetable_to_intervals (List_T intervals, Uinttable_T sitetable,
					int mincount, int minsupport, bool need_canonical_p) {
  Genomicpos_T *keys;
  int n, i;
  Site_T site;

  List_T p;
  Splice_T splice;

  n = Uinttable_length(sitetable);
  keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    site = (Site_T) Uinttable_get(sitetable,keys[i]);
    if (site == NULL) {
      fprintf(stderr,"key is %u, value is NULL\n",keys[i]);
      abort();
    } else {
      for (p = site->intervals; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (Splice_validp(splice) == true && Splice_count(splice) >= mincount && Splice_maxminsupport(splice) >= minsupport) {
	  if (need_canonical_p == false || Splice_canonicalp(splice) == true) {
	    intervals = List_push(intervals,(void *) Interval_new(Splice_donorpos(splice),Splice_acceptorpos(splice),/*type*/Splice_count(splice)));
	  }
	}
      }
    }
  }
  FREE(keys);

  return intervals;
}
#endif


static void
add_pairings (long int *cum, List_T extents) {
  List_T p;
  Extent_T extent;

  for (p = extents; p != NULL; p = List_next(p)) {
    extent = (Extent_T) List_head(p);
    cum[extent->low] += extent->count;
    /* This is at extent->high + 1, because we want to keep cum at extent->high */
    cum[extent->high + 1] -= extent->count;
  }

  return;
}


#if 0
static void
print_runlengths (Uinttable_T sitetable, char *chr) {
  List_T p;
  Extent_T extent;
  Genomicpos_T chrlength = 0, lastpos, pos;
  long int *cum, level;
  Genomicpos_T *keys;
  Site_T site;
  int n, i;

  if ((n = Uinttable_length(sitetable)) > 0) {
    keys = Uinttable_keys(sitetable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      for (p = site->intervals; p != NULL; p = List_next(p)) {
	extent = (Extent_T) List_head(p);
	if (extent->high + 1 > chrlength) {
	  chrlength = extent->high + 1;
	}
      }
    }

    cum = (long int *) CALLOC(chrlength+1,sizeof(long int));
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      add_pairings(cum,site->intervals);
    }

    FREE(keys);

    /* dump_cum(cum,chrlength); */

    /* Print runlengths */
    lastpos = 1U;
    level = 0;
    for (pos = 1; pos <= chrlength; pos++) {
      if (cum[pos] != 0) {
	/* printf("cum at pos %u is %d\n",pos,cum[pos]); */
	printf(">%ld %s:%u..%u\n",level,chr,lastpos,pos-1);
	level += cum[pos];
	lastpos = pos;
      }
    }

    /* Should print the final level of zero */
    printf(">%ld %s:%u..%u\n",level,chr,lastpos,pos-1);
    if (level != 0) {
      fprintf(stderr,"Ended with a non-zero level\n");
      abort();
    }

    FREE(cum);
  }

  return;
}
#endif


static long int *
store_runlength (Uinttable_T sitetable, Genomicpos_T chrlength) {
  Genomicpos_T obs_chrlength;
  long int *cum, level;
  List_T p;
  Extent_T extent;
  Genomicpos_T lastpos, pos, x;
  Genomicpos_T *keys;
  Site_T site;
  int n, i;

  if (sitetable == NULL || (n = Uinttable_length(sitetable)) == 0) {
    return (long int *) NULL;

  } else {
    obs_chrlength = 0U;
    keys = Uinttable_keys(sitetable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      for (p = site->intervals; p != NULL; p = List_next(p)) {
	extent = (Extent_T) List_head(p);
	if (extent->high + 1 > obs_chrlength) {
	  obs_chrlength = extent->high + 1;
	}
      }
    }

    if (obs_chrlength > chrlength) {
      fprintf(stderr,"observed chrlength %u > chrlength %u\n",obs_chrlength,chrlength);
      chrlength = obs_chrlength;
    }

    cum = (long int *) CALLOC(chrlength+1,sizeof(long int));
    for (i = 0; i < n; i++) {
      site = (Site_T) Uinttable_get(sitetable,keys[i]);
      add_pairings(cum,site->intervals);
    }

    FREE(keys);

    /* Convert to runlengths */
    lastpos = 1U;
    level = 0;
    /* cum[0] = 0; -- Should not be necessary with CALLOC */
    for (pos = 1; pos <= obs_chrlength; pos++) {
      if (cum[pos] != 0) {
	/* printf("cum at pos %u is %d\n",pos,cum[pos]); */
	/* printf(">%ld %s:%u..%u\n",level,chr,lastpos,pos-1); */
	for (x = lastpos; x < pos; x++) {
	  cum[x] = level;
	}

	level += cum[pos];
	lastpos = pos;
      }
    }

    /* Should print the final level of zero */
    if (level != 0) {
      fprintf(stderr,"Ended with a non-zero level\n");
      abort();
    } else {
      /* printf(">%ld %s:%u..%u\n",level,chr,lastpos,pos-1); */
      for (x = lastpos; x < pos; x++) {
	cum[x] = level;
      }
    }

    return cum;
  }
}


#if 0
static long int *
add_runlengths (long int *x, long int *y, Genomicpos_T chrlength) {
  long int *sum;
  Genomicpos_T pos;

  sum = (long int *) CALLOC(chrlength+1,sizeof(long int));
  for (pos = 0; pos <= chrlength; pos++) {
    sum[pos] = x[pos] + y[pos];
  }
  return sum;
}
#endif



/************************************************************************
 *   Input
 ************************************************************************/

static char
combine_strands (char strand1, char strand2) {

  if (strand1 == '?' || strand2 == '?') {
    return '?';

  } else if (strand1 == ' ' && strand2 == ' ') {
    return ' ';

  } else if (strand1 == '+' && strand2 == '-') {
    /* Skip.  Inconsistent directions */
    return '?';

  } else if (strand1 == '-' && strand2 == '+') {
    /* Skip.  Inconsistent directions */
    return '?';

#if 0
  } else if (strcmp(chr1,chr2)) {
    /* Skip.  Inconsistent chromosomes */
    return '?';
#endif

  } else if (strand1 == ' ') {
    return strand2;

  } else if (strand2 == ' ') {
    return strand1;

  } else {
    return strand1;
  }
}



#ifdef BAM_INPUT


/* Copied to Multimap_resolve */
T
Gstruct_bam_input (int *mean_readlength, int *mean_insertlength, List_T bamfiles,
		   char *chromosome, Genomicpos_T chrstart, Genomicpos_T chrend,
		   Tableint_T ngoodhits_low_table, Tableint_T ngoodhits_high_table,
		   IIT_T genes_iit, Genomicpos_T shortsplicedist, Genomicpos_T max_pairlength,
		   Genome_T genome, IIT_T chromosome_iit,
		   char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		   bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		   bool trust_sam_p, bool need_canonical_p, char *bam_lacks_chr, int bam_lacks_chr_length) {
  T new = (T) MALLOC(sizeof(*new));

  List_T p;
  char *bamfile;
  Bamreader_T bamreader;
  unsigned int chrlength;

  Bamline_T bamline, bamline_low;
  int nhits, nhits1, nhits2, flag;
  bool crosshybp1, crosshybp2;
  char *chr, *chrptr, *chrptr_low, strand;
  Chrom_T chrom;

  long int total_readlength = 0;
  int npairs = 0;
  int *insertlength_freq, max_insertlength_freq;
  int insertlength, best_insertlength;


  new->fwd_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->rev_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->null_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  new->primary_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->crosshyb_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  new->fwd_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->rev_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  new->bamstore_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->chrlength_table = Tableuint_new(100,Chrom_compare_table,Chrom_hash_table);

  insertlength_freq = (int *) CALLOC(MAX_INSERTLENGTH+1,sizeof(int));


  for (p = bamfiles; p != NULL; p = List_next(p)) {
    bamfile = (char *) List_head(p);
    bamreader = Bamread_new(bamfile);
    if (chromosome != NULL) {
      Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
    }
    while ((bamline = Bamread_next_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					   need_unique_p,need_primary_p,ignore_duplicates_p,
					   /*need_concordant_p*/true)) != NULL) {

      if ((chr = Bamline_chr(bamline)) == NULL) {
	/* Skip */
	Bamline_free(&bamline);

      } else {
	if (bam_lacks_chr == NULL) {
	  chrptr = chr;
	} else {
	  chrptr = (char *) CALLOC(bam_lacks_chr_length+strlen(chr)+1,sizeof(char));
	  sprintf(chrptr,"%s%s",bam_lacks_chr,chr);
	}

	/* Store chrlength */
	chrom = Chrom_from_string(chrptr,/*mitochondrial_string*/NULL,/*order*/0U);
	if (Tableuint_get(new->chrlength_table,(void *) chrom) == 0) {
	  chrlength = Bamread_chrlength(bamreader,chrptr);
	  Tableuint_put(new->chrlength_table,(void *) chrom,chrlength);
	} else {
	  Chrom_free(&chrom);
	}

	/* Create site tables, even if splices don't exist */
	chrom = Chrom_from_string(chrptr,/*mitochondrial_string*/NULL,/*order*/0U);
	if (Table_get(new->fwd_splices_chrtable,(void *) chrom) == NULL) {
	  Table_put(new->fwd_splices_chrtable,(void *) chrom,(void *) Uinttable_new(65522));
	} else {
	  Chrom_free(&chrom);
	}
	chrom = Chrom_from_string(chrptr,/*mitochondrial_string*/NULL,/*order*/0U);
	if (Table_get(new->rev_splices_chrtable,(void *) chrom) == NULL) {
	  Table_put(new->rev_splices_chrtable,(void *) chrom,(void *) Uinttable_new(65522));
	} else {
	  Chrom_free(&chrom);
	}

	flag = Bamline_flag(bamline);
	if ((flag & PAIRED_READ) == 0) {
	  nhits = Bamline_nhits(bamline);
	  strand = Bamline_strand(bamline,genome,chromosome_iit);

	  if (Bamline_good_unique_p(bamline) == true) {
	    crosshybp2 = false;
	  } else {
	    crosshybp2 = true;
	  }

	  total_readlength = Bamline_readlength(bamline);
	  parse_splices(Bamline_nhits(bamline),chrptr,Bamline_chrpos_low(bamline),
			crosshybp2,Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
			Bamline_cigar_querylength(bamline),Bamline_lowend_p(bamline),strand,
			new->fwd_splices_chrtable,new->rev_splices_chrtable,genome,chromosome_iit,
			shortsplicedist,trust_sam_p,need_canonical_p);

	  if (strand == '+') {
	    Site_add_at_low(new->fwd_extents_chrtable,chrptr,
			    Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),
			    max_pairlength,nhits);
	  } else if (strand == '-') {
	    Site_add_at_low(new->rev_extents_chrtable,chrptr,
			    Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),
			    max_pairlength,nhits);
	  } else if (strand == ' ') {
	    Site_add_at_low(new->null_extents_chrtable,chrptr,
			    Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),
			    max_pairlength,nhits);
	    /* Also compute null_high and null_low */
	  }

	  if (crosshybp2 == false) {
	    Site_add_at_low(new->primary_extents_chrtable,chrptr,
			    Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),
			    max_pairlength,nhits);
	  } else {
	    Site_add_at_low(new->crosshyb_extents_chrtable,chrptr,
			    Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),
			    max_pairlength,nhits);
	  }

	} else if (Bamline_mate_chrpos_low(bamline) > Bamline_chrpos_low(bamline)) {
	  Bamstore_add_at_low(new->bamstore_chrtable,chrptr,Bamline_chrpos_low(bamline),
			      bamline);
	} else {
	  bamline_low = Bamstore_get(new->bamstore_chrtable,chrptr,Bamline_mate_chrpos_low(bamline),
				     Bamline_acc(bamline),Bamline_chrpos_low(bamline));

	  if (bamline_low != NULL) {
	    nhits1 = Bamline_nhits(bamline_low);
	    nhits2 = Bamline_nhits(bamline);
	    if (nhits1 < nhits2) {
	      nhits = nhits1;
	    } else {
	      nhits = nhits2;
	    }

	    strand = combine_strands(Bamline_strand(bamline_low,genome,chromosome_iit),
				     Bamline_strand(bamline,genome,chromosome_iit));

#if 0
	    crosshybp1 = determine_crosshybp(Bamline_nhits(bamline_low),Bamline_chr(bamline_low),Bamline_chrpos_low(bamline_low),
					     Bamline_chrpos_high(bamline_low),Bamline_acc(bamline_low),
					     ngoodhits_low_table,genes_iit);
	    crosshybp2 = determine_crosshybp(Bamline_nhits(bamline),chrptr,Bamline_chrpos_low(bamline),
					     Bamline_chrpos_high(bamline),Bamline_acc(bamline),
					     ngoodhits_high_table,genes_iit);
#else
	    if (Bamline_good_unique_p(bamline_low) == true) {
	      crosshybp1 = false;
	    } else {
	      crosshybp1 = true;
	    }

	    if (Bamline_good_unique_p(bamline) == true) {
	      crosshybp2 = false;
	    } else {
	      crosshybp2 = true;
	    }
#endif

	    total_readlength += Bamline_readlength(bamline_low);
	    total_readlength += Bamline_readlength(bamline);
	    if ((insertlength = Bamline_insert_length(bamline)) < 0) {
	      insertlength = -insertlength;
	    }
	    if (insertlength <= MAX_INSERTLENGTH) {
	      insertlength_freq[insertlength] += 1;
	    }
	    npairs += 1;

	    if (Bamline_chr(bamline_low) == NULL) {
	      chrptr_low = (char *) NULL;
	    } else if (bam_lacks_chr == NULL) {
	      chrptr_low = Bamline_chr(bamline_low);
	    } else {
	      chrptr_low = (char *) CALLOC(bam_lacks_chr_length+strlen(chrptr_low)+1,sizeof(char));
	      sprintf(chrptr_low,"%s%s",bam_lacks_chr,Bamline_chr(bamline_low));
	    }

	    parse_splices(Bamline_nhits(bamline_low),chrptr_low,Bamline_chrpos_low(bamline_low),
			  crosshybp1,Bamline_cigar_types(bamline_low),Bamline_cigar_npositions(bamline_low),
			  Bamline_cigar_querylength(bamline_low),Bamline_lowend_p(bamline_low),strand,
			  new->fwd_splices_chrtable,new->rev_splices_chrtable,genome,chromosome_iit,
			  shortsplicedist,trust_sam_p,need_canonical_p);

	    parse_splices(Bamline_nhits(bamline),chrptr,Bamline_chrpos_low(bamline),
			  crosshybp2,Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
			  Bamline_cigar_querylength(bamline),Bamline_lowend_p(bamline),strand,
			  new->fwd_splices_chrtable,new->rev_splices_chrtable,genome,chromosome_iit,
			  shortsplicedist,trust_sam_p,need_canonical_p);

	    if (strand == '+') {
	      Site_add_at_low(new->fwd_extents_chrtable,chrptr,
			      Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline),
			      max_pairlength,nhits);
	    } else if (strand == '-') {
	      Site_add_at_low(new->rev_extents_chrtable,chrptr,
			      Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline),
			      max_pairlength,nhits);
	    } else if (strand == ' ') {
	      Site_add_at_low(new->null_extents_chrtable,chrptr,
			      Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline),
			      max_pairlength,nhits);
	      /* Also compute null_high and null_low */
	    }

	    if (crosshybp1 == false && crosshybp2 == false) {
	      Site_add_at_low(new->primary_extents_chrtable,chrptr,
			      Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline),
			      max_pairlength,nhits);
	    } else {
	      Site_add_at_low(new->crosshyb_extents_chrtable,chrptr,
			      Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline),
			      max_pairlength,nhits);
	    }

	    Bamline_free(&bamline_low);
	  }
	  Bamline_free(&bamline);
	}
      }
    }

    if (chromosome != NULL) {
      Bamread_unlimit_region(bamreader);
    }
    Bamread_free(&bamreader);
  }


  new->fwd_chroms = (Chrom_T *) Table_keys(new->fwd_splices_chrtable,NULL);
  if ((new->fwd_nchroms = Table_length(new->fwd_splices_chrtable)) > 0) {
    qsort(new->fwd_chroms,new->fwd_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  new->rev_chroms = (Chrom_T *) Table_keys(new->rev_splices_chrtable,NULL);
  if ((new->rev_nchroms = Table_length(new->rev_splices_chrtable)) > 0) {
    qsort(new->rev_chroms,new->rev_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  new->fwd_chromi = new->rev_chromi = 0;

  new->npairs = npairs;

  if (npairs == 0) {
    *mean_readlength = 0;
    *mean_insertlength = 0;
  } else {
    *mean_readlength = (int) rint((double) total_readlength/(double) npairs/2.0);

    best_insertlength = 0;
    max_insertlength_freq = insertlength_freq[0];

    for (insertlength = 1; insertlength <= MAX_INSERTLENGTH; insertlength++) {
      if (insertlength_freq[insertlength] > max_insertlength_freq) {
	best_insertlength = insertlength;
	max_insertlength_freq = insertlength_freq[insertlength];
      }
    }

    *mean_insertlength = best_insertlength;
    /* fprintf(stderr,"Mean readlength: %d\n",*mean_readlength); */
    /* fprintf(stderr,"Mean insertlength: %d\n",*mean_insertlength); */
  }

  FREE(insertlength_freq);

  return new;
}


static long int
get_gene_tally (IIT_T genes_iit, char *chr, Genomicpos_T coordstart, Genomicpos_T coordend) {
  long int maxsum = 0;
  int *sumptr;
  int *matches;
  int nmatches, i;

  matches = IIT_get(&nmatches,genes_iit,chr,coordstart,coordend,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    sumptr = (int *) IIT_data(genes_iit,matches[i]);
    if (sumptr[0] > maxsum) {
      maxsum = sumptr[0];
    }
  }

  FREE(matches);

  return maxsum;
}


/* Taken from Gstruct_bam_input */
void
Gstruct_recount_multimappers (Tableint_T *ngoodhits_low_table, Tableint_T *ngoodhits_high_table, Bamreader_T bamreader,
			      IIT_T genes_iit, int maximum_nhits) {
  Bamline_T bamline, bamline_low;
  int nhits, count;
  Table_T bamstore_chrtable;
  char *chr, *acc, *key;

  Uinttable_T bamstore_table;
  Chrom_T *chroms, chrom;
  int n, i;


  *ngoodhits_low_table = Tableint_new(TEN_MILLION,Table_string_compare,Table_string_hash);
  *ngoodhits_high_table = Tableint_new(TEN_MILLION,Table_string_compare,Table_string_hash);

  bamstore_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/35,maximum_nhits,
					 /*need_unique_p*/false,/*need_primary_p*/false,/*ignore_duplicates_p*/false,
					 /*need_concordant_p*/false)) != NULL) {
    if ((chr = Bamline_chr(bamline)) == NULL) {
      /* Skip */
      Bamline_free(&bamline);

    } else if ((nhits = Bamline_nhits(bamline)) == 1) {
      /* Unique, so skip */
      Bamline_free(&bamline);

    } else if (Bamline_concordantp(bamline) == false || (Bamline_flag(bamline) & PAIRED_READ) == 0) {
      /* Handle now */
      
      acc = Bamline_acc(bamline);
      if (Bamline_firstend_p(bamline) == true) {
	if (IIT_contained(genes_iit,chr,/*coordstart*/Bamline_chrpos_low(bamline),
			  /*coordend*/Bamline_chrpos_high(bamline)) == true) {
	  if ((count = Tableint_get(*ngoodhits_low_table,(void *) acc)) > 0) {
	    key = acc;
	  } else {
	    key = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key,acc);
	  }
	  Tableint_put(*ngoodhits_low_table,(void *) key,count + 1);

#if 0
	  printf("%s %s:%u..%u good => %d\n",
		 acc,chr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),count+1);
#endif
	}

      } else {
	if (IIT_contained(genes_iit,chr,/*coordstart*/Bamline_chrpos_low(bamline),
			  /*coordend*/Bamline_chrpos_high(bamline)) == true) {
	  if ((count = Tableint_get(*ngoodhits_high_table,(void *) acc)) > 0) {
	    key = acc;
	  } else {
	    key = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key,acc);
	  }
	  Tableint_put(*ngoodhits_high_table,(void *) key,count + 1);

#if 0
	  printf("%s %s:%u..%u good => %d\n",
		 acc,chr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),count+1);
#endif
	}
      }

      Bamline_free(&bamline);

    } else if (Bamline_lowend_p(bamline) == true) {
      /* Wait for high end */
      Bamstore_add_at_low(bamstore_chrtable,Bamline_chr(bamline),Bamline_chrpos_low(bamline),
			  bamline);

    } else {
      /* This is the high end */
      acc = Bamline_acc(bamline);
      bamline_low = Bamstore_get(bamstore_chrtable,chr,Bamline_mate_chrpos_low(bamline),
				 acc,Bamline_chrpos_low(bamline));
      if (bamline_low == NULL) {
	fprintf(stderr,"Hmm...low end not found for %s at %s:%u\n",acc,Bamline_chr(bamline),Bamline_chrpos_low(bamline));
      } else {
	if (IIT_contained(genes_iit,chr,/*coordstart*/Bamline_chrpos_low(bamline_low),
			  /*coordend*/Bamline_chrpos_high(bamline_low)) == true) {
	  if ((count = Tableint_get(*ngoodhits_low_table,(void *) acc)) > 0) {
	    key = acc;
	  } else {
	    key = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key,acc);
	  }
	  Tableint_put(*ngoodhits_low_table,(void *) key,count + 1);

#if 0
	  printf("%s %s:%u..%u good => %d => %d\n",
		 acc,chr,Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline_low),count+1,
		 Tableint_get(*ngoodhits_low_table,(void *) acc));
#endif
	}

	if (IIT_contained(genes_iit,chr,/*coordstart*/Bamline_chrpos_low(bamline),
			  /*coordend*/Bamline_chrpos_high(bamline)) == true) {
	  if ((count = Tableint_get(*ngoodhits_high_table,(void *) acc)) > 0) {
	    key = acc;
	  } else {
	    key = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key,acc);
	  }
	  Tableint_put(*ngoodhits_high_table,(void *) key,count + 1);

#if 0
	  printf("%s %s:%u..%u good => %d => %d\n",
		 acc,chr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline),count+1,
		 Tableint_get(*ngoodhits_high_table,(void *) acc));
#endif
	}

	Bamline_free(&bamline_low);
      }

      Bamline_free(&bamline);
    }
  }


  if ((n = Table_length(bamstore_chrtable)) > 0) {
    chroms = (Chrom_T *) Table_keys(bamstore_chrtable,NULL);
    for (i = 0; i < n; i++) {
      chrom = chroms[i];
      bamstore_table = (Uinttable_T) Table_get(bamstore_chrtable,(void *) chrom);
      Bamstore_table_free(&bamstore_table);
      Uinttable_free(&bamstore_table);
    }
    for (i = 0; i < n; i++) {
      Chrom_free(&(chroms[i]));
    }
    FREE(chroms);
  }
  Table_free(&bamstore_chrtable);

  return;
}



static void
parse_gene (List_T *intronlist, List_T *middle_exonlist, List_T *end_exonlist,
	    int genei, char *gene, int sign,
	    Table_T fwd_splices_chrtable, Table_T rev_splices_chrtable,
	    char *chr, Genome_T genome, IIT_T chromosome_iit) {
  Genomicpos_T exonstart, exonend, firstpos, secondpos;
  char *p;
  Uintlist_T exonstarts = NULL, exonends = NULL, donors, acceptors, a, b;
  char truestrand;
  int nexons;


  p = gene;

  /* Skip gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) == 2) {
      exonstarts = Uintlist_push(exonstarts,exonstart);
      exonends = Uintlist_push(exonends,exonend);
    }

    /* Skip part just read */
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && *p != '\n' && isspace(*p)) { p++; }

    /* Read to end of line */
    if (*p == '\0' || *p == '\n') {
      /* info = (char *) CALLOC(1,sizeof(char)); */
      /* info[0] = '\0'; */
    } else {
      /* q = p; */
      while (*p != '\0' && *p != '\n') { p++; }

      /* info = (char *) CALLOC(p - q + 1,sizeof(char)); */
      /* strncpy(info,q,p-q); */
      if (*p == '\n') { p++; }
    }

    /* *info_list = List_push(*info_list,(void *) info); */
  }

  exonstarts = Uintlist_reverse(exonstarts);
  acceptors = Uintlist_next(exonstarts);

  exonends = Uintlist_reverse(exonends);
  donors = exonends;
  

  if (sign > 0) {
    truestrand = '+';
  } else {
    truestrand = '-';
  }

  /* Need to check b != NULL and not a != NULL */
  for (a = donors, b = acceptors; b != NULL; a = Uintlist_next(a), b = Uintlist_next(b)) {
    if (sign > 0) {
      firstpos = Uintlist_head(a);
      secondpos = Uintlist_head(b);
    } else {
      firstpos = Uintlist_head(b);
      secondpos = Uintlist_head(a);
    }

    add_splice(firstpos,secondpos,chr,/*crosshybp*/false,
	       /*nhits*/1,/*support1*/100,/*support2*/100,
	       /*overhang1*/0,/*overhang2*/0,/*querylength*/0,/*lowend_p*/true,
	       fwd_splices_chrtable,rev_splices_chrtable,
	       genome,chromosome_iit,truestrand,/*shortsplicedist*/1000000,
	       /*trust_sam_p*/true,/*need_canonical_p*/false,/*knownp*/true);
  }


  if ((nexons = Uintlist_length(exonstarts)) == 1) {
    exonstart = Uintlist_head(exonstarts);
    exonend = Uintlist_head(exonends);
    *end_exonlist = List_push(*end_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));

  } else if (nexons == 2) {
    exonstart = Uintlist_head(exonstarts);
    exonend = Uintlist_head(exonends);
    *end_exonlist = List_push(*end_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));

    exonstart = Uintlist_head(Uintlist_next(exonstarts));
    *intronlist = List_push(*intronlist,(void *) Interval_new(exonend,exonstart,/*type*/genei));

    exonend = Uintlist_head(Uintlist_next(exonends));
    *end_exonlist = List_push(*end_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));

  } else {
    exonstart = Uintlist_head(exonstarts);
    exonend = Uintlist_head(exonends);
    *end_exonlist = List_push(*end_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));

    for (a = Uintlist_next(exonstarts), b = Uintlist_next(exonends); Uintlist_next(a) != NULL; a = Uintlist_next(a), b = Uintlist_next(b)) {
      exonstart = Uintlist_head(a);
      *intronlist = List_push(*intronlist,(void *) Interval_new(exonend,exonstart,/*type*/genei));

      exonend = Uintlist_head(b);
      *middle_exonlist = List_push(*middle_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));
    }

    exonstart = Uintlist_head(a);
    *intronlist = List_push(*intronlist,(void *) Interval_new(exonend,exonstart,/*type*/genei));

    exonend = Uintlist_head(b);
    *end_exonlist = List_push(*end_exonlist,(void *) Interval_new(exonstart,exonend,/*type*/genei));
  }

  Uintlist_free(&exonstarts);
  Uintlist_free(&exonends);

  return;
}



T
Gstruct_knowngenes_input (IIT_T *introns_iit, IIT_T *middle_exons_iit, IIT_T *end_exons_iit,
			  IIT_T knowngenes_iit, Genome_T genome, IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));
  int divno;
  char *chr, *copy, Buffer[30];

  char *gene, *restofheader;
  Interval_T *array, interval;
  int *matches, nmatches, n, i;
  int sign;
  bool allocp;

  Table_T intron_intervaltable, intron_labeltable,
    middle_intervaltable, middle_labeltable, end_intervaltable, end_labeltable;
  List_T divlist, intronlist, middle_exonlist, end_exonlist,
    intron_labellist, middle_labellist, end_labellist, typelist, p, q;
  Intlist_T genei, r;
  char *divstring, *typestring, *intronlabel, *exonlabel, *genelabel;
  int introni = 0, exoni = 0, nchars;

  List_T unique, duplicates;


  new->fwd_extents_chrtable = (Table_T) NULL;
  new->rev_extents_chrtable = (Table_T) NULL;
  new->null_extents_chrtable = (Table_T) NULL;

  new->primary_extents_chrtable = (Table_T) NULL;
  new->crosshyb_extents_chrtable = (Table_T) NULL;

  new->fwd_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->rev_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  new->bamstore_chrtable = (Table_T) NULL;
  new->chrlength_table = (Tableuint_T) NULL;

  intron_intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  intron_labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  middle_intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  middle_labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  end_intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  end_labeltable = Table_new(65522,Table_string_compare,Table_string_hash);

  divlist = (List_T) NULL;
  for (divno = 1; divno < IIT_ndivs(knowngenes_iit); divno++) {
    chr = IIT_divstring(knowngenes_iit,divno);
    copy = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(copy,chr);
    divlist = List_push(divlist,(void *) copy);

    intronlist = middle_exonlist = end_exonlist = (List_T) NULL;
    intron_labellist = middle_labellist = end_labellist = (List_T) NULL;
    matches = IIT_get(&nmatches,knowngenes_iit,chr,/*x*/0,/*y*/-1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      gene = IIT_annotation(&restofheader,knowngenes_iit,matches[i],&allocp);
      interval = IIT_interval(knowngenes_iit,matches[i]);
      sign = Interval_sign(interval);
      parse_gene(&intronlist,&middle_exonlist,&end_exonlist,/*genei*/matches[i],gene,sign,
		 new->fwd_splices_chrtable,new->rev_splices_chrtable,
		 chr,genome,chromosome_iit);
      if (allocp) FREE(restofheader);
    }
    FREE(matches);


    /* Remove duplicates from introns */
    if ((n = List_length(intronlist)) == 0) {
      unique = (List_T) NULL;
      assert(intron_labellist == NULL);
    } else {
      array = (Interval_T *) List_to_array(intronlist,NULL);
      qsort(array,n,sizeof(Interval_T),Interval_cmp);
      List_free(&intronlist);

      unique = List_push(NULL,(void *) array[0]);
      genei = Intlist_push(NULL,Interval_type(array[0]));

      duplicates = (List_T) NULL;
      for (i = 1; i < n; i++) {
	if (Interval_cmp(&(array[i]),&(array[i-1])) == 0) {
	  duplicates = List_push(duplicates,(void *) array[i]);
	} else {
	  unique = List_push(unique,(void *) array[i]);
	  genei = Intlist_push(genei,Interval_type(array[i]));
	}
      }
      FREE(array);
    }
    for (p = duplicates; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&duplicates);

    intron_labellist = (List_T) NULL;
    for (r = genei; r != NULL; r = Intlist_next(r)) {
      sprintf(Buffer,"Intron%d",introni++);
      nchars = strlen(Buffer);

      intronlabel = (char *) CALLOC(nchars+1,sizeof(char));
      strcpy(intronlabel,Buffer);

      intron_labellist = List_push(intron_labellist,(void *) intronlabel);
    }
    Intlist_free(&genei);

    Table_put(intron_intervaltable,(void *) copy,unique);
    Table_put(intron_labeltable,(void *) copy,intron_labellist);



    /* Remove duplicates from middle exons */
    if ((n = List_length(middle_exonlist)) == 0) {
      unique = (List_T) NULL;
      assert(middle_labellist == NULL);
    } else {
      array = (Interval_T *) List_to_array(middle_exonlist,NULL);
      qsort(array,n,sizeof(Interval_T),Interval_cmp);
      List_free(&middle_exonlist);

      unique = List_push(NULL,(void *) array[0]);
      genei = Intlist_push(NULL,Interval_type(array[0]));

      duplicates = (List_T) NULL;
      for (i = 1; i < n; i++) {
	if (Interval_cmp(&(array[i]),&(array[i-1])) == 0) {
	  duplicates = List_push(duplicates,(void *) array[i]);
	} else {
	  unique = List_push(unique,(void *) array[i]);
	  genei = Intlist_push(genei,Interval_type(array[i]));
	}
      }
      FREE(array);
    }
    for (p = duplicates; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&duplicates);

    middle_labellist = (List_T) NULL;
    for (r = genei; r != NULL; r = Intlist_next(r)) {
      sprintf(Buffer,"Exon%d",exoni++);
      nchars = strlen(Buffer);

      exonlabel = (char *) CALLOC(nchars+1,sizeof(char));
      strcpy(exonlabel,Buffer);

      middle_labellist = List_push(middle_labellist,(void *) exonlabel);
    }
    Intlist_free(&genei);

    Table_put(middle_intervaltable,(void *) copy,unique);
    Table_put(middle_labeltable,(void *) copy,middle_labellist);


    /* Remove duplicates from end exons */
    if ((n = List_length(end_exonlist)) == 0) {
      unique = (List_T) NULL;
      assert(end_labellist == NULL);
    } else {
      array = (Interval_T *) List_to_array(end_exonlist,NULL);
      qsort(array,n,sizeof(Interval_T),Interval_cmp);
      List_free(&end_exonlist);

      unique = List_push(NULL,(void *) array[0]);
      genei = Intlist_push(NULL,Interval_type(array[0]));

      duplicates = (List_T) NULL;
      for (i = 1; i < n; i++) {
	if (Interval_cmp(&(array[i]),&(array[i-1])) == 0) {
	  duplicates = List_push(duplicates,(void *) array[i]);
	} else {
	  unique = List_push(unique,(void *) array[i]);
	  genei = Intlist_push(genei,Interval_type(array[i]));
	}
      }
      FREE(array);
    }
    for (p = duplicates; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&duplicates);

    end_labellist = (List_T) NULL;
    for (r = genei; r != NULL; r = Intlist_next(r)) {
      sprintf(Buffer,"Exon%d",exoni++);
      nchars = strlen(Buffer);

      exonlabel = (char *) CALLOC(nchars+1,sizeof(char));
      strcpy(exonlabel,Buffer);

      end_labellist = List_push(end_labellist,(void *) exonlabel);
    }
    Intlist_free(&genei);


    Table_put(end_intervaltable,(void *) copy,unique);
    Table_put(end_labeltable,(void *) copy,end_labellist);
  }

  new->fwd_chroms = (Chrom_T *) Table_keys(new->fwd_splices_chrtable,NULL);
  if ((new->fwd_nchroms = Table_length(new->fwd_splices_chrtable)) > 0) {
    qsort(new->fwd_chroms,new->fwd_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  new->rev_chroms = (Chrom_T *) Table_keys(new->rev_splices_chrtable,NULL);
  if ((new->rev_nchroms = Table_length(new->rev_splices_chrtable)) > 0) {
    qsort(new->rev_chroms,new->rev_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  new->fwd_chromi = new->rev_chromi = 0;



  /* Make genes_iit */
  divlist = List_reverse(divlist);

  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(divlist,(void *) divstring);

  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(NULL,(void *) typestring);

  *introns_iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,intron_intervaltable,
			    intron_labeltable,/*datatable*/NULL,/*divsort*/NO_SORT,
			    /*version*/IIT_LATEST_VERSION,/*presortedp*/false);

  *middle_exons_iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,middle_intervaltable,
				 middle_labeltable,/*datatable*/NULL,/*divsort*/NO_SORT,
				 /*version*/IIT_LATEST_VERSION,/*presortedp*/false);

  *end_exons_iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,end_intervaltable,
			      end_labeltable,/*datatable*/NULL,/*divsort*/NO_SORT,
			      /*version*/IIT_LATEST_VERSION,/*presortedp*/false);

  FREE(typestring);
  List_free(&typelist);

  for (q = divlist; q != NULL; q = List_next(q)) {
    divstring = (char *) List_head(q);

    intron_labellist = (List_T) Table_get(intron_labeltable,(void *) divstring);
    for (p = intron_labellist; p != NULL; p = List_next(p)) {
      intronlabel = (char *) List_head(p);
      FREE(intronlabel);
    }
    List_free(&intron_labellist);

    intronlist = (List_T) Table_get(intron_intervaltable,(void *) divstring);
    for (p = intronlist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&intronlist);


    middle_labellist = (List_T) Table_get(middle_labeltable,(void *) divstring);
    for (p = middle_labellist; p != NULL; p = List_next(p)) {
      exonlabel = (char *) List_head(p);
      FREE(exonlabel);
    }
    List_free(&middle_labellist);

    middle_exonlist = (List_T) Table_get(middle_intervaltable,(void *) divstring);
    for (p = middle_exonlist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&middle_exonlist);


    end_labellist = (List_T) Table_get(end_labeltable,(void *) divstring);
    for (p = end_labellist; p != NULL; p = List_next(p)) {
      exonlabel = (char *) List_head(p);
      FREE(exonlabel);
    }
    List_free(&end_labellist);

    end_exonlist = (List_T) Table_get(end_intervaltable,(void *) divstring);
    for (p = end_exonlist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&end_exonlist);

    FREE(divstring);
  }
  List_free(&divlist);

  Table_free(&end_labeltable);
  Table_free(&end_intervaltable);
  Table_free(&middle_labeltable);
  Table_free(&middle_intervaltable);
  Table_free(&intron_labeltable);
  Table_free(&intron_intervaltable);

  return new;
}



#else

static List_T
lines_gc (List_T *lines) {
  char *line;
  void *item;

  while (*lines != NULL) {
    *lines = List_pop(*lines,&item);
    line = (char *) item;
    FREE(line);
  }
  return NULL;
}


static void
parse_paired (char strand1, char strand2, char *chr1, char *chr2,
	      Genomicpos_T chrpos_low_1, Genomicpos_T chrpos_high_1,
	      Genomicpos_T chrpos_low_2, Genomicpos_T chrpos_high_2,
	      int nhits, Table_T fwd_extents_chrtable, Table_T rev_extents_chrtable,
	      Table_T null_extents_chrtable, int max_pairlength) {
  char strand;
  Genomicpos_T chrpos_low, chrpos_high;

  strand = combine_strands(strand1,strand2);

  debug(printf("strand1 is %c.  strand2 is %c.  strand is %c.\n",strand1,strand2,strand));

  debug(printf("Got %c%u..%u and %c%u..%u\n",
	       strand1,chrpos_low_1,chrpos_high_1,strand2,chrpos_low_2,chrpos_high_2));
  if (chrpos_low_1 < chrpos_low_2) {
    chrpos_low = chrpos_low_1;
  } else {
    chrpos_low = chrpos_low_2;
  }
  
  if (chrpos_high_1 > chrpos_high_2) {
    chrpos_high = chrpos_high_1;
  } else {
    chrpos_high = chrpos_high_2;
  }

  if (strand == '+') {
    Site_add_at_low(fwd_extents_chrtable,chr1,chrpos_low,chrpos_high,
		    max_pairlength,nhits);
  } else if (strand == '-') {
    Site_add_at_low(rev_extents_chrtable,chr1,chrpos_low,chrpos_high,
		    max_pairlength,nhits);
  } else if (strand == ' ') {
    Site_add_at_low(null_extents_chrtable,chr1,chrpos_low,chrpos_high,
		    max_pairlength,nhits);
    /* Also compute null_high and null_low */
  }

  return;
}



static void
parse_line (char strand, char *chr, Genomicpos_T chrpos_low, Genomicpos_T *chrpos_high,
	    char *cigar, Genome_T genome, IIT_T chromosome_iit, int nhits,
	    Table_T fwd_splices_chrtable, Table_T rev_splices_chrtable,
	    Genomicpos_T shortsplicedist, bool trust_sam_p, bool need_canonical_p) {
  Genomicpos_T chrpos, firstpos, secondpos;
  int cigar_readlength;
  int type;
  Intlist_T types, p;
  Uintlist_T npositions, q;

  bool splice_present_p = false;

  types = Samread_parse_cigar(&npositions,&cigar_readlength,cigar);

  chrpos = chrpos_low;
  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos += Uintlist_head(q);

    } else if (type == 'N') {
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;
      splice_present_p = true;

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      chrpos += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
    debug(printf("  type = %c, chrpos = %u\n",type,chrpos));
  }

  if (splice_present_p == true) {
    parse_splices(nhits,chr,chrpos_low,types,npositions,strand,
		  fwd_splices_chrtable,rev_splices_chrtable,
		  genome,chromosome_iit,shortsplicedist,
		  trust_sam_p,need_canonical_p);
  }


  Intlist_free(&types);
  Uintlist_free(&npositions);

  *chrpos_high = chrpos - 1U;

  /* FREE(*chr); */

  return;
}


static void
process_lines (Table_T fwd_extents_chrtable, Table_T rev_extents_chrtable,
	       Table_T fwd_splices_chrtable, Table_T rev_splices_chrtable,
	       Table_T null_extents_chrtable, List_T lines_first, List_T lines_second,
	       int nhits1, int nhits2, int minimum_mapq, bool need_primary_p,
	       Genomicpos_T max_pairlength, Genome_T genome, IIT_T chromosome_iit,
	       Genomicpos_T shortsplicedist, bool trust_sam_p, bool need_canonical_p) {
  List_T ptr1, ptr2;
  char *line1, *line2;
  int nhits;

  char *auxinfo;
  char strand, strand1, strand2;
  char *chr1, *chr2, *mate_chr;
  Genomicpos_T chrpos_low_1, chrpos_high_1, chrpos_low_2, chrpos_high_2, mate_chrpos_low;

  unsigned int flag1, flag2;
  char *acc, *read, *quality_string, *cigar1, *cigar2;
  int readlength;
  int mapq1, mapq2;


  debug0(printf("Entering process_lines with nhits %d\n",nhits));

  ptr1 = lines_first;
  ptr2 = lines_second;
  while (ptr1 != NULL && ptr2 != NULL) {
    line1 = (char *) List_head(ptr1);
    line2 = (char *) List_head(ptr2);
    if (nhits1 < nhits2) {
      nhits = nhits1;
    } else {
      nhits = nhits2;
    }

    auxinfo = Samread_parse_line(&acc,&flag1,&mapq1,&chr1,&chrpos_low_1,&cigar1,
				 &mate_chr,&mate_chrpos_low,&readlength,&read,&quality_string,line1);
    strand1 = Samread_splice_strand(auxinfo);
    FREE(acc);
    FREE(read);
    FREE(quality_string);
    FREE(mate_chr);

    auxinfo = Samread_parse_line(&acc,&flag2,&mapq2,&chr2,&chrpos_low_2,&cigar2,
				 &mate_chr,&mate_chrpos_low,&readlength,&read,&quality_string,line2);
    strand2 = Samread_splice_strand(auxinfo);
    FREE(acc);
    FREE(read);
    FREE(quality_string);
    FREE(mate_chr);

    strand = combine_strands(strand1,strand2);

    if (need_primary_p == true && ((flag1 & NOT_PRIMARY) || (flag2 & NOT_PRIMARY))) {
      /* Skip */
    } else if (mapq1 < minimum_mapq || mapq2 < minimum_mapq) {
      /* Skip */
    } else {
      parse_line(strand,chr1,chrpos_low_1,&chrpos_high_1,cigar1,genome,chromosome_iit,
		 nhits1,fwd_splices_chrtable,rev_splices_chrtable,shortsplicedist,
		 trust_sam_p,need_canonical_p);
      parse_line(strand,chr2,chrpos_low_2,&chrpos_high_2,cigar2,genome,chromosome_iit,
		 nhits2,fwd_splices_chrtable,rev_splices_chrtable,shortsplicedist,
		 trust_sam_p,need_canonical_p);
      parse_paired(strand1,strand2,chr1,chr2,chrpos_low_1,chrpos_high_1,
		   chrpos_low_2,chrpos_high_2,nhits,
		   fwd_extents_chrtable,rev_extents_chrtable,null_extents_chrtable,
		   max_pairlength);
    }

    FREE(cigar2);
    FREE(cigar1);

    FREE(chr2);
    FREE(chr1);

    ptr1 = List_next(ptr1);
    ptr2 = List_next(ptr2);
  }

  return;
}


T
Gstruct_sam_input (int minimum_mapq, Genomicpos_T max_pairlength,
		   Genome_T genome, IIT_T chromosome_iit, Genomicpos_T shortsplicedist,
		   bool need_unique_p, bool need_primary_p, bool trust_sam_p, bool need_canonical_p) {
  T new = (T) MALLOC(sizeof(*new));
  List_T lines_first = NULL, lines_second = NULL;
  char line[1024000], *copy;
  int nhits_first = 0, nhits_second = 0;
  char *lastacc, *acc;
  unsigned int flag;
  bool pairedp, concordantp;

  new->fwd_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->rev_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->null_extents_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  new->fwd_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  new->rev_splices_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);


  lastacc = (char *) CALLOC(1,sizeof(char));
  lastacc[0] = '\0';

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '@') {
      /* Skip */
    } else {
      acc = Samread_get_acc(&flag,line);
      if (strcmp(acc,lastacc)) {
	if (lastacc[0] != '\0') {
	  if (concordantp == true) {
	    if (need_unique_p == false || (nhits_first == 1 && nhits_second == 1)) {
	      process_lines(new->fwd_extents_chrtable,new->rev_extents_chrtable,
			    new->fwd_splices_chrtable,new->rev_splices_chrtable,
			    new->null_extents_chrtable,lines_first,lines_second,
			    nhits_first,nhits_second,minimum_mapq,need_primary_p,
			    max_pairlength,genome,chromosome_iit,
			    shortsplicedist,trust_sam_p,need_canonical_p);
	    }
	  }
	  lines_first = lines_gc(&lines_first);
	  lines_second = lines_gc(&lines_second);
	  nhits_first = 0;
	  nhits_second = 0;
	}
	FREE(lastacc);
	lastacc = acc;
      } else {
	FREE(acc);
      }

      if (flag & QUERY_UNMAPPED) {
	/* Skip line */
      } else {
	copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
	strcpy(copy,line);
	if (!(flag & PAIRED_READ)) {
	  pairedp = false;
	  concordantp = false;
	  lines_first = List_push(lines_first,(void *) copy);
	  nhits_first++;
	} else if (flag & FIRST_READ_P) {
	  pairedp = true;
	  concordantp = (flag & PAIRED_MAPPING) ? true : false;
	  lines_first = List_push(lines_first,(void *) copy);
	  nhits_first++;
	} else if (flag & SECOND_READ_P) {
	  pairedp = true;
	  concordantp = (flag & PAIRED_MAPPING) ? true : false;
	  lines_second = List_push(lines_second,(void *) copy);
	  nhits_second++;
	} else {
	  pairedp = true;
	  fprintf(stderr,"Flag %u is paired (%u), but contains neither first_read nor second_read flag\n",
		  flag,flag & PAIRED_MAPPING);
	  abort();
	}
      }
    }
  }

  if (lastacc[0] != '\0') {
    if (concordantp == true) {
      if (need_unique_p == false || (nhits_first == 1 && nhits_second == 1)) {
	process_lines(new->fwd_extents_chrtable,new->rev_extents_chrtable,
		      new->fwd_splices_chrtable,new->rev_splices_chrtable,
		      new->null_extents_chrtable,lines_first,lines_second,
		      nhits_first,nhits_second,minimum_mapq,need_primary_p,
		      max_pairlength,genome,chromosome_iit,
		      shortsplicedist,trust_sam_p,need_canonical_p);
      }
    }
    lines_first = lines_gc(&lines_first);
    lines_second = lines_gc(&lines_second);
  }
  FREE(lastacc);

  new->fwd_chroms = (Chrom_T *) Table_keys(new->fwd_splices_chrtable,NULL);
  if ((new->fwd_nchroms = Table_length(new->fwd_splices_chrtable)) > 0) {
    qsort(new->fwd_chroms,new->fwd_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  new->rev_chroms = (Chrom_T *) Table_keys(new->rev_splices_chrtable,NULL);
  if ((new->rev_nchroms = Table_length(new->rev_splices_chrtable)) > 0) {
    qsort(new->rev_chroms,new->rev_nchroms,sizeof(Chrom_T),Chrom_compare_chrom);
  }

  return new;
}

#endif


char *
Gstruct_next_chr (List_T *splices, long int **fwd_extents, long int **rev_extents, long int **null_extents,
		  long int **primary_extents, long int **crosshyb_extents,
		  Genomicpos_T *chrlength, T this,
		  int max_exonlength, int mincount_end_alt, int minsupport, bool need_canonical_p,
		  char *bam_lacks_chr, int bam_lacks_chr_length) {
  char *chr, *chrptr;
  Chrom_T chrom;
  int cmp;
  Splice_T *array;
  int nsplices;


  if (this->fwd_chromi >= this->fwd_nchroms && this->rev_chromi >= this->rev_nchroms) {
    *splices = (List_T) NULL;
    *fwd_extents = (long int *) NULL;
    *rev_extents = (long int *) NULL;
    *null_extents = (long int *) NULL;
    *primary_extents = (long int *) NULL;
    *crosshyb_extents = (long int *) NULL;
    *chrlength = 0;
    return (char *) NULL;

  } else if (this->rev_chromi >= this->rev_nchroms) {
    chrom = this->fwd_chroms[this->fwd_chromi++];
    chr = Chrom_string(chrom);
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }
    *chrlength = Tableuint_get(this->chrlength_table,(void *) chrom);

    *primary_extents = store_runlength((Uinttable_T) Table_get(this->primary_extents_chrtable,(void *) chrom),*chrlength);
    *crosshyb_extents = store_runlength((Uinttable_T) Table_get(this->crosshyb_extents_chrtable,(void *) chrom),*chrlength);

    *fwd_extents = store_runlength((Uinttable_T) Table_get(this->fwd_extents_chrtable,(void *) chrom),*chrlength);
    *rev_extents = (long int *) NULL;
    *null_extents = store_runlength((Uinttable_T) Table_get(this->null_extents_chrtable,(void *) chrom),*chrlength);
    *splices = accumulate_splices((List_T) NULL,(Uinttable_T) Table_get(this->fwd_splices_chrtable,(void *) chrom));
    Splice_filter_list(*splices,mincount_end_alt,minsupport,need_canonical_p);

    return chr;

  } else if (this->fwd_chromi >= this->fwd_nchroms) {
    chrom = this->rev_chroms[this->rev_chromi++];
    chr = Chrom_string(chrom);
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }
    *chrlength = Tableuint_get(this->chrlength_table,(void *) chrom);

    *primary_extents = store_runlength((Uinttable_T) Table_get(this->primary_extents_chrtable,(void *) chrom),*chrlength);
    *crosshyb_extents = store_runlength((Uinttable_T) Table_get(this->crosshyb_extents_chrtable,(void *) chrom),*chrlength);

    *fwd_extents = (long int *) NULL;
    *rev_extents = store_runlength((Uinttable_T) Table_get(this->rev_extents_chrtable,(void *) chrom),*chrlength);
    *null_extents = store_runlength((Uinttable_T) Table_get(this->null_extents_chrtable,(void *) chrom),*chrlength);
    *splices = accumulate_splices((List_T) NULL,(Uinttable_T) Table_get(this->rev_splices_chrtable,(void *) chrom));
    Splice_filter_list(*splices,mincount_end_alt,minsupport,need_canonical_p);

    return chr;

  } else if ((cmp = Chrom_compare_chrom(&(this->fwd_chroms[this->fwd_chromi]),&(this->rev_chroms[this->rev_chromi]))) < 0) {
    chrom = this->fwd_chroms[this->fwd_chromi++];
    chr = Chrom_string(chrom);
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }
    *chrlength = Tableuint_get(this->chrlength_table,(void *) chrom);

    *primary_extents = store_runlength((Uinttable_T) Table_get(this->primary_extents_chrtable,(void *) chrom),*chrlength);
    *crosshyb_extents = store_runlength((Uinttable_T) Table_get(this->crosshyb_extents_chrtable,(void *) chrom),*chrlength);

    *fwd_extents = store_runlength((Uinttable_T) Table_get(this->fwd_extents_chrtable,(void *) chrom),*chrlength);
    *rev_extents = (long int *) NULL;
    *null_extents = store_runlength((Uinttable_T) Table_get(this->null_extents_chrtable,(void *) chrom),*chrlength);
    *splices = accumulate_splices((List_T) NULL,(Uinttable_T) Table_get(this->fwd_splices_chrtable,(void *) chrom));
    Splice_filter_list(*splices,mincount_end_alt,minsupport,need_canonical_p);

    return chr;

  } else if (cmp > 0) {
    chrom = this->rev_chroms[this->rev_chromi++];
    chr = Chrom_string(chrom);
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }
    *chrlength = Tableuint_get(this->chrlength_table,(void *) chrom);

    *primary_extents = store_runlength((Uinttable_T) Table_get(this->primary_extents_chrtable,(void *) chrom),*chrlength);
    *crosshyb_extents = store_runlength((Uinttable_T) Table_get(this->crosshyb_extents_chrtable,(void *) chrom),*chrlength);

    *fwd_extents = (long int *) NULL;
    *rev_extents = store_runlength((Uinttable_T) Table_get(this->rev_extents_chrtable,(void *) chrom),*chrlength);
    *null_extents = store_runlength((Uinttable_T) Table_get(this->null_extents_chrtable,(void *) chrom),*chrlength);
    *splices = accumulate_splices((List_T) NULL,(Uinttable_T) Table_get(this->rev_splices_chrtable,(void *) chrom));
    Splice_filter_list(*splices,mincount_end_alt,minsupport,need_canonical_p);

    return chr;

  } else {
    chrom = this->fwd_chroms[this->fwd_chromi++];
    chrom = this->rev_chroms[this->rev_chromi++];
    chr = Chrom_string(chrom);
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }
    *chrlength = Tableuint_get(this->chrlength_table,(void *) chrom);

    *primary_extents = store_runlength((Uinttable_T) Table_get(this->primary_extents_chrtable,(void *) chrom),*chrlength);
    *crosshyb_extents = store_runlength((Uinttable_T) Table_get(this->crosshyb_extents_chrtable,(void *) chrom),*chrlength);

    *fwd_extents = store_runlength((Uinttable_T) Table_get(this->fwd_extents_chrtable,(void *) chrom),*chrlength);
    *rev_extents = store_runlength((Uinttable_T) Table_get(this->rev_extents_chrtable,(void *) chrom),*chrlength);
    *null_extents = store_runlength((Uinttable_T) Table_get(this->null_extents_chrtable,(void *) chrom),*chrlength);

    *splices = (List_T) NULL;
    *splices = accumulate_splices(*splices,(Uinttable_T) Table_get(this->fwd_splices_chrtable,(void *) chrom));
    *splices = accumulate_splices(*splices,(Uinttable_T) Table_get(this->rev_splices_chrtable,(void *) chrom));

    nsplices = List_length(*splices);
    array = (Splice_T *) List_to_array(*splices,NULL);

    /* Sets status to non-valid for various splices */
    if (*fwd_extents != NULL && *rev_extents != NULL) {
      Splice_resolve(array,nsplices,*fwd_extents,*rev_extents,*chrlength,max_exonlength,chr);
    }
#if 0
    /* Splice_resolve handles inconsistent directions using extents.  Splice_turn causes problems. */
    Splice_turn(array,nsplices,max_exonlength,chr);
#endif
    Splice_filter_array(array,nsplices,mincount_end_alt,minsupport,need_canonical_p);
    FREE(array);

    return chr;
  }
}



List_T
Gstruct_get_known (T this, char *chr) {
  List_T splices;
  Chrom_T chrom;

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);

  splices = accumulate_splices((List_T) NULL,(Uinttable_T) Table_get(this->fwd_splices_chrtable,(void *) chrom));
  splices = accumulate_splices(splices,(Uinttable_T) Table_get(this->rev_splices_chrtable,(void *) chrom));
  Chrom_free(&chrom);

  return splices;
}


