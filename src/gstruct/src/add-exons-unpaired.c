static char rcsid[] = "$Id: add-exons-unpaired.c 107255 2013-09-09 20:30:43Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>		/* For tolower, toupper, islower */
#include <math.h>		/* For qsort */

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "list.h"
#include "listdef.h"
#include "table.h"

#include "iit-read.h"
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define uabs(x,y) ((x < y) ? y - x : x - y)
#define ABSDIST_SUFFICIENT 10
#define FALSE_SPLICE 3

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;

static char *novelsites_iitfile = NULL;
static char *ignore_transcripts_file = NULL;
static int nflanking = 5;
static bool distantp = false;


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Known genes */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  {"novel", required_argument, 0, 'n'}, /* novelsites_iitfile */
  {"ignore", required_argument, 0, 'I'}, /* ignore_transcripts_file */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */
  {"distant", no_argument, 0, 0},	    /* distantp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"ADD-EXONS-UNPAIRED -- Finds readthrough splices in concordant GSNAP output\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: add-exons-unpaired [options] [-o <output>] <bam file...>\n\
\n\
Options:\n\
");

    return;
}


typedef struct Site_T *Site_T;
struct Site_T {
  int dist;			/* signed */
  Genomicpos_T pos;
  char *gene;
  char *acc;
  int exoni;
  int nexons;
  bool novelp;
  bool flankingp;
};

static Site_T
Site_new (int dist, Genomicpos_T pos, char *gene, int gene_namelength, char *acc,
	  int exoni, int nexons, bool novelp, bool flankingp) {
  Site_T new = (Site_T) MALLOC(sizeof(*new));
  char *p;
  int k;

  new->dist = dist;
  new->pos = pos;

  new->gene = (char *) CALLOC(gene_namelength + 1,sizeof(char));
  strncpy(new->gene,gene,gene_namelength);
  for (p = new->gene, k = 0; k < gene_namelength; p++, k++) {
    if (*p == '-') {
      *p = '.';
    }
  }

  new->acc = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
  strcpy(new->acc,acc);

  new->exoni = exoni;
  new->nexons = nexons;
  new->novelp = novelp;
  new->flankingp = flankingp;

  return new;
}

static void
Site_free (Site_T *old) {
  FREE((*old)->gene);
  FREE((*old)->acc);
  FREE(*old);
  return;
}


static int
string_cmp (const void *a, const void *b) {
  char *x = * (char **) a;
  char *y = * (char **) b;

  return strcmp(x,y);
}



#if 0
static int
Site_acc_cmp (const void *a, const void *b) {
  Site_T x = * (Site_T *) a;
  Site_T y = * (Site_T *) b;

  return strcmp(x->acc,y->acc);
}
#endif

static int
Site_pos_acc_cmp (const void *a, const void *b) {
  Site_T x = * (Site_T *) a;
  Site_T y = * (Site_T *) b;

  if (x->pos < y->pos) {
    return -1;
  } else if (y->pos < x->pos) {
    return +1;
  } else {
    return strcmp(x->acc,y->acc);
  }
}

#if 0
static List_T
Site_acc_sort (List_T sites) {
  List_T sorted = NULL;
  Site_T *array;
  int n, i;

  array = (Site_T *) List_to_array_n(&n,sites);
  List_free(&sites);

  qsort(array,n,sizeof(Site_T),Site_acc_cmp);
  for (i = n-1; i >= 0; i--) {
    sorted = List_push(sorted,(void *) array[i]);
  }
  FREE(array);

  return sorted;
}
#endif

static List_T
Site_pos_acc_sort (List_T sites) {
  List_T sorted = NULL;
  Site_T *array;
  int n, i;

  array = (Site_T *) List_to_array_n(&n,sites);

  qsort(array,n,sizeof(Site_T),Site_pos_acc_cmp);
  for (i = n-1; i >= 0; i--) {
    sorted = List_push(sorted,(void *) array[i]);
  }
  FREE(array);

  return sorted;
}

static void
Site_gc (List_T *sites) {
  Site_T site;
  List_T p;

  for (p = *sites; p != NULL; p = List_next(p)) {
    site = (Site_T) List_head(p);
    Site_free(&site);
  }
  List_free(&(*sites));
  return;
}

#if 0
static Genomicpos_T
Site_minimum_pos (int *dist, List_T sites) {
  Genomicpos_T minpos = -1U;
  Site_T site;
  List_T p;

  for (p = sites; p != NULL; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (site->pos < minpos) {
      *dist = site->dist;
      minpos = site->pos;
    }
  }

  return minpos;
}
#endif

#if 0
static Genomicpos_T
Site_maximum_pos (int *dist, List_T sites) {
  Genomicpos_T maxpos = 0;
  Site_T site;
  List_T p;

  for (p = sites; p != NULL; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (site->pos > maxpos) {
      *dist = site->dist;
      maxpos = site->pos;
    }
  }

  return maxpos;
}
#endif


static void
Site_print_genes (List_T start, List_T end) {
  Site_T site;
  char **array;
  int n, i;
  List_T printed = NULL, p, q;
  bool seenp;

  for (p = start; p != end; p = List_next(p)) {
    site = (Site_T) List_head(p);
    seenp = false;
    for (q = printed; q != NULL; q = List_next(q)) {
      if (!strcmp(site->gene,(char *) List_head(q))) {
	seenp = true;
      }
    }
    if (seenp == false) {
      printed = List_push(printed,(void *) site->gene);
    }
  }

  array = (char **) List_to_array_n(&n,printed);
  List_free(&printed);

  qsort(array,n,sizeof(char *),string_cmp);
  printf("%s",array[0]);
  for (i = 1; i < n; i++) {
    printf(",%s",array[i]);
  }
  FREE(array);

  return;
}

/* To be consistent with find-novel-splicesites, need to be comma-delimited */
static void
Site_print_donor_exons (List_T start, List_T end) {
  Site_T site;
  List_T p;

  for (p = start; p != end; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (p != start) {
      printf(",");
    }
    if (site->flankingp == true) {
      printf("*");
    }
    /* printf("%s_exondon%d/%d",site->acc,site->exoni,site->nexons); */
    if (site->exoni == 0 || site->nexons == 0) {
      printf("%s",site->acc);
    } else {
      printf("%s_exon%d/%d",site->acc,site->exoni,site->nexons);
    }
  }

  return;
}

static void
Site_print_acceptor_exons (List_T start, List_T end) {
  Site_T site;
  List_T p;

  for (p = start; p != end; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (p != start) {
      printf(",");
    }
    if (site->flankingp == true) {
      printf("*");
    }
    /* printf("%s_exonacc%d/%d",site->acc,site->exoni,site->nexons); */
    if (site->exoni == 0 || site->nexons == 0) {
      printf("%s",site->acc);
    } else {
      printf("%s_exon%d/%d",site->acc,site->exoni,site->nexons);
    }
  }

  return;
}



static int
get_exons (Uintlist_T *exonstarts, Uintlist_T *exonends, int *gene_namelength, char *annot) {
  Genomicpos_T exonstart, exonend;
  char *p;


  *exonstarts = (Uintlist_T) NULL;
  *exonends = (Uintlist_T) NULL;

  p = annot;
  /* Expecting gene to be first word in annotation */
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }
  *gene_namelength = p - annot;
  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      *exonstarts = Uintlist_push(*exonstarts,exonstart);
      *exonends = Uintlist_push(*exonends,exonend);
    }
    while (*p != '\0' && *p != '\n') {
      p++;
    }
    if (*p == '\n') p++;
  }

  *exonstarts = Uintlist_reverse(*exonstarts);
  *exonends = Uintlist_reverse(*exonends);

  if (exonstart <= exonend) {
    return +1;
  } else {
    return -1;
  }
}


static Site_T
get_best_donor_known (Genomicpos_T *min_absdist, char *acc, char *annot, Genomicpos_T chrpos,
		      bool flankingp) {
  Uintlist_T exonstarts, exonends, p, q;
  Genomicpos_T best_splicesite, exonstart, exonend, prev_exonend = 0U;
  int dist;
  int sign;
  int gene_namelength;
  int best_exoni, exoni = 0, nexons;

  sign = get_exons(&exonstarts,&exonends,&gene_namelength,annot);
  nexons = Uintlist_length(exonstarts);

  if (sign > 0) {
    for (p = exonstarts, q = exonends; p != NULL; p = Uintlist_next(p), q = Uintlist_next(q)) {
      exonstart = Uintlist_head(p);
      exonend = Uintlist_head(q);
      exoni++;

      if (exoni == 1) {
	*min_absdist = uabs(chrpos,exonend);
	best_exoni = exoni;
	best_splicesite = exonend;
      }

      if (Uintlist_next(p) == NULL) {
	/* We have a gene end, not exon end */
	if (chrpos >= exonstart && chrpos <= exonend &&
	    chrpos <= exonstart + FALSE_SPLICE && prev_exonend != 0U) {
	  /* False alignment into next exon */
	  *min_absdist = 0;
	  best_exoni = exoni - 1;
	  best_splicesite = prev_exonend;
	}

      } else if (chrpos >= exonstart && chrpos <= exonend) {
	*min_absdist = 0;
	if (chrpos <= exonstart + FALSE_SPLICE && prev_exonend != 0U) {
	  /* False alignment into next exon */
	  best_exoni = exoni - 1;
	  best_splicesite = prev_exonend;
	} else {
	  best_exoni = exoni;
	  best_splicesite = exonend;
	}

      } else if (uabs(chrpos,exonend) < *min_absdist) {
	*min_absdist = uabs(chrpos,exonend);
	best_exoni = exoni;
	best_splicesite = exonend;
      }

      prev_exonend = exonend;
    }
      
  } else {
    for (p = exonstarts, q = exonends; p != NULL; p = Uintlist_next(p), q = Uintlist_next(q)) {
      exonstart = Uintlist_head(p);
      exonend = Uintlist_head(q);
      exoni++;

      if (exoni == 1) {
	*min_absdist = uabs(exonend,chrpos);
	best_exoni = exoni;
	best_splicesite = exonend;
      }

      if (Uintlist_next(p) == NULL) {
	/* We have a gene end, not exon end */
	if (chrpos >= exonend && chrpos <= exonstart &&
	    exonstart <= chrpos + FALSE_SPLICE && prev_exonend != 0U) {
	  /* False alignment into next exon */
	  *min_absdist = 0;
	  best_exoni = exoni - 1;
	  best_splicesite = prev_exonend;
	}

      } else if (chrpos >= exonend && chrpos <= exonstart) {
	*min_absdist = 0;
	if (exonstart <= chrpos + FALSE_SPLICE && prev_exonend != 0U) {
	  /* False alignment into next exon */
	  best_exoni = exoni - 1;
	  best_splicesite = prev_exonend;
	} else {
	  best_exoni = exoni;
	  best_splicesite = exonend;
	}

      } else if (uabs(exonend,chrpos) < *min_absdist) {
	*min_absdist = uabs(exonend,chrpos);
	best_exoni = exoni;
	best_splicesite = exonend;
      }

      prev_exonend = exonend;
    }
  }

  Uintlist_free(&exonstarts);
  Uintlist_free(&exonends);

  if (exoni == 1) {
    return NULL;
  } else {
    if (*min_absdist == 0) {
      dist = 0;
    } else if (sign > 0) {
      dist = chrpos - best_splicesite;
    } else {
      dist = best_splicesite - chrpos;
    }
    if (*min_absdist <= ABSDIST_SUFFICIENT) {
      *min_absdist = 0;
    }

    debug(printf("Best donor is acc %s, position %u, exon %d out of %d, with dist %d\n",
		 acc,best_splicesite,best_exoni,nexons,dist));
    return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		    best_exoni,nexons,/*novelp*/false,flankingp);
  }
}


static Site_T
get_best_acceptor_known (Genomicpos_T *min_absdist, char *acc, char *annot, Genomicpos_T chrpos,
			 bool flankingp) {
  Uintlist_T exonstarts, exonends, p, q;
  Genomicpos_T best_splicesite, exonstart, exonend, prev_exonstart = 0U;
  int dist;
  int sign;
  int gene_namelength;
  int best_exoni, exoni, nexons;

  sign = get_exons(&exonstarts,&exonends,&gene_namelength,annot);

  /* Want to progress from distal to medial */
  exonstarts = Uintlist_reverse(exonstarts);
  exonends = Uintlist_reverse(exonends);

  nexons = Uintlist_length(exonstarts);
  exoni = nexons + 1;
  if (sign > 0) {
    for (p = exonstarts, q = exonends; p != NULL; p = Uintlist_next(p), q = Uintlist_next(q)) {
      exonstart = Uintlist_head(p);
      exonend = Uintlist_head(q);
      exoni--;

      if (exoni == nexons) {
	*min_absdist = uabs(chrpos,exonstart);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }

      if (Uintlist_next(p) == NULL) {
	/* We have a gene start, not exon start */
	if (chrpos >= exonstart && chrpos <= exonend &&
	    exonend <= chrpos + FALSE_SPLICE && prev_exonstart != 0U) {
	  /* False alignment into next exon */
	  *min_absdist = 0;
	  best_exoni = exoni + 1;
	  best_splicesite = prev_exonstart;
	}

      } else if (chrpos >= exonstart && chrpos <= exonend) {
	*min_absdist = 0;
	if (exonend <= chrpos + FALSE_SPLICE && prev_exonstart != 0U) {
	  /* False alignment into next exon */
	  best_exoni = exoni + 1;
	  best_splicesite = prev_exonstart;
	} else {
	  best_exoni = exoni;
	  best_splicesite = exonstart;
	}

      } else if (uabs(chrpos, exonstart) < *min_absdist) {
	*min_absdist = uabs(chrpos,exonstart);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }

      prev_exonstart = exonstart;
    }

  } else {
    for (p = exonstarts, q = exonends; p != NULL; p = Uintlist_next(p), q = Uintlist_next(q)) {
      exonstart = Uintlist_head(p);
      exonend = Uintlist_head(q);
      exoni--;

      if (exoni == nexons) {
	*min_absdist = uabs(exonstart,chrpos);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }

      if (Uintlist_next(p) == NULL) {
	/* We have a gene start, not exon start */
	if (chrpos >= exonend && chrpos <= exonstart &&
	    chrpos <= exonend + FALSE_SPLICE && prev_exonstart != 0U) {
	  /* False alignment into next exon */
	  *min_absdist = 0;
	  best_exoni = exoni + 1;
	  best_splicesite = prev_exonstart;
	}
	  
      } else if (chrpos >= exonend && chrpos <= exonstart) {
	*min_absdist = 0;
	if (chrpos <= exonend + FALSE_SPLICE && prev_exonstart != 0U) {
	  /* False alignment into next exon */
	  best_exoni = exoni + 1;
	  best_splicesite = prev_exonstart;
	} else {
	  best_exoni = exoni;
	  best_splicesite = exonstart;
	}
	
      } else if (uabs(exonstart,chrpos) < *min_absdist) {
	*min_absdist = uabs(exonstart,chrpos);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }

      prev_exonstart = exonstart;
    }
  }

  Uintlist_free(&exonstarts);
  Uintlist_free(&exonends);

  if (exoni == nexons) {
    return NULL;
  } else {
    if (*min_absdist == 0) {
      dist = 0;
    } else if (sign > 0) {
      dist = chrpos - best_splicesite;
    } else {
      dist = best_splicesite - chrpos;
    }
    if (*min_absdist <= ABSDIST_SUFFICIENT) {
      *min_absdist = 0;
    }

    debug(printf("Best acceptor is acc %s, position %u, exon %d out of %d, with dist %d\n",
		 acc,best_splicesite,best_exoni,nexons,dist));
    return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		    best_exoni,nexons,/*novelp*/false,flankingp);
  }
}



static Site_T
get_donor_novel (Genomicpos_T *min_absdist, char *acc, char *annot, Interval_T interval, Genomicpos_T chrpos,
		 bool flankingp) {
  Genomicpos_T best_splicesite;
  Genomicpos_T splicestart, spliceend;
  int dist;
  int sign;
  char *p;
  int gene_namelength;

  p = annot;
  /* Expecting gene to be first word in annotation */
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }
  gene_namelength = p - annot;

  if ((sign = Interval_sign(interval)) > 0) {
    splicestart = Interval_low(interval);
    spliceend = Interval_high(interval);
  } else {
    splicestart = Interval_high(interval);
    spliceend = Interval_low(interval);
  }

  if (sign > 0) {
    if (chrpos >= splicestart && chrpos <= spliceend) {
      *min_absdist = 0;
      best_splicesite = spliceend;
    } else {
      *min_absdist = uabs(chrpos,spliceend);
      best_splicesite = spliceend;
    }
    
  } else {
    if (chrpos >= spliceend && chrpos <= splicestart) {
      *min_absdist = 0;
      best_splicesite = spliceend;
    } else {
      *min_absdist = uabs(spliceend,chrpos);
      best_splicesite = spliceend;
    }
  }

  if (*min_absdist == 0) {
    dist = 0;
  } else if (sign > 0) {
    dist = chrpos - best_splicesite;
  } else {
    dist = best_splicesite - chrpos;
  }
  if (*min_absdist <= ABSDIST_SUFFICIENT) {
    *min_absdist = 0;
  }

  return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		  /*best_exoni*/0,/*nexons*/0,/*novelp*/true,flankingp);
}


static Site_T
get_acceptor_novel (Genomicpos_T *min_absdist, char *acc, char *annot, Interval_T interval, Genomicpos_T chrpos,
		    bool flankingp) {
  Genomicpos_T best_splicesite;
  Genomicpos_T splicestart, spliceend;
  int dist;
  int sign;
  char *p;
  int gene_namelength;

  p = annot;
  /* Expecting gene to be first word in annotation */
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }
  gene_namelength = p - annot;

  if ((sign = Interval_sign(interval)) > 0) {
    splicestart = Interval_low(interval);
    spliceend = Interval_high(interval);
  } else {
    splicestart = Interval_high(interval);
    spliceend = Interval_low(interval);
  }

  if (sign > 0) {
    if (chrpos >= splicestart && chrpos <= spliceend) {
      *min_absdist = 0;
      best_splicesite = splicestart;
    } else {
      *min_absdist = uabs(chrpos,splicestart);
      best_splicesite = splicestart;
    }
    
  } else {
    if (chrpos >= spliceend && chrpos <= splicestart) {
      *min_absdist = 0;
      best_splicesite = splicestart;
    } else {
      *min_absdist = uabs(splicestart,chrpos);
      best_splicesite = splicestart;
    }
  }

  if (*min_absdist == 0) {
    dist = 0;
  } else if (sign > 0) {
    dist = chrpos - best_splicesite;
  } else {
    dist = best_splicesite - chrpos;
  }
  if (*min_absdist <= ABSDIST_SUFFICIENT) {
    *min_absdist = 0;
  }

  return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		  /*best_exoni*/0,/*nexons*/0,/*novelp*/true,flankingp);
}



static List_T
get_donor_sites (int *nbadchrs, Genomicpos_T *min_absdist, char *divstring,
		 char genestrand, Genomicpos_T outside, Genomicpos_T inside, IIT_T map_iit,
		 Table_T ignore_transcripts_table, IIT_T novelsites_iit, int donor_typeint) {
  List_T donor_sites = NULL;
  int nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int *matches;
  Genomicpos_T coordstart, coordend;
  Genomicpos_T absdist;
  Interval_T interval;
  Site_T site;
  int sign;
  char *acc, *annot, *restofheader;
  bool allocp;
  int i;

  if (outside < inside) {
    coordstart = outside;
    coordend = inside;
  } else {
    coordstart = inside;
    coordend = outside;
  }
  debug(printf("Entering get_donor_sites with %c%s:%u..%u\n",genestrand,divstring,coordstart,coordend));

  if (genestrand == '-') {
    sign = -1;
  } else {
    sign = +1;
  }

  donor_sites = (List_T) NULL;

  if (map_iit != NULL) {
    if (IIT_divint(map_iit,divstring) < 0) {
      (*nbadchrs)++;
    }

    /* Okay for this to be signed (but not add-exons-transloc) because we are looking only for known sites */
    matches = IIT_get_signed(&nmatches,map_iit,divstring,coordstart,coordend,sign,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      acc = IIT_label(map_iit,matches[i],&allocp);
      annot = IIT_annotation(&restofheader,map_iit,matches[i],&allocp);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_best_donor_known(&absdist,acc,annot,inside,/*flankingp*/false)) == NULL) {
	/* Skip */
      } else if (donor_sites == NULL) {
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&donor_sites);
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	donor_sites = List_push(donor_sites,(void *) site);
      } else {
	Site_free(&site);
      }

      if (allocp) {
	FREE(restofheader);
      }
    }
    FREE(matches);

    if (donor_sites == NULL && nflanking > 0) {
      debug(printf("Getting flanking\n"));
      IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,divstring,
		       coordstart,coordend,nflanking,sign);
      for (i = 0; i < nleftflanks; i++) {
	acc = IIT_label(map_iit,leftflanks[i],&allocp);
	annot = IIT_annotation(&restofheader,map_iit,leftflanks[i],&allocp);

	if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	  /* Skip */
	} else if ((site = get_best_donor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	  /* Skip */
	} else if (donor_sites == NULL) {
	  donor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist < *min_absdist) {
	  Site_gc(&donor_sites);
	  donor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist == *min_absdist) {
	  donor_sites = List_push(donor_sites,(void *) site);
	} else {
	  Site_free(&site);
	}

	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(leftflanks);

      for (i = 0; i < nrightflanks; i++) {
	acc = IIT_label(map_iit,rightflanks[i],&allocp);
	annot = IIT_annotation(&restofheader,map_iit,rightflanks[i],&allocp);

	if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	  /* Skip */
	} else if ((site = get_best_donor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	  /* Skip */
	} else if (donor_sites == NULL) {
	  donor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist < *min_absdist) {
	  Site_gc(&donor_sites);
	  donor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist == *min_absdist) {
	  donor_sites = List_push(donor_sites,(void *) site);
	} else {
	  Site_free(&site);
	}

	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(rightflanks);
    }
  }

  if (novelsites_iit != NULL) {
    matches = IIT_get_typed(&nmatches,novelsites_iit,divstring,coordstart,coordend,donor_typeint,/*sortp*/false);
    IIT_get_flanking_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,novelsites_iit,divstring,
			   coordstart,coordend,/*nflanking*/1,donor_typeint,/*sign*/0);
    for (i = 0; i < nmatches; i++) {
      acc = IIT_label(novelsites_iit,matches[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,matches[i],&allocp);
      interval = IIT_interval(novelsites_iit,matches[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_donor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/false)) == NULL) {
	/* Skip */
      } else if (donor_sites == NULL) {
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&donor_sites);
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	donor_sites = List_push(donor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(matches);

    for (i = 0; i < nleftflanks; i++) {
      acc = IIT_label(novelsites_iit,leftflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,leftflanks[i],&allocp);
      interval = IIT_interval(novelsites_iit,leftflanks[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_donor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (donor_sites == NULL) {
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&donor_sites);
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	donor_sites = List_push(donor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(leftflanks);

    for (i = 0; i < nrightflanks; i++) {
      acc = IIT_label(novelsites_iit,rightflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,rightflanks[i],&allocp);
      interval = IIT_interval(novelsites_iit,rightflanks[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_donor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (donor_sites == NULL) {
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&donor_sites);
	donor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	donor_sites = List_push(donor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(rightflanks);

  }

  return donor_sites;
}

static List_T
get_acceptor_sites (int *nbadchrs, Genomicpos_T *min_absdist, char *divstring,
		    char genestrand, Genomicpos_T outside, Genomicpos_T inside, IIT_T map_iit,
		    Table_T ignore_transcripts_table,
		    IIT_T novelsites_iit, int acceptor_typeint) {
  List_T acceptor_sites = NULL;
  int nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int *matches;
  Genomicpos_T coordstart, coordend;
  Genomicpos_T absdist;
  Interval_T interval;
  Site_T site;
  int sign;
  char *acc, *annot, *restofheader;
  bool allocp;
  int i;

  if (outside < inside) {
    coordstart = outside;
    coordend = inside;
  } else {
    coordstart = inside;
    coordend = outside;
  }
  debug(printf("Entering get_acceptor_sites with %c%s:%u..%u\n",genestrand,divstring,coordstart,coordend));

  if (genestrand == '-') {
    sign = -1;
  } else {
    sign = +1;
  }

  acceptor_sites = (List_T) NULL;

  if (map_iit != NULL) {
    if (IIT_divint(map_iit,divstring) < 0) {
      (*nbadchrs)++;
    }

    /* Okay for this to be signed (but not add-exons-transloc) because we are looking only for known sites */
    matches = IIT_get_signed(&nmatches,map_iit,divstring,coordstart,coordend,sign,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      acc = IIT_label(map_iit,matches[i],&allocp);
      annot = IIT_annotation(&restofheader,map_iit,matches[i],&allocp);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_best_acceptor_known(&absdist,acc,annot,inside,/*flankingp*/false)) == NULL) {
	/* Skip */
      } else if (acceptor_sites == NULL) {
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&acceptor_sites);
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	acceptor_sites = List_push(acceptor_sites,(void *) site);
      } else {
	Site_free(&site);
      }

      if (allocp) {
	FREE(restofheader);
      }
    }
    FREE(matches);

    if (acceptor_sites == NULL && nflanking > 0) {
      debug(printf("Getting flanking\n"));
      IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,divstring,
		       coordstart,coordend,nflanking,sign);
      for (i = 0; i < nleftflanks; i++) {
	acc = IIT_label(map_iit,leftflanks[i],&allocp);
	annot = IIT_annotation(&restofheader,map_iit,leftflanks[i],&allocp);

	if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	  /* Skip */
	} else if ((site = get_best_acceptor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	  /* Skip */
	} else if (acceptor_sites == NULL) {
	  acceptor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist < *min_absdist) {
	  Site_gc(&acceptor_sites);
	  acceptor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist == *min_absdist) {
	  acceptor_sites = List_push(acceptor_sites,(void *) site);
	} else {
	  Site_free(&site);
	}

	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(leftflanks);

      for (i = 0; i < nrightflanks; i++) {
	acc = IIT_label(map_iit,rightflanks[i],&allocp);
	annot = IIT_annotation(&restofheader,map_iit,rightflanks[i],&allocp);

	if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	  /* Skip */
	} else if ((site = get_best_acceptor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	  /* Skip */
	} else if (acceptor_sites == NULL) {
	  acceptor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist < *min_absdist) {
	  Site_gc(&acceptor_sites);
	  acceptor_sites = List_push(NULL,(void *) site);
	  *min_absdist = absdist;
	} else if (absdist == *min_absdist) {
	  acceptor_sites = List_push(acceptor_sites,(void *) site);
	} else {
	  Site_free(&site);
	}

	if (allocp) {
	  FREE(restofheader);
	}
      }
      FREE(rightflanks);
    }
  }

  if (novelsites_iit != NULL) {
    matches = IIT_get_typed(&nmatches,novelsites_iit,divstring,coordstart,coordend,acceptor_typeint,/*sortp*/false);
    IIT_get_flanking_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,novelsites_iit,divstring,
			   coordstart,coordend,/*nflanking*/1,acceptor_typeint,/*sign*/0);
    for (i = 0; i < nmatches; i++) {
      acc = IIT_label(novelsites_iit,matches[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,matches[i],&allocp);
      interval = IIT_interval(novelsites_iit,matches[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_acceptor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/false)) == NULL) {
	/* Skip */
      } else if (acceptor_sites == NULL) {
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&acceptor_sites);
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	acceptor_sites = List_push(acceptor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(matches);

    for (i = 0; i < nleftflanks; i++) {
      acc = IIT_label(novelsites_iit,leftflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,leftflanks[i],&allocp);
      interval = IIT_interval(novelsites_iit,leftflanks[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_acceptor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (acceptor_sites == NULL) {
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&acceptor_sites);
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	acceptor_sites = List_push(acceptor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(leftflanks);

    for (i = 0; i < nrightflanks; i++) {
      acc = IIT_label(novelsites_iit,rightflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,novelsites_iit,rightflanks[i],&allocp);
      interval = IIT_interval(novelsites_iit,rightflanks[i]);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_acceptor_novel(&absdist,acc,annot,interval,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (acceptor_sites == NULL) {
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist < *min_absdist) {
	Site_gc(&acceptor_sites);
	acceptor_sites = List_push(NULL,(void *) site);
	*min_absdist = absdist;
      } else if (absdist == *min_absdist) {
	acceptor_sites = List_push(acceptor_sites,(void *) site);
      } else {
	Site_free(&site);
      }
    }
    FREE(rightflanks);
  }

  return acceptor_sites;
}


static bool
sites_close_p (List_T donor_sites, int min_donor_absdist,
	       List_T acceptor_sites, int min_acceptor_absdist) {
  if (donor_sites != NULL && acceptor_sites != NULL &&
      min_donor_absdist <= 10 && min_acceptor_absdist <= 10) {
    return true;
  } else {
    return false;
  }
}


static void
print_sites (List_T donor_sites, char *donor_chr, char donor_genestrand,
	     Genomicpos_T donor_outside, Genomicpos_T donor_inside,
	     List_T acceptor_sites, char *acceptor_chr, char acceptor_genestrand,
	     Genomicpos_T acceptor_outside, Genomicpos_T acceptor_inside) {
  List_T p, q, r, s;
  Genomicpos_T donorpos, acceptorpos;
  int dist;

  donor_sites = Site_pos_acc_sort(donor_sites);
  acceptor_sites = Site_pos_acc_sort(acceptor_sites);

  p = donor_sites;
  while (p != NULL) {
    donorpos = ((Site_T) p->first)->pos;
    q = p->rest;
    while (q != NULL && ((Site_T) q->first)->pos == donorpos) {
      q = q->rest;
    }
    
    r = acceptor_sites;
    while (r != NULL) {
      acceptorpos = ((Site_T) r->first)->pos;
      s = r->rest;
      while (s != NULL && ((Site_T) s->first)->pos == acceptorpos) {
	s = s->rest;
      }

      printf("%s\t%c\t%u\t%u\t",donor_chr,donor_genestrand,donor_outside,donor_inside);
      printf("%u\t",donorpos);
      if ((dist = ((Site_T) p->first)->dist) == 0) {
	printf("0\t");
      } else if (dist > 0) {
	printf("+%d\t",dist);
      } else {
	printf("%d\t",dist);
      }

      printf("%s\t%c\t%u\t%u\t",acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside);
      printf("%u\t",acceptorpos);
      if ((dist = ((Site_T) r->first)->dist) == 0) {
	printf("0\t");
      } else if (dist > 0) {
	printf("+%d\t",dist);
      } else {
	printf("%d\t",dist);
      }

      Site_print_genes(p,q);
      printf("\t");
      Site_print_genes(r,s);
      printf("\t");

      Site_print_donor_exons(p,q);
      printf("\t");
      Site_print_acceptor_exons(r,s);
      printf("\n");
      
      r = s;
    }
    
    p = q;
  }

  List_free(&donor_sites);
  List_free(&acceptor_sites);

  return;
}



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL, *mapdir = NULL;
  IIT_T map_iit = NULL, novelsites_iit = NULL;
  int donor_typeint, acceptor_typeint;
  Table_T ignore_transcripts_table = NULL;
  char **keys;
  int nkeys, i;

  char *iitfile;

  FILE *fp;
  char Buffer[1024], acc[1024], chr1[1000], chr2[1000], cigar1[1000], cigar2[1000], *copy;
  char genestrand1, genestrand2;
  char readstrand1, readstrand2;
  Genomicpos_T outside1, inside1, outside2, inside2;
  Genomicpos_T min_donor_absdist_1, min_acceptor_absdist_1, min_donor_absdist_2, min_acceptor_absdist_2;
  List_T donor_sites_1, acceptor_sites_1, donor_sites_2, acceptor_sites_2;
  int nbadchrs = 0;

  bool case1p, case2p, case3p, case4p, case5p, case6p, case7p, case8p;


  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:M:m:n:I:u:",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"distant")) {
	distantp = true;
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'bam_fasta --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;

    case 'n': novelsites_iitfile = optarg; break;
    case 'I': ignore_transcripts_file = optarg; break;
    case 'u': nflanking = atoi(optarg); break;

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      

  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag\n");
    print_program_usage();
    exit(9);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    FREE(dbversion);
    FREE(dbroot);
  }

  mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
  if (map_iitfile == NULL) {
    /* Skip */
  } else if ((map_iit = IIT_read(map_iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
				 /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+strlen(map_iitfile)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,map_iitfile);
    if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
			    /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
    }
    FREE(iitfile);
  }
  FREE(mapdir);
  FREE(fileroot);
  FREE(genomesubdir);

  if (ignore_transcripts_file != NULL) {
    ignore_transcripts_table = Table_new(100,Table_string_compare,Table_string_hash);
    fp = fopen(ignore_transcripts_file,"r");
    while (fgets(Buffer,1024,fp) != NULL) {
      if (sscanf(Buffer,"%s",acc) == 1) {
	copy = (char *) CALLOC(strlen(acc)+1,sizeof(char));
	strcpy(copy,acc);
	Table_put(ignore_transcripts_table,copy,/*value*/(void *) 1);
      }
    }
    fclose(fp);
  }

  if (novelsites_iitfile != NULL) {
    novelsites_iit = IIT_read(novelsites_iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true);
    donor_typeint = IIT_typeint(novelsites_iit,"donor");
    acceptor_typeint = IIT_typeint(novelsites_iit,"acceptor");
  }

  while (fgets(Buffer,1024,stdin) != NULL) {
    if (sscanf(Buffer,"%s %c %u %u %c %s %s %c %u %u %c %s",
	       chr1,&readstrand1,&outside1,&inside1,&genestrand1,cigar1,chr2,&readstrand2,&outside2,&inside2,&genestrand2,cigar2) == 12) {
      debug(printf("Read %s\n",Buffer));

      if (genestrand1 == '.' && genestrand2 == '.') {
	case1p = case2p = true;
	if (distantp == false) {
	  case3p = case4p = case5p = case6p = case7p = case8p = false;
	} else {
	  case3p = case4p = case5p = case6p = case7p = case8p = true;
	}
      } else if (genestrand1 == '.') {
	if (genestrand2 == '+') {
	  case1p = true;
	  case2p = false;
	  if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	  } else {
	    case4p = case6p = case7p = true;
	    case3p = case5p = case8p = false;
	  }
	} else if (genestrand2 == '-') {
	  case1p = false;
	  case2p = true;
	  if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	  } else {
	    case3p = case5p = case8p = true;
	    case4p = case6p = case7p = false;
	  }
	} else {
	  abort();
	}
      } else if (genestrand2 == '.') {
	if (genestrand1 == '+') {
	  case1p = true;
	  case2p = false;
	  if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	  } else {
	    case4p = case5p = case8p = true;
	    case3p = case6p = case7p = false;
	  }
	} else if (genestrand1 == '-') {
	  case1p = false;
	  case2p = true;
	  if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	  } else {
	    case3p = case6p = case7p = true;
	    case4p = case5p = case8p = false;
	  }
	}
      } else if (genestrand1 == '+' && genestrand2 == '+') {
	case1p = true;
	case2p = false;
	if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	} else {
	  case4p = true;
	  case3p = case5p = case6p = case7p = case8p = false;
	}
      } else if (genestrand1 == '-' && genestrand2 == '-') {
	case1p = false;
	case2p = true;
	if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	} else {
	  case3p = true;
	  case4p = case5p = case6p = case7p = case8p = false;
	}
      } else if (genestrand1 == '+' && genestrand2 == '-') {
	case1p = case2p = false;
	if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	    /* No need for "\n", because Buffer already has it */
	    fprintf(stderr,"Inconsistent gene strands for local pair: %s",Buffer);
	} else {
	  case5p = case8p = true;
	  case3p = case4p = case6p = case7p = false;
	}
      } else if (genestrand1 == '-' && genestrand2 == '+') {
	case1p = case2p = false;
	if (distantp == false) {
	    case3p = case4p = case5p = case6p = case7p = case8p = false;
	    /* No need for "\n", because Buffer already has it */
	    fprintf(stderr,"Inconsistent gene strands for local pair: %s",Buffer);
	} else {
	  case6p = case7p = true;
	  case3p = case4p = case5p = case8p = false;
	}
      } else {
	abort();
      }

      if (case1p == true) {
	/* Case 1: Gene on plus strand: Donor plus, Acceptor plus */
	donor_sites_1 = get_donor_sites(&nbadchrs,&min_donor_absdist_1,chr1,/*genestrand*/'+',
					outside1,inside1,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_2 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_2,chr2,/*genestrand*/'+',
					      outside2,inside2,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_1,min_donor_absdist_1,acceptor_sites_2,min_acceptor_absdist_2) == true) {
	  print_sites(donor_sites_1,chr1,/*donor_genestrand*/'+',outside1,inside1,
		      acceptor_sites_2,chr2,/*acceptor_genstrand*/'+',outside2,inside2);
	}
	Site_gc(&acceptor_sites_2);
	Site_gc(&donor_sites_1);
      }

      if (case2p == true) {
	/* Case 2: Gene on minus strand: Acceptor minus, Donor minus */
	donor_sites_2 = get_donor_sites(&nbadchrs,&min_donor_absdist_2,chr2,/*genestrand*/'-',
					outside2,inside2,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_1 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_1,chr1,/*genestrand*/'-',
					      outside1,inside1,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_2,min_donor_absdist_2,acceptor_sites_1,min_acceptor_absdist_1) == true) {
	  print_sites(donor_sites_2,chr2,/*donor_genestrand*/'-',outside2,inside2,
		      acceptor_sites_1,chr1,/*acceptor_genestrand*/'-',outside1,inside1);
	}
	Site_gc(&acceptor_sites_1);
	Site_gc(&donor_sites_2);
      }

      if (case3p == true) {
	/* Case 3: Donor minus, Acceptor minus */
	donor_sites_1 = get_donor_sites(&nbadchrs,&min_donor_absdist_1,chr1,/*genestrand*/'-',
					outside1,inside1,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_2 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_2,chr2,/*genestrand*/'-',
					      outside2,inside2,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_1,min_donor_absdist_1,acceptor_sites_2,min_acceptor_absdist_2) == true) {
	  print_sites(donor_sites_1,chr1,/*donor_genestrand*/'-',outside1,inside1,
		      acceptor_sites_2,chr2,/*acceptor_genestrand*/'-',outside2,inside2);
	}
	Site_gc(&acceptor_sites_2);
	Site_gc(&donor_sites_1);
      }

      if (case4p == true) {
	/* Case 4: Acceptor plus, Donor plus */
	donor_sites_2 = get_donor_sites(&nbadchrs,&min_donor_absdist_2,chr2,/*genestrand*/'+',
					outside2,inside2,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_1 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_1,chr1,/*genestrand*/'+',
					      outside1,inside1,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_2,min_donor_absdist_2,acceptor_sites_1,min_acceptor_absdist_1) == true) {
	  print_sites(donor_sites_2,chr2,/*donor_genestrand*/'+',outside2,inside2,
		      acceptor_sites_1,chr1,/*acceptor_genestrand*/'+',outside1,inside1);
	}
	Site_gc(&acceptor_sites_1);
	Site_gc(&donor_sites_2);
      }

      if (case5p == true) {
	/* Case 5: Donor plus, Acceptor minus */
	donor_sites_1 = get_donor_sites(&nbadchrs,&min_donor_absdist_1,chr1,/*genestrand*/'+',
					outside1,inside1,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_2 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_2,chr2,/*genestrand*/'-',
					      outside2,inside2,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_1,min_donor_absdist_1,acceptor_sites_2,min_acceptor_absdist_2) == true) {
	  print_sites(donor_sites_1,chr1,/*donor_genestrand*/'+',outside1,inside1,
		      acceptor_sites_2,chr2,/*acceptor_genestrand*/'-',outside2,inside2);
	}
	Site_gc(&acceptor_sites_2);
	Site_gc(&donor_sites_1);
      }

      if (case6p == true) {
	/* Case 6: Acceptor minus, Donor plus */
	donor_sites_2 = get_donor_sites(&nbadchrs,&min_donor_absdist_2,chr2,/*genestrand*/'+',
					outside2,inside2,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_1 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_1,chr1,/*genestrand*/'-',
					      outside1,inside1,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_2,min_donor_absdist_2,acceptor_sites_1,min_acceptor_absdist_1) == true) {
	  print_sites(donor_sites_2,chr2,/*donor_genestrand*/'+',outside2,inside2,
		      acceptor_sites_1,chr1,/*acceptor_genestrand*/'-',outside1,inside1);
	}
	Site_gc(&acceptor_sites_1);
	Site_gc(&donor_sites_2);
      }

      if (case7p == true) {
	/* Case 7: Donor minus, Acceptor plus */
	donor_sites_1 = get_donor_sites(&nbadchrs,&min_donor_absdist_1,chr1,/*genestrand*/'-',
					outside1,inside1,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_2 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_2,chr2,/*genestrand*/'+',
					      outside2,inside2,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_1,min_donor_absdist_1,acceptor_sites_2,min_acceptor_absdist_2) == true) {
	  print_sites(donor_sites_1,chr1,/*donor_genestrand*/'-',outside1,inside1,
		      acceptor_sites_2,chr2,/*acceptor_genestrand*/'+',outside2,inside2);
	}
	Site_gc(&acceptor_sites_2);
	Site_gc(&donor_sites_1);
      }

      if (case8p == true) {
	/* Case 8: Acceptor plus, Donor minus */
	donor_sites_2 = get_donor_sites(&nbadchrs,&min_donor_absdist_2,chr2,/*genestrand*/'-',
					outside2,inside2,map_iit,ignore_transcripts_table,
					novelsites_iit,donor_typeint);
	acceptor_sites_1 = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist_1,chr1,/*genestrand*/'+',
					      outside1,inside1,map_iit,ignore_transcripts_table,
					      novelsites_iit,acceptor_typeint);
	if (sites_close_p(donor_sites_2,min_donor_absdist_2,acceptor_sites_1,min_acceptor_absdist_1) == true) {
	  print_sites(donor_sites_2,chr2,/*donor_genestrand*/'-',outside2,inside2,
		      acceptor_sites_1,chr1,/*acceptor_genestrand*/'+',outside1,inside1);
	}
	Site_gc(&acceptor_sites_1);
	Site_gc(&donor_sites_2);
      }
    }
  }

  if (nbadchrs > 0) {
    fprintf(stderr,"%d lines skipped because chromosome not found in map file\n",nbadchrs);
  }

  if (novelsites_iit != NULL) {
    IIT_free(&novelsites_iit);
  }

  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }

  if (ignore_transcripts_table != NULL) {
    nkeys = Table_length(ignore_transcripts_table);
    keys = (char **) Table_keys(ignore_transcripts_table,NULL);
    for (i = 0; i < nkeys; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
    Table_free(&ignore_transcripts_table);
  }

  return 0;
}

