static char rcsid[] = "$Id: add-exons-spliced.c 143395 2014-08-05 16:02:26Z twu $";
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

#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define uabs(x,y) ((x < y) ? y - x : x - y)

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;

static char *ignore_transcripts_file = NULL;
static int nflanking = 5;

static char *sample = "sample";
static char *read_source = "local"; /* local, distant, localalt, or distantalt */
static char *gene_prefix = "X";	    /* Used when no annotation is provided */


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Known genes */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  {"ignore", required_argument, 0, 'I'}, /* ignore_transcripts_file */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */

  {"sample", required_argument, 0, 's'}, /* sample */
  {"source", required_argument, 0, 't'}, /* read_source */

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
Site_print_genes (FILE *fp, List_T list) {
  Site_T site;
  char **array;
  int n, i;
  List_T printed = NULL, p, q;
  bool seenp;

  for (p = list; p != NULL; p = List_next(p)) {
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
  fprintf(fp,"%s",array[0]);
  for (i = 1; i < n; i++) {
    fprintf(fp,",%s",array[i]);
  }
  FREE(array);

  return;
}



static void
Site_print_donor_exons (FILE *fp, List_T list) {
  Site_T site;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (p != list) {
      fprintf(fp,",");
    }
    if (site->novelp == false) {
      fprintf(fp,"%s_exon%d/%d",site->acc,site->exoni,site->nexons);
    } else {
      fprintf(fp,"*%s_exon%d/%d",site->acc,site->exoni,site->nexons);
      if (site->dist > 0) {
	fprintf(fp,"+%d",site->dist);
      } else {
	fprintf(fp,"%d",site->dist);
      }
    }
  }

  return;
}

static void
Site_print_acceptor_exons (FILE *fp, List_T list) {
  Site_T site;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    site = (Site_T) List_head(p);
    if (p != list) {
      fprintf(fp,",");
    }
    if (site->novelp == false) {
      fprintf(fp,"%s_exon%d/%d",site->acc,site->exoni,site->nexons);
    } else {
      fprintf(fp,"*%s_exon%d/%d",site->acc,site->exoni,site->nexons);
      if (site->dist > 0) {
	fprintf(fp,"+%d",site->dist);
      } else {
	fprintf(fp,"%d",site->dist);
      }
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
  Genomicpos_T best_splicesite, exonstart, exonend;
  int dist;
  int sign;
  int gene_namelength;
  int best_exoni, exoni = 0, nexons;
  bool novelp;

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

      } else if (chrpos == exonend) {
	*min_absdist = 0;
	best_exoni = exoni;
	best_splicesite = exonend;

      } else if (uabs(chrpos,exonend) < *min_absdist) {
	*min_absdist = uabs(chrpos,exonend);
	best_exoni = exoni;
	best_splicesite = exonend;
      }
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

      } else if (chrpos == exonend) {
	*min_absdist = 0;
	best_exoni = exoni;
	best_splicesite = exonend;

      } else if (uabs(exonend,chrpos) < *min_absdist) {
	*min_absdist = uabs(exonend,chrpos);
	best_exoni = exoni;
	best_splicesite = exonend;
      }
    }
  }

  Uintlist_free(&exonstarts);
  Uintlist_free(&exonends);

  if (exoni == 1) {
    return NULL;
  } else {
    if (*min_absdist == 0) {
      dist = 0;
      novelp = false;
    } else if (sign > 0) {
      dist = chrpos - best_splicesite;
      novelp = true;
    } else {
      dist = best_splicesite - chrpos;
      novelp = true;
    }

    debug(printf("Best donor is acc %s, position %u, exon %d out of %d, with dist %d\n",
		 acc,best_splicesite,best_exoni,nexons,dist));
    return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		    best_exoni,nexons,novelp,flankingp);
  }
}


static Site_T
get_best_acceptor_known (Genomicpos_T *min_absdist, char *acc, char *annot, Genomicpos_T chrpos,
			 bool flankingp) {
  Uintlist_T exonstarts, exonends, p, q;
  Genomicpos_T best_splicesite, exonstart, exonend;
  int dist;
  int sign;
  int gene_namelength;
  int best_exoni, exoni, nexons;
  bool novelp;

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

      } else if (chrpos == exonstart) {
	*min_absdist = 0;
	best_exoni = exoni;
	best_splicesite = exonstart;

      } else if (uabs(chrpos, exonstart) < *min_absdist) {
	*min_absdist = uabs(chrpos,exonstart);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }
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
	  
      } else if (chrpos == exonstart) {
	*min_absdist = 0;
	best_exoni = exoni;
	best_splicesite = exonstart;
	
      } else if (uabs(exonstart,chrpos) < *min_absdist) {
	*min_absdist = uabs(exonstart,chrpos);
	best_exoni = exoni;
	best_splicesite = exonstart;
      }
    }
  }

  Uintlist_free(&exonstarts);
  Uintlist_free(&exonends);

  if (exoni == nexons) {
    return NULL;
  } else {
    if (*min_absdist == 0) {
      dist = 0;
      novelp = false;
    } else if (sign > 0) {
      dist = chrpos - best_splicesite;
      novelp = true;
    } else {
      dist = best_splicesite - chrpos;
      novelp = true;
    }

    debug(printf("Best acceptor is acc %s, position %u, exon %d out of %d, with dist %d\n",
		 acc,best_splicesite,best_exoni,nexons,dist));
    return Site_new(dist,best_splicesite,/*gene*/annot,gene_namelength,acc,
		    best_exoni,nexons,novelp,flankingp);
  }
}


static List_T
get_donor_sites (int *nbadchrs, Genomicpos_T *min_absdist,
		 List_T *donor_sites_left, List_T *donor_sites_right,
		 char *divstring, char genestrand,
		 Genomicpos_T outside, Genomicpos_T inside, IIT_T map_iit,
		 Table_T ignore_transcripts_table) {
  List_T donor_sites = NULL;
  int nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int *matches;
  Genomicpos_T coordstart, coordend;
  Genomicpos_T min_absdist_left, min_absdist_right, absdist;
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
  *donor_sites_left = (List_T) NULL;
  *donor_sites_right = (List_T) NULL;
  *min_absdist = -1U;

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
    IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,divstring,
		     coordstart,coordend,nflanking,sign);
    debug(printf("Getting flanking with %d on left and %d on right\n",nleftflanks,nrightflanks));

    for (i = 0; i < nleftflanks; i++) {
      acc = IIT_label(map_iit,leftflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,map_iit,leftflanks[i],&allocp);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_best_donor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (*donor_sites_left == NULL) {
	*donor_sites_left = List_push(NULL,(void *) site);
	min_absdist_left = absdist;
      } else if (absdist < min_absdist_left) {
	Site_gc(&*donor_sites_left);
	*donor_sites_left = List_push(NULL,(void *) site);
	min_absdist_left = absdist;
      } else if (absdist == min_absdist_left) {
	*donor_sites_left = List_push(*donor_sites_left,(void *) site);
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
      } else if (*donor_sites_right == NULL) {
	*donor_sites_right = List_push(NULL,(void *) site);
	min_absdist_right = absdist;
      } else if (absdist < min_absdist_right) {
	Site_gc(&*donor_sites_right);
	*donor_sites_right = List_push(NULL,(void *) site);
	min_absdist_right = absdist;
      } else if (absdist == min_absdist_right) {
	*donor_sites_right = List_push(*donor_sites_right,(void *) site);
      } else {
	Site_free(&site);
      }

      if (allocp) {
	FREE(restofheader);
      }
    }
    FREE(rightflanks);
  }

  return donor_sites;
}

static List_T
get_acceptor_sites (int *nbadchrs, Genomicpos_T *min_absdist,
		    List_T *acceptor_sites_left, List_T *acceptor_sites_right,
		    char *divstring, char genestrand,
		    Genomicpos_T outside, Genomicpos_T inside, IIT_T map_iit,
		    Table_T ignore_transcripts_table) {
  List_T acceptor_sites = NULL;
  int nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int *matches;
  Genomicpos_T coordstart, coordend;
  Genomicpos_T min_absdist_left, min_absdist_right, absdist;
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
  *acceptor_sites_left = (List_T) NULL;
  *acceptor_sites_right = (List_T) NULL;
  *min_absdist = -1U;

  if (IIT_divint(map_iit,divstring) < 0) {
    if (*nbadchrs == 10) {
      fprintf(stderr,"...more than 10 warning messages.  No longer printing them\n");
    } else if (*nbadchrs > 10) {
      /* Skip */
    } else {
      fprintf(stderr,"chromosome %s not found in map_iit file\n",divstring);
    }
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
    IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,divstring,
		     coordstart,coordend,nflanking,sign);
    debug(printf("Getting flanking with %d on left and %d on right\n",nleftflanks,nrightflanks));

    for (i = 0; i < nleftflanks; i++) {
      acc = IIT_label(map_iit,leftflanks[i],&allocp);
      annot = IIT_annotation(&restofheader,map_iit,leftflanks[i],&allocp);

      if (ignore_transcripts_table != NULL && Table_get(ignore_transcripts_table,(void *) acc) != 0) {
	/* Skip */
      } else if ((site = get_best_acceptor_known(&absdist,acc,annot,inside,/*flankingp*/true)) == NULL) {
	/* Skip */
      } else if (*acceptor_sites_left == NULL) {
	*acceptor_sites_left = List_push(NULL,(void *) site);
	min_absdist_left = absdist;
      } else if (absdist < min_absdist_left) {
	Site_gc(&*acceptor_sites_left);
	*acceptor_sites_left = List_push(NULL,(void *) site);
	min_absdist_left = absdist;
      } else if (absdist == min_absdist_left) {
	*acceptor_sites_left = List_push(*acceptor_sites_left,(void *) site);
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
      } else if (*acceptor_sites_right == NULL) {
	*acceptor_sites_right = List_push(NULL,(void *) site);
	min_absdist_right = absdist;
      } else if (absdist < min_absdist_right) {
	Site_gc(&*acceptor_sites_right);
	*acceptor_sites_right = List_push(NULL,(void *) site);
	min_absdist_right = absdist;
      } else if (absdist == min_absdist_right) {
	*acceptor_sites_right = List_push(*acceptor_sites_right,(void *) site);
      } else {
	Site_free(&site);
      }

      if (allocp) {
	FREE(restofheader);
      }
    }
    FREE(rightflanks);
  }

  return acceptor_sites;
}



static void
print_sites (FILE *fp, Genomicpos_T min_donor_absdist, List_T donor_sites,
	     List_T donor_sites_left, List_T donor_sites_right,
	     char *donor_chr, char donor_genestrand,
	     Genomicpos_T donor_outside, Genomicpos_T donor_inside,
	     Genomicpos_T min_acceptor_absdist, List_T acceptor_sites,
	     List_T acceptor_sites_left, List_T acceptor_sites_right,
	     char *acceptor_chr, char acceptor_genestrand,
	     Genomicpos_T acceptor_outside, Genomicpos_T acceptor_inside,
	     int count) {
  int dist;

  if (donor_sites == NULL && donor_sites_left == NULL && donor_sites_right == NULL) {
    return;
  }
  if (acceptor_sites == NULL && acceptor_sites_left == NULL && acceptor_sites_right == NULL) {
    return;
  }

  fprintf(fp,"%s\t%c\t%u\t%u\t",donor_chr,donor_genestrand,donor_outside,donor_inside);

  donor_sites = Site_pos_acc_sort(donor_sites);
  donor_sites_left = Site_pos_acc_sort(donor_sites_left);
  donor_sites_right = Site_pos_acc_sort(donor_sites_right);

  if (donor_sites != NULL) {
    if ((dist = ((Site_T) donor_sites->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (donor_sites_left != NULL && donor_sites_right != NULL) {
    if ((dist = ((Site_T) donor_sites_left->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
    fprintf(fp,",");
    if ((dist = ((Site_T) donor_sites_right->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (donor_sites_left != NULL) {
    if ((dist = ((Site_T) donor_sites_left->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (donor_sites_right != NULL) {
    if ((dist = ((Site_T) donor_sites_right->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  }
  fprintf(fp,"\t");

  acceptor_sites = Site_pos_acc_sort(acceptor_sites);
  acceptor_sites_left = Site_pos_acc_sort(acceptor_sites_left);
  acceptor_sites_right = Site_pos_acc_sort(acceptor_sites_right);

  fprintf(fp,"%s\t%c\t%u\t%u\t",acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside);

  if (acceptor_sites != NULL) {
    if ((dist = ((Site_T) acceptor_sites->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (acceptor_sites_left != NULL && acceptor_sites_right != NULL) {
    if ((dist = ((Site_T) acceptor_sites_left->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
    fprintf(fp,",");
    if ((dist = ((Site_T) acceptor_sites_right->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (acceptor_sites_left != NULL) {
    if ((dist = ((Site_T) acceptor_sites_left->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  } else if (acceptor_sites_right != NULL) {
    if ((dist = ((Site_T) acceptor_sites_right->first)->dist) == 0) {
      fprintf(fp,"0");
    } else if (dist > 0) {
      fprintf(fp,"+%d",dist);
    } else {
      fprintf(fp,"%d",dist);
    }
  }
  fprintf(fp,"\t");

  if (min_donor_absdist != 0U) {
    fprintf(fp,"Novel@");
  }
  if (donor_sites != NULL) {
    Site_print_genes(fp,donor_sites);
  } else if (donor_sites_left != NULL && donor_sites_right != NULL) {
    Site_print_genes(fp,donor_sites_left);
    fprintf(fp,"&");
    Site_print_genes(fp,donor_sites_right);
  } else if (donor_sites_left != NULL) {
    Site_print_genes(fp,donor_sites_left);
  } else if (donor_sites_right != NULL) {
    Site_print_genes(fp,donor_sites_right);
  }
  fprintf(fp,"\t");

  if (min_acceptor_absdist != 0U) {
    fprintf(fp,"Novel@");
  }
  if (acceptor_sites != NULL) {
    Site_print_genes(fp,acceptor_sites);
  } else if (acceptor_sites_left != NULL && acceptor_sites_right != NULL) {
    Site_print_genes(fp,acceptor_sites_left);
    fprintf(fp,"&");
    Site_print_genes(fp,acceptor_sites_right);
  } else if (acceptor_sites_left != NULL) {
    Site_print_genes(fp,acceptor_sites_left);
  } else if (acceptor_sites_right != NULL) {
    Site_print_genes(fp,acceptor_sites_right);
  }
  fprintf(fp,"\t");

  if (donor_sites != NULL) {
    Site_print_donor_exons(fp,donor_sites);
  } else if (donor_sites_left != NULL && donor_sites_right != NULL) {
    Site_print_donor_exons(fp,donor_sites_left);
    fprintf(fp,",");
    Site_print_donor_exons(fp,donor_sites_right);
  } else if (donor_sites_left != NULL) {
    Site_print_donor_exons(fp,donor_sites_left);
  } else if (donor_sites_right != NULL) {
    Site_print_donor_exons(fp,donor_sites_right);
  }
  fprintf(fp,"\t");

  if (acceptor_sites != NULL) {
    Site_print_acceptor_exons(fp,acceptor_sites);
  } else if (acceptor_sites_left != NULL && acceptor_sites_right != NULL) {
    Site_print_acceptor_exons(fp,acceptor_sites_left);
    fprintf(fp,",");
    Site_print_acceptor_exons(fp,acceptor_sites_right);
  } else if (acceptor_sites_left != NULL) {
    Site_print_acceptor_exons(fp,acceptor_sites_left);
  } else if (acceptor_sites_right != NULL) {
    Site_print_acceptor_exons(fp,acceptor_sites_right);
  }
  fprintf(fp,"\t");

  fprintf(fp,"%d",count);
  fprintf(fp,"\n");

  List_free(&donor_sites);
  List_free(&donor_sites_left);
  List_free(&donor_sites_right);
  List_free(&acceptor_sites);
  List_free(&acceptor_sites_left);
  List_free(&acceptor_sites_right);

  return;
}


static void
print_site_novel (int *genei, FILE *fp, char *donor_chr, char donor_genestrand,
		  Genomicpos_T donor_outside, Genomicpos_T donor_inside,
		  char *acceptor_chr, char acceptor_genestrand,
		  Genomicpos_T acceptor_outside, Genomicpos_T acceptor_inside,
		  int count) {
  int dist;
  int geneX = ++(*genei);
  int geneY = ++(*genei);

  fprintf(fp,"%s\t%c\t%u\t%u\t",donor_chr,donor_genestrand,donor_outside,donor_inside);
  fprintf(fp,"NA\t");

  fprintf(fp,"%s\t%c\t%u\t%u\t",acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside);
  fprintf(fp,"NA\t");


  fprintf(fp,"Novel@Gene%s%d\t",gene_prefix,geneX);
  fprintf(fp,"Novel@Gene%s%d\t",gene_prefix,geneY);

  fprintf(fp,"*Gene%s%d_exon_0/0\t",gene_prefix,geneX);
  fprintf(fp,"*Gene%s%d_exon_0/0\t",gene_prefix,geneY);

  fprintf(fp,"%d",count);
  fprintf(fp,"\n");

  return;
}




int
main (int argc, char *argv[]) {
  FILE *fp, *known_fp, *novel_fp;
  char *iitfile, *filename;

  char *genomesubdir = NULL, *fileroot = NULL, *mapdir = NULL;
  IIT_T map_iit = NULL;
  Table_T ignore_transcripts_table = NULL;
  char **keys;
  int nkeys, i;
  int genei = 0;

  char Buffer[1024], acc[1000], donor_chr[1000], acceptor_chr[1000], *copy;
  char donor_genestrand, acceptor_genestrand;
  Genomicpos_T donor_outside, donor_inside, acceptor_outside, acceptor_inside;
  Genomicpos_T min_donor_absdist, min_acceptor_absdist;
  List_T donor_sites, donor_sites_left, donor_sites_right,
    acceptor_sites, acceptor_sites_left, acceptor_sites_right;
  int nbadchrs = 0, count;


  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:M:m:s:t:I:u:g:",
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

    case 'I': ignore_transcripts_file = optarg; break;
    case 'u': nflanking = atoi(optarg); break;

    case 's': sample = optarg; break;
    case 't': read_source = optarg; break;

    case 'g': gene_prefix = optarg; break;

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


  filename = (char *) CALLOC(strlen(sample)+strlen(".")+strlen(read_source)+strlen(".known.exons")+1,
			     sizeof(char));
  sprintf(filename,"%s.%s.known.exons",sample,read_source);
  known_fp = fopen(filename,"w");
  FREE(filename);

  filename = (char *) CALLOC(strlen(sample)+strlen(".")+strlen(read_source)+strlen(".novel.exons")+1,
			     sizeof(char));
  sprintf(filename,"%s.%s.novel.exons",sample,read_source);
  novel_fp = fopen(filename,"w");
  FREE(filename);

  while (fgets(Buffer,1024,stdin) != NULL) {
    if (Buffer[0] == '#') {
      /* Skip line */
    } else if (sscanf(Buffer,"%s %c %u %u %s %c %u %u %d",
	       donor_chr,&donor_genestrand,&donor_outside,&donor_inside,
	       acceptor_chr,&acceptor_genestrand,&acceptor_outside,&acceptor_inside,
	       &count) == 9) {
      debug(printf("Read %s\n",Buffer));

      if (map_iit == NULL) {
	print_site_novel(&genei,novel_fp,donor_chr,donor_genestrand,donor_outside,donor_inside,
			 acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside,
			 count);

      } else {
	donor_sites = get_donor_sites(&nbadchrs,&min_donor_absdist,&donor_sites_left,&donor_sites_right,
				      donor_chr,donor_genestrand,donor_outside,donor_inside,
				      map_iit,ignore_transcripts_table);
	acceptor_sites = get_acceptor_sites(&nbadchrs,&min_acceptor_absdist,&acceptor_sites_left,&acceptor_sites_right,
					    acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside,
					    map_iit,ignore_transcripts_table);

	debug(printf("min_donor_absdist %u and min_acceptor_absdist %u\n",
		     min_donor_absdist,min_acceptor_absdist));
	if (min_donor_absdist == 0 && min_acceptor_absdist == 0) {
	  print_sites(known_fp,min_donor_absdist,donor_sites,
		      donor_sites_left,donor_sites_right,
		      donor_chr,donor_genestrand,donor_outside,donor_inside,
		      min_acceptor_absdist,acceptor_sites,
		      acceptor_sites_left,acceptor_sites_right,
		      acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside,
		      count);
	} else {
	  print_sites(novel_fp,min_donor_absdist,donor_sites,
		      donor_sites_left,donor_sites_right,
		      donor_chr,donor_genestrand,donor_outside,donor_inside,
		      min_acceptor_absdist,acceptor_sites,
		      acceptor_sites_left,acceptor_sites_right,
		      acceptor_chr,acceptor_genestrand,acceptor_outside,acceptor_inside,
		      count);
	}

	Site_gc(&acceptor_sites_right);
	Site_gc(&acceptor_sites_left);
	Site_gc(&acceptor_sites);
	Site_gc(&donor_sites_right);
	Site_gc(&donor_sites_left);
	Site_gc(&donor_sites);
      }
    }
  }

  fclose(novel_fp);
  fclose(known_fp);

  if (nbadchrs > 0) {
    fprintf(stderr,"%d lines skipped because chromosome not found in map file\n",nbadchrs);
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

