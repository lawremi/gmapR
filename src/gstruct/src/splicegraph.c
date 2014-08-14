static char rcsid[] = "$Id: splicegraph.c 136513 2014-05-16 17:58:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splicegraph.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For isdigit */
#include <unistd.h>		/* For getopt */
#include <math.h>		/* For sqrt */
#include <limits.h>		/* For INT_MAX */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include "mem.h"
#include "fopen.h"

#include "uintlist.h"
#include "chrnum.h"
#include "interval.h"
#include "iit-read.h"
#include "uinttable.h"

#include "splice.h"
#include "cappaths.h"
#include "listdef.h"
#include "intlist.h"


#define RETAINED_BUFFER 100U	/* Must be bigger than window size used to compute window_diff */
#define PSEUDOCOUNT 0
#define MIN_COUNT 10


typedef int Score_T;

/* For dynamic programming */
/* if (site->chrpos == 51015289) x */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Adding exons */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Adding introns */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Terminals */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Homology */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Alternate */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* flat_chrpos_low and flat_chrpos_high */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* According to Sakharkar, 2004, max exon length in human is 11,923
   bp.  However, that is probably an UTR */


#ifdef COMPUTE_DONOR_BOUND
#define ALTSITE_DISTANCE 20
#endif


static bool print_mainalt_p = true;

static double slope_threshold = 0.1;
static int length_threshold = 50;

static int highlow_length = 70;

static int halfdata_genebounds = 20;
static int halfdata_cappaths = 10;

static int short_form_support = 5;
static int long_form_support = 10;



#if 0
/************************************************************************
 *   Sites
 ************************************************************************/

typedef struct Site_T *Site_T;
struct Site_T {
  Genomicpos_T chrpos;
};

static void
Site_free (Site_T *old) {
  FREE(*old);
  return;
}


static Site_T
Site_new (Genomicpos_T chrpos) {
  Site_T new = (Site_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  return new;
}


static int
Site_ascending_cmp (const void *a, const void *b) {
  Site_T x = * (Site_T *) a;
  Site_T y = * (Site_T *) b;

  if (x->chrpos < y->chrpos) {
    return -1;
  } else if (y->chrpos < x->chrpos) {
    return 1;
  } else {
    return 0;
  }
}

static int
Site_descending_cmp (const void *a, const void *b) {
  Site_T x = * (Site_T *) a;
  Site_T y = * (Site_T *) b;

  if (x->chrpos > y->chrpos) {
    return -1;
  } else if (y->chrpos > x->chrpos) {
    return 1;
  } else {
    return 0;
  }
}
#endif


typedef enum {NONE, MAIN, ALT_DONOR, ALT_ACCEPTOR,
	      SKIPPED_EXON, SKIPPED_DOUBLE_EXON, EXTRA_EXON, SWITCH_EXON, SWITCH_INTRON,
	      RETAINED_INTRON, INITIAL_ALT_EXON, INITIAL_SHORT, INITIAL_LONG,
	      TERMINAL_ALT_SITE, TERMINAL_ALT_EXON, TERMINAL_SHORT, TERMINAL_LONG,
	      TERMINAL_RETAINED} Introntype_T;

static char *
Introntype_string (Introntype_T introntype) {
  switch (introntype) {
  case NONE: return "none";
  case MAIN: return "ref";
  case ALT_DONOR: return "A5SS"; /* Also in MISO */
  case ALT_ACCEPTOR: return "A3SS"; /* Also in MISO */
  case SKIPPED_EXON: return "SE"; /* Also in MISO */
  case SKIPPED_DOUBLE_EXON: return "S2E";
  case EXTRA_EXON: return "EE";
  case SWITCH_EXON: return "MXE"; /* Also in MISO */
  case SWITCH_INTRON: return "MXI";
  case RETAINED_INTRON: return "RI"; /* Also in MISO */
  case INITIAL_ALT_EXON: return "AFE"; /* Also in MISO */
  case INITIAL_SHORT: return "IS";
  case INITIAL_LONG: return "IL";
  case TERMINAL_ALT_SITE: return "A3SS"; /* Just treat as a standard alt 3' splice site */
  case TERMINAL_ALT_EXON: return "ALE"; /* Also in MISO */
  case TERMINAL_SHORT: return "TS";
  case TERMINAL_LONG: return "TL";
  case TERMINAL_RETAINED: return "RI"; /* Just call it a retained intron, instead of terminal retained */
  }
  return "";
}



typedef struct Donor_T *Donor_T;
typedef struct Acceptor_T *Acceptor_T;
typedef struct Intron_T *Intron_T;
typedef struct Exon_T *Exon_T;
typedef struct Donor_alt_T *Donor_alt_T;
typedef struct Acceptor_alt_T *Acceptor_alt_T;
typedef struct Path_T *Path_T;


struct Donor_alt_T {
  bool knownp;
  Acceptor_T prev_site;
  Introntype_T type;
  Exon_T alt_exon;

  Intron_T main_intron;
  Intron_T alt_intron;
};


static void
Donor_alt_free (Donor_alt_T *old) {
  FREE(*old);
  return;
}

static Donor_alt_T
Donor_alt_new (Acceptor_T prev_site, Introntype_T type, Exon_T alt_exon,
	       Intron_T main_intron, Intron_T alt_intron, bool knownp) {
  Donor_alt_T new = (Donor_alt_T) MALLOC(sizeof(*new));
  
  new->knownp = knownp;
  new->prev_site = prev_site;
  new->type = type;
  new->alt_exon = alt_exon;
  new->main_intron = main_intron;
  new->alt_intron = alt_intron;

  return new;
}


struct Acceptor_alt_T {
  bool knownp;
  Donor_T prev_site;
  Introntype_T type;

  Intron_T main_intron;
  Intron_T main_intron_2;
  Intron_T main_intron_3;
  Intron_T alt_intron;		/* For alt_prev, alt_prev_extra, and alt_prev_switch */
  Intron_T alt_intron_2;	/* For alt_prev_extra and alt_prev_switch */
};

static void
Acceptor_alt_free (Acceptor_alt_T *old) {
  FREE(*old);
  return;
}


static Acceptor_alt_T
Acceptor_alt_new (Donor_T prev_site, Introntype_T type,
		  Intron_T main_intron, Intron_T main_intron_2, Intron_T main_intron_3,
		  Intron_T alt_intron, Intron_T alt_intron_2, bool knownp) {
  Acceptor_alt_T new = (Acceptor_alt_T) MALLOC(sizeof(*new));

  new->knownp = knownp;
  new->prev_site = prev_site;
  new->type = type;
  new->main_intron = main_intron;
  new->main_intron_2 = main_intron_2;
  new->main_intron_3 = main_intron_3;
  new->alt_intron = alt_intron;
  new->alt_intron_2 = alt_intron_2;

  return new;
}

struct Donor_T {
  Genomicpos_T chrpos;
  bool genebound_p;
  bool switchp;
  
  Uinttable_T exons_usedp;

  List_T exons;

  bool utr_p;
  Score_T points;
  Score_T score;

  Genomicpos_T usedp;

  Acceptor_T prev;
  List_T alt_structs;
  List_T alt_sites;
  List_T alt_switches;

  Genomicpos_T known_prev_acceptorpos;
  Exon_T main_exon;		/* Points backward to previous acceptor */
};


struct Acceptor_T {
  Genomicpos_T chrpos;
  bool genebound_p;
  bool switchp;
  
  Uinttable_T splices_table;

  List_T introns;

  bool utr_p;
  Score_T points;
  Score_T score;

  bool known_terminalp;
  bool terminalp;

  Path_T known_mainpath;		/* Found from processing known genes */
  Path_T known_termpath;		/* Found from processing known genes */
  List_T known_altpaths;		/* Found from processing known genes */
  List_T known_alt_termpaths;		/* Found from processing known genes */
  List_T known_ret_termpaths;		/* Found from processing known genes */
  Genomicpos_T usedp;
  bool paths_usedp;

  Donor_T prev;
  List_T alt_structs;
  List_T alt_sites;
  List_T alt_switches;

  Intron_T main_intron;
  Exon_T main_exon;		/* Points forward to next donor */

  Genomicpos_T retained_geneend;
};


/* ups get default utr_p of false.  When they don't have any exons, we
   assign them utr_p true.  downs get default utr_p of true.  When an
   intron comes into one, we assign it utr_p false. */


static Donor_T
Donor_new (Genomicpos_T chrpos, Genomicpos_T known_prev_acceptorpos) {
  Donor_T new = (Donor_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->genebound_p = false;
  new->switchp = false;

  new->exons_usedp = Uinttable_new(100);

  new->utr_p = false;		/* false for donors, true for acceptors */
  new->points = 0;
  new->score = 0;

  new->usedp = 0U;

  new->prev = (Acceptor_T) NULL;
  new->exons = (List_T) NULL;

  new->alt_structs = (List_T) NULL;
  new->alt_sites = (List_T) NULL;
  new->alt_switches = (List_T) NULL;

  new->known_prev_acceptorpos = known_prev_acceptorpos;
  new->main_exon = (Exon_T) NULL;

  return new;
}

static Acceptor_T
Acceptor_new (Genomicpos_T chrpos) {
  Acceptor_T new = (Acceptor_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->genebound_p = false;
  new->switchp = false;

  new->splices_table = Uinttable_new(100);

  new->utr_p = true;		/* true for acceptors, false for donors */
  new->points = 0;
  new->score = 0;

  new->known_terminalp = false;
  new->terminalp = true;

  new->known_mainpath = (Path_T) NULL;
  new->known_termpath = (Path_T) NULL;
  new->known_altpaths = (List_T) NULL;
  new->known_alt_termpaths = (List_T) NULL;
  new->known_ret_termpaths = (List_T) NULL;
  new->usedp = 0U;

  new->prev = (Donor_T) NULL;
  new->introns = (List_T) NULL;

  new->alt_structs = (List_T) NULL;
  new->alt_sites = (List_T) NULL;
  new->alt_switches = (List_T) NULL;

  new->main_intron = (Intron_T) NULL;
  new->main_exon = (Exon_T) NULL;

  new->retained_geneend = 0U;

  return new;
}



static int
Acceptor_score_cmp (const void *a, const void *b) {
  Acceptor_T x = * (Acceptor_T *) a;
  Acceptor_T y = * (Acceptor_T *) b;

  if (x->known_mainpath != NULL && y->known_mainpath == NULL) {
    return -1;
  } else if (y->known_mainpath != NULL && x->known_mainpath == NULL) {
    return +1;
  } else if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else {
    return 0;
  }
}



/************************************************************************
 *   Introns
 ************************************************************************/

/* alt_exon means the same as alt_acceptor */
/* alt_intron means the same as alt_donor */

static bool
Pathtype_initialp (Introntype_T type) {
  if (type == INITIAL_ALT_EXON) {
    return true;
  } else if (type == INITIAL_SHORT) {
    return true;
  } else if (type == INITIAL_LONG) {
    return true;
  } else {
    return false;
  }
}

static bool
Pathtype_terminalp (Introntype_T type) {
  if (type == TERMINAL_ALT_SITE) {
    return true;
  } else if (type == TERMINAL_ALT_EXON) {
    return true;
  } else if (type == TERMINAL_SHORT) {
    return true;
  } else if (type == TERMINAL_LONG) {
    return true;
  } else {
    return false;
  }
}
	      


#if 0
typedef enum {NO_TERMINAL_DIFF, TERMINAL_ALT_SITE, TERMINAL_ALT_EXON, TERMINAL_SHORT, TERMINAL_LONG} Terminaltype_T;

static char *
Terminaltype_string (Terminaltype_T terminaltype) {
  switch (terminaltype) {
  case NO_TERMINAL_DIFF: return "";
  case TERMINAL_ALT_SITE: return "alt_acceptor";
  case TERMINAL_ALT_EXON: return "alt_last";
  case TERMINAL_SHORT: return "terminal_short";
  case TERMINAL_LONG: return "terminal_log";
  }
  return "";
}
#endif


struct Intron_T {
  Splice_T splice;
  bool validp;
  bool known_retained_p;

  Introntype_T introntype;
  int altcount1;
  int altcount2;

  Donor_T prev_site;
  Acceptor_T site;

  Genomicpos_T donorpos;
  Genomicpos_T acceptorpos;

  bool forwardp;		/* For conflict resolution */
  Genomicpos_T low;
  Genomicpos_T high;
  Score_T pathscore;		/* For conflict resolution */

  int *genei;			/* For known introns found in introns_iit, but using types (genei), not matches[i] */
  int ngenei;
};


static Intron_T
Intron_new (Donor_T prev_site, Acceptor_T site, Splice_T splice,
	    IIT_T introns_iit, char *chr, bool forwardp) {
  Intron_T new = (Intron_T) MALLOC(sizeof(*new));

  new->splice = splice;
  new->validp = true;
  new->known_retained_p = false;

  new->introntype = NONE;
  new->altcount1 = 0;
  new->altcount2 = 0;

  new->prev_site = prev_site;
  new->site = site;
  new->forwardp = forwardp;
  new->pathscore = 0.0;

  new->donorpos = prev_site->chrpos;
  new->acceptorpos = site->chrpos;

  if (new->forwardp == true) {
    new->low = prev_site->chrpos;
    new->high = site->chrpos;
  } else {
    new->low = site->chrpos;
    new->high = prev_site->chrpos;
  }

  if (introns_iit == NULL) {
    new->genei = (int *) NULL;
    new->ngenei = 0;
  } else {
    new->genei = IIT_get_exact_types_multiple(&new->ngenei,introns_iit,chr,new->low,new->high);
  }

  return new;
}

static void
Intron_free (Intron_T *old) {
  if ((*old)->genei != NULL) {
    FREE((*old)->genei);
  }
  FREE(*old);
  return;
}


static void
Intron_path_print (List_T introns) {
  List_T p;
  Intron_T intron;

  for (p = introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    printf("%u %u\n",intron->prev_site->chrpos,intron->site->chrpos);
  }
  printf("\n");

  return;
}



static List_T
Intron_path_trim (List_T introns) {
  Intron_T intron;

  if (introns == NULL) {
    return (List_T) NULL;
  } else {
    while (introns != NULL && Splice_valid_end_p(((Intron_T) List_head(introns))->splice) == false) {
      introns = List_pop(introns,(void **) &intron);
      /* printf("#Trimming intron at %u..%u because not valid end\n",intron->prev_site->chrpos,intron->site->chrpos); */
      /* Intron_free(&intron); */
    }

    introns = List_reverse(introns);
    while (introns != NULL && Splice_valid_end_p(((Intron_T) List_head(introns))->splice) == false) {
      introns = List_pop(introns,(void **) &intron);
      /* printf("#Trimming intron at %u..%u because not valid end\n",intron->prev_site->chrpos,intron->site->chrpos); */
      /* Intron_free(&intron); */
    }

    return List_reverse(introns);
  }
}


static bool
Intron_invalid_end_p (List_T introns) {
  Intron_T intron;
  List_T p;

  for (p = introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    if (Splice_valid_end_p(intron->splice) == false) {
      return true;
    }
  }

  return false;
}



static int
Intron_low_cmp (const void *a, const void *b) {
  Intron_T x = * (Intron_T *) a;
  Intron_T y = * (Intron_T *) b;

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

static int
Intron_high_cmp (const void *a, const void *b) {
  Intron_T x = * (Intron_T *) a;
  Intron_T y = * (Intron_T *) b;

  if (x->high > y->high) {
    return -1;
  } else if (y->high > x->high) {
    return +1;
  } else if (x->low > y->low) {
    return -1;
  } else if (y->low > x->low) {
    return +1;
  } else {
    return 0;
  }
}

static Intron_T
Acceptor_find_intron (Acceptor_T acceptor) {
  Intron_T intron;
  List_T p;

  for (p = acceptor->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    if (intron->prev_site == acceptor->prev) {
      return intron;
    }
  }

  fprintf(stderr,"At acceptor %u, intron not found for donor at %u\n",acceptor->chrpos,acceptor->prev->chrpos);
  abort();
  return (Intron_T) NULL;
}


static bool
Donor_alt_equiv (List_T alts, Intron_T alt_intron) {
  List_T p;
  Donor_alt_T this;

  debug(printf("\nTesting donor equivalence of intron at %u..%u",
	       alt_intron->prev_site->chrpos,alt_intron->site->chrpos));
  for (p = alts; p != NULL; p = List_next(p)) {
    this = (Donor_alt_T) List_head(p);
    if (this->alt_intron == alt_intron) {
      debug(printf(" => true\n"));
      return true;
    }
  }

  debug(printf(" => false\n"));
  return false;
}

static bool
Acceptor_alt_equiv (List_T alts, Intron_T alt_intron, Intron_T alt_intron_2) {
  List_T p;
  Acceptor_alt_T this;

  debug(printf("\nTesting acceptor equivalence of intron at %u..%u",
	       alt_intron->prev_site->chrpos,alt_intron->site->chrpos));
  if (alt_intron_2 != NULL) {
    debug(printf(" and %u..%u",alt_intron_2->prev_site->chrpos,alt_intron_2->site->chrpos));
  }

  for (p = alts; p != NULL; p = List_next(p)) {
    this = (Acceptor_alt_T) List_head(p);
    if (this->alt_intron == alt_intron && this->alt_intron_2 == alt_intron_2) {
      debug(printf(" => true\n"));
      return true;
    }
  }

  debug(printf(" => false\n"));
  return false;
}


/************************************************************************
 *   Exons
 ************************************************************************/

struct Exon_T {
  bool validp;			/* May not be necessary */
  
  Acceptor_T prev_site;
  Donor_T site;

  Score_T incrscore;
};


static Genomicpos_T
Exon_length (Exon_T this) {
  if (this->prev_site->chrpos < this->site->chrpos) {
    return this->site->chrpos - this->prev_site->chrpos;
  } else {
    return this->prev_site->chrpos - this->site->chrpos;
  }
}


static Exon_T
Exon_new (Acceptor_T prev_site, Donor_T site, Score_T incrscore) {
  Exon_T new = (Exon_T) MALLOC(sizeof(*new));
  
  new->validp = true;
  new->prev_site = prev_site;
  new->site = site;
  new->incrscore = incrscore;
  return new;
}

static void
Exon_free (Exon_T *old) {
  FREE(*old);
  return;
}


/************************************************************************
 *   End of sites
 ************************************************************************/

static void
Donor_alt_print (Donor_alt_T this) {
  printf("%s\n",Introntype_string(this->type));
  printf("%u %u\n",this->alt_intron->prev_site->chrpos,this->alt_intron->prev_site->chrpos);
  printf("%u %u\n",this->alt_intron->site->chrpos,this->alt_intron->site->chrpos);
  return;
}

static void
Acceptor_alt_print (Acceptor_alt_T this) {
  printf("%s\n",Introntype_string(this->type));
  if (this->alt_intron_2 != NULL) {
    printf("%u %u\n",this->alt_intron->prev_site->chrpos,this->alt_intron->prev_site->chrpos);
    printf("%u %u\n",this->alt_intron->site->chrpos,this->alt_intron_2->prev_site->chrpos);
    printf("%u %u\n",this->alt_intron_2->site->chrpos,this->alt_intron_2->site->chrpos);

  } else {
    printf("%u %u\n",this->alt_intron->prev_site->chrpos,this->alt_intron->prev_site->chrpos);
    printf("%u %u\n",this->alt_intron->site->chrpos,this->alt_intron->site->chrpos);
  }

  return;
}


static void
Donor_free (Donor_T *old) {
  Exon_T exon;
  Donor_alt_T alt;
  List_T p;

  Uinttable_free(&(*old)->exons_usedp);

  for (p = (*old)->alt_switches; p != NULL; p = List_next(p)) {
    alt = (Donor_alt_T) List_head(p);
    Donor_alt_free(&alt);
  }
  List_free(&(*old)->alt_switches);

  for (p = (*old)->alt_sites; p != NULL; p = List_next(p)) {
    alt = (Donor_alt_T) List_head(p);
    Donor_alt_free(&alt);
  }
  List_free(&(*old)->alt_sites);

  for (p = (*old)->alt_structs; p != NULL; p = List_next(p)) {
    alt = (Donor_alt_T) List_head(p);
    Donor_alt_free(&alt);
  }
  List_free(&(*old)->alt_structs);

  for (p = (*old)->exons; p != NULL; p = List_next(p)) {
    exon = (Exon_T) List_head(p);
    Exon_free(&exon);
  }
  List_free(&((*old)->exons));

  FREE(*old);
  return;
}

static void
Acceptor_free (Acceptor_T *old) {
  Intron_T intron;
  Acceptor_alt_T alt;
  List_T p;

  Uinttable_free(&(*old)->splices_table);

  for (p = (*old)->alt_switches; p != NULL; p = List_next(p)) {
    alt = (Acceptor_alt_T) List_head(p);
    Acceptor_alt_free(&alt);
  }
  List_free(&(*old)->alt_switches);

  for (p = (*old)->alt_sites; p != NULL; p = List_next(p)) {
    alt = (Acceptor_alt_T) List_head(p);
    Acceptor_alt_free(&alt);
  }
  List_free(&(*old)->alt_sites);

  for (p = (*old)->alt_structs; p != NULL; p = List_next(p)) {
    alt = (Acceptor_alt_T) List_head(p);
    Acceptor_alt_free(&alt);
  }
  List_free(&(*old)->alt_structs);

  for (p = (*old)->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    Intron_free(&intron);
  }
  List_free(&((*old)->introns));

  FREE(*old);
  return;
}




/************************************************************************/

static bool
homologousp (Genome_T genome, Genomicpos_T chroffset, Genomicpos_T mainend,
	     Genomicpos_T altstart, Genomicpos_T altend) {
  Genomicpos_T exonlength, mainstart;
  char *main_buffer, *alt_buffer;
  int nidentical, i;

  if (altstart < altend) {
    exonlength = altend - altstart;
    mainstart = mainend - exonlength;
    main_buffer = (char *) CALLOC(exonlength+1,sizeof(char));
    alt_buffer = (char *) CALLOC(exonlength+1,sizeof(char));
    Genome_fill_buffer_simple(genome,chroffset+mainstart,exonlength,main_buffer);
    Genome_fill_buffer_simple(genome,chroffset+altstart,exonlength,alt_buffer);

  } else {
    exonlength = altstart - altend;
    mainstart = mainend + exonlength;
    main_buffer = (char *) CALLOC(exonlength+1,sizeof(char));
    alt_buffer = (char *) CALLOC(exonlength+1,sizeof(char));
    Genome_fill_buffer_simple(genome,chroffset+mainend,exonlength,main_buffer);
    Genome_fill_buffer_simple(genome,chroffset+altend,exonlength,alt_buffer);
  }
  
  nidentical = 0;
  for (i = 0; i < (int) exonlength; i++) {
    if (main_buffer[i] == alt_buffer[i]) {
      nidentical++;
    }
  }

  debug4(printf("Main: %s\n",main_buffer));
  debug4(printf("Alt: %s\n",alt_buffer));
  debug4(printf("nidentical = %d, ntotal = %d, pct = %f\n",
		nidentical,exonlength,(double) identical/(double) exonlength));

  FREE(alt_buffer);
  FREE(main_buffer);

  if ((double) nidentical/(double) exonlength > 0.90) {
    return true;
  } else {
    return false;
  }
}


static Acceptor_T
get_donor_main (Exon_T *mainexon, Score_T *mainscore, Donor_T site) {
  Acceptor_T mainprev = NULL, prev_site;
  Score_T score;
  Exon_T exon;
  List_T p;
  Genomicpos_T known_prev_acceptorpos;
  bool saw_knownp = false;

  if ((*mainexon = site->main_exon) != NULL) {
    /* Must have been set previously by known transcripts */
    prev_site = (*mainexon)->prev_site;
    *mainscore = prev_site->score + site->points + (*mainexon)->incrscore;
    return prev_site;

  } else {
    *mainexon = (Exon_T) NULL;
    *mainscore = 0;

    debug(printf("Looking for main at donor: %u (points %d)\n",site->chrpos,site->points));
    known_prev_acceptorpos = (Genomicpos_T) site->known_prev_acceptorpos;

    for (p = site->exons; p != NULL; p = List_next(p)) {
      exon = (Exon_T) List_head(p);
      debug(printf("  Evaluating acceptor at %u",exon->prev_site->chrpos));

      if (exon->prev_site->chrpos == known_prev_acceptorpos) {
	debug(printf(" ... known_main"));
	prev_site = exon->prev_site;

	*mainexon = exon;
	*mainscore = prev_site->score + site->points + exon->incrscore;
	mainprev = prev_site;
	saw_knownp = true;

      } else if (saw_knownp == true) {
	/* Cannot supplant KNOWN_MAIN */
	debug(printf(" ... not known_main"));

      } else if (exon->validp == true) {	/* exons are always valid */
	debug(printf(" ... valid"));

	prev_site = exon->prev_site;
	if (prev_site->prev != NULL) {
	  /* Cannot point to a acceptor that has no donor */
	  debug(printf(" ... has a prev"));
	
	  debug(printf(" ... score %d + points %d + exon %d = %d",
		       prev_site->score,site->points,exon->incrscore,prev_site->score+site->points+exon->incrscore));
	  if ((score = prev_site->score + site->points + exon->incrscore) > *mainscore) {
	    debug(printf("  **"));
	    *mainexon = exon;
	    *mainscore = score;
	    mainprev = prev_site;
	  }
	}

      }

      debug(printf("\n"));
    }

    debug(printf("\n"));
    return mainprev;
  }
}


static List_T
get_donor_alts (List_T *altexons, Donor_T site, Acceptor_T excludeprev) {
  List_T altprevs = NULL, p;
  Acceptor_T prev_site;
  Exon_T exon;

  altprevs = (List_T) NULL;
  *altexons = (List_T) NULL;

  debug(printf("Looking for alts at donor: %u\n",site->chrpos));
  for (p = site->exons; p != NULL; p = List_next(p)) {
    exon = (Exon_T) List_head(p);

    debug(printf("  Evaluating acceptor at %u",exon->prev_site->chrpos));

    if ((prev_site = exon->prev_site) == excludeprev) {
      debug(printf(" ... excluded"));

    } else if (prev_site->prev == NULL) {
      debug(printf(" ... cannot point to an acceptor that has no donor"));

    } else if (exon->validp == true) {	/* exons are always valid */
      debug(printf(" ... score %d + points %d + exon %d = %d",
		   prev_site->score,site->points,exon->incrscore,prev_site->score+site->points+exon->incrscore));
      altprevs = List_push(altprevs,(void *) prev_site);
      *altexons = List_push(*altexons,(void *) exon);
    }

    debug(printf("\n"));
  }

  debug(printf("\n"));
  return altprevs;
}


static void
process_donor_struct_obs (Donor_T site) {
  int mainscore;
  Acceptor_T mainprev;
  Exon_T mainexon;

  /* Main route */
  if ((mainprev = get_donor_main(&mainexon,&mainscore,site)) != NULL) {
    site->main_exon = mainexon;
    site->score = mainscore;
    site->prev = mainprev;
    if (mainprev->known_terminalp == true) {
      /* Don't alter known terminal */
    } else {
      mainprev->terminalp = false;
    }
  }

  return;
}

static void
process_donor_struct_known (Donor_T site) {
  int mainscore;
  Acceptor_T mainprev;
  Exon_T mainexon;

  /* Main route */
  if ((mainprev = get_donor_main(&mainexon,&mainscore,site)) != NULL) {
    site->main_exon = mainexon;
    site->score = mainscore;
    site->prev = mainprev;
    mainprev->terminalp = false;
  }

  return;
}


static bool
same_acceptor_chrpos (Donor_T altsite, Donor_T mainsite, int nlevels);


static bool
same_donor_chrpos (Acceptor_T altsite, Acceptor_T mainsite, int nlevels) {
  List_T p, q;
  Acceptor_alt_T alt1, alt2;

  if (nlevels == 0) {
    if (altsite->chrpos == mainsite->chrpos) {
      return true;
    } else {
      return false;
    }
  } else if (same_acceptor_chrpos(altsite->prev,mainsite->prev,nlevels-1) == true) {
    return true;
  } else {
    for (p = altsite->alt_structs; p != NULL; p = List_next(p)) {
      alt1 = (Acceptor_alt_T) List_head(p);
      if (same_acceptor_chrpos(alt1->prev_site,mainsite->prev,nlevels-1) == true) {
	return true;
      }
    }
    for (p = altsite->alt_sites; p != NULL; p = List_next(p)) {
      alt1 = (Acceptor_alt_T) List_head(p);
      if (same_acceptor_chrpos(alt1->prev_site,mainsite->prev,nlevels-1) == true) {
	return true;
      }
    }
    for (p = mainsite->alt_structs; p != NULL; p = List_next(p)) {
      alt2 = (Acceptor_alt_T) List_head(p);
      if (same_acceptor_chrpos(altsite->prev,alt2->prev_site,nlevels-1) == true) {
	return true;
      }
    }
    for (p = mainsite->alt_sites; p != NULL; p = List_next(p)) {
      alt2 = (Acceptor_alt_T) List_head(p);
      if (same_acceptor_chrpos(altsite->prev,alt2->prev_site,nlevels-1) == true) {
	return true;
      }
    }

    for (p = altsite->alt_structs; p != NULL; p = List_next(p)) {
      alt1 = (Acceptor_alt_T) List_head(p);
      for (q = altsite->alt_structs; q != NULL; q = List_next(q)) {
	alt2 = (Acceptor_alt_T) List_head(q);
	if (same_acceptor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
      for (q = altsite->alt_sites; q != NULL; q = List_next(q)) {
	alt2 = (Acceptor_alt_T) List_head(q);
	if (same_acceptor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
    }
    
    for (p = altsite->alt_sites; p != NULL; p = List_next(p)) {
      alt1 = (Acceptor_alt_T) List_head(p);
      for (q = altsite->alt_structs; q != NULL; q = List_next(q)) {
	alt2 = (Acceptor_alt_T) List_head(q);
	if (same_acceptor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
      for (q = altsite->alt_sites; q != NULL; q = List_next(q)) {
	alt2 = (Acceptor_alt_T) List_head(q);
	if (same_acceptor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
    }

    return false;
  }
}

static bool
same_acceptor_chrpos (Donor_T altsite, Donor_T mainsite, int nlevels) {
  List_T p, q;
  Donor_alt_T alt1, alt2;

  if (nlevels == 0) {
    if (altsite->chrpos == mainsite->chrpos) {
      return true;
    } else {
      return false;
    }
  } else if (same_donor_chrpos(altsite->prev,mainsite->prev,nlevels-1) == true) {
    return true;
  } else {
    for (p = altsite->alt_structs; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      if (same_donor_chrpos(alt1->prev_site,mainsite->prev,nlevels-1) == true) {
	return true;
      }
    }
    for (p = altsite->alt_sites; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      if (same_donor_chrpos(alt1->prev_site,mainsite->prev,nlevels-1) == true) {
	return true;
      }
    }
    for (p = mainsite->alt_structs; p != NULL; p = List_next(p)) {
      alt2 = (Donor_alt_T) List_head(p);
      if (same_donor_chrpos(altsite->prev,alt2->prev_site,nlevels-1) == true) {
	return true;
      }
    }
    for (p = mainsite->alt_sites; p != NULL; p = List_next(p)) {
      alt2 = (Donor_alt_T) List_head(p);
      if (same_donor_chrpos(altsite->prev,alt2->prev_site,nlevels-1) == true) {
	return true;
      }
    }

    for (p = altsite->alt_structs; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      for (q = altsite->alt_structs; q != NULL; q = List_next(q)) {
	alt2 = (Donor_alt_T) List_head(q);
	if (same_donor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
      for (q = altsite->alt_sites; q != NULL; q = List_next(q)) {
	alt2 = (Donor_alt_T) List_head(q);
	if (same_donor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
    }
    
    for (p = altsite->alt_sites; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      for (q = altsite->alt_structs; q != NULL; q = List_next(q)) {
	alt2 = (Donor_alt_T) List_head(q);
	if (same_donor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
      for (q = altsite->alt_sites; q != NULL; q = List_next(q)) {
	alt2 = (Donor_alt_T) List_head(q);
	if (same_donor_chrpos(alt1->prev_site,alt2->prev_site,nlevels-1) == true) {
	  return true;
	}
      }
    }

    return false;
  }
}



static void
process_donor_sites_obs (Donor_T site) {
  List_T altprevs, altexons, p, q;

  Acceptor_T mainprev, altprev;
  Intron_T mainintron, altintron;
  Exon_T mainexon, altexon;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainexon = site->main_exon;

    /* Alternate route for donor */
    altprevs = get_donor_alts(&altexons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altexons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Acceptor_T) List_head(p);
      altexon = (Exon_T) List_head(q);

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  /* altprev->chrpos != mainprev->chrpos && */
	  same_donor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
	  Donor_alt_equiv(site->alt_sites,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {
	if (altprev->known_terminalp == true) {
	  /* Don't alter known terminal */
	} else {
	  altprev->terminalp = false;
	}

	mainintron = Acceptor_find_intron(mainexon->prev_site);
	altintron = Acceptor_find_intron(altexon->prev_site);
	/* main_density = Splice_density(mainintron->splice); */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already known */
	} else if (altintron->introntype == NONE) {
	  /* Alternate acceptor splice site */
	  altintron->introntype = ALT_ACCEPTOR;
	  site->alt_sites = 
	    List_push(site->alt_sites,
		      (void *) Donor_alt_new(altprev,/*introntype*/ALT_ACCEPTOR,altexon,
					     mainintron,altintron,/*knownp*/false));
	}
      }
    }

    List_free(&altexons);
    List_free(&altprevs);
  }

  return;
}


static void
process_donor_sites_known (Donor_T site) {
  List_T altprevs, altexons, p, q;

  Acceptor_T mainprev, altprev;
  Exon_T mainexon, altexon;

#ifdef DEBUG5
  Intron_T mainintron, altintron;
#endif


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainexon = site->main_exon;

    debug5(printf("Donor is %u, with main acceptor at position %u\n",
		  site->chrpos,mainprev->chrpos));

    /* Alternate route for donor */
    altprevs = get_donor_alts(&altexons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altexons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Acceptor_T) List_head(p);
      altexon = (Exon_T) List_head(q);
      debug5(printf(" Considering alternate acceptor at position %u\n",altprev->chrpos));

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  /* altprev->chrpos != mainprev->chrpos && */
	  same_donor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
	  Donor_alt_equiv(site->alt_sites,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {
	if (altprev->known_terminalp == true) {
	  /* Don't alter known terminal */
	} else {
	  altprev->terminalp = false;
	}

	/* Alternate acceptor splice site */
	site->alt_sites = 
	  List_push(site->alt_sites,
		    (void *) Donor_alt_new(altprev,/*introntype*/ALT_ACCEPTOR,altexon,
					   /*mainintron*/Acceptor_find_intron(mainexon->prev_site),
					   /*altintron*/Acceptor_find_intron(altexon->prev_site),
					   /*knownp*/true));

	debug5(mainintron = Acceptor_find_intron(mainexon->prev_site));
	debug5(altintron = Acceptor_find_intron(altexon->prev_site));
	debug5(printf("Found alt_acceptor off main\n"));
	debug5(printf("mainintron: %u..%u\n",mainintron->prev_site->chrpos,mainintron->site->chrpos));
	debug5(printf("altintron: %u..%u\n",altintron->prev_site->chrpos,altintron->site->chrpos));
      }
    }

    List_free(&altexons);
    List_free(&altprevs);
  }

  return;
}


static void
process_donor_switch_obs (Donor_T site) {
  List_T altprevs, altexons, p, q;

  Acceptor_T mainprev, altprev;
  Exon_T mainexon, altexon;
  Intron_T mainintron, altintron;

  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainexon = site->main_exon;

    /* Alternate route for donor */
    altprevs = get_donor_alts(&altexons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altexons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Acceptor_T) List_head(p);
      altexon = (Exon_T) List_head(q);

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  altprev->prev->prev == NULL && mainprev->prev->prev == NULL &&
	  altprev->chrpos != mainprev->chrpos &&
	  same_donor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
	  Donor_alt_equiv(site->alt_switches,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {
	
	mainintron = Acceptor_find_intron(mainexon->prev_site);
	altintron = Acceptor_find_intron(altexon->prev_site);
	/* main_density = Splice_density(mainintron->splice); */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already known */
	} else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Switch intron at beginning of gene */
	  altintron->introntype = SWITCH_INTRON;
	  site->alt_switches =
	    List_push(site->alt_switches,
		      (void *) Donor_alt_new(altprev,/*introntype*/SWITCH_INTRON,altexon,
					     mainintron,altintron,/*knownp*/false));

	  /* altprev->terminalp = false;  ??? */
	}

      } else if (altprev->prev != NULL && mainprev->prev != NULL &&
		 altprev->prev->prev != NULL && mainprev->prev->prev != NULL &&
		 altprev->chrpos != mainprev->chrpos &&
		 same_donor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
		 same_acceptor_chrpos(altprev->prev,mainprev->prev,/*nlevels*/0) == true &&
		 Donor_alt_equiv(site->alt_switches,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {

	mainintron = Acceptor_find_intron(mainexon->prev_site);
	altintron = Acceptor_find_intron(altexon->prev_site);
	/* main_density = Splice_density(mainintron->splice); */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already known */
	} else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Switch intron in middle of gene */
	  altintron->introntype = SWITCH_INTRON;
	  site->alt_switches =
	    List_push(site->alt_switches,
		      (void *) Donor_alt_new(altprev,/*introntype*/SWITCH_INTRON,altexon,
					     mainintron,altintron,/*knownp*/false));
	  
	  /* altprev->terminalp = false;  ??? */
	}
	
      }
    }

    List_free(&altexons);
    List_free(&altprevs);
  }

  return;
}

static void
process_donor_switch_known (Donor_T site) {
  List_T altprevs, altexons, p, q;

  Acceptor_T mainprev, altprev;
  Exon_T mainexon, altexon;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainexon = site->main_exon;

    /* Alternate route for donor */
    altprevs = get_donor_alts(&altexons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altexons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Acceptor_T) List_head(p);
      altexon = (Exon_T) List_head(q);

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  altprev->prev->prev == NULL && mainprev->prev->prev == NULL &&
	  altprev->chrpos != mainprev->chrpos &&
	  same_donor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
	  Donor_alt_equiv(site->alt_switches,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {

	/* Switch intron at beginning of gene */
	site->alt_switches =
	  List_push(site->alt_switches,
		    (void *) Donor_alt_new(altprev,/*introntype*/SWITCH_INTRON,altexon,
					   /*mainintron*/Acceptor_find_intron(mainexon->prev_site),
					   /*altintron*/Acceptor_find_intron(altexon->prev_site),
					   /*knownp*/true));
	
      } else if (altprev->prev != NULL && mainprev->prev != NULL &&
		 altprev->prev->prev != NULL && mainprev->prev->prev != NULL &&
		 altprev->chrpos != mainprev->chrpos &&
		 same_donor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
		 same_acceptor_chrpos(altprev->prev,mainprev->prev,/*nlevels*/0) == true &&
		 Donor_alt_equiv(site->alt_switches,/*altintron*/Acceptor_find_intron(altexon->prev_site)) == false) {

	/* Switch intron in middle of gene */
	site->alt_switches =
	  List_push(site->alt_switches,
		    (void *) Donor_alt_new(altprev,/*introntype*/SWITCH_INTRON,altexon,
					   /*mainintron*/Acceptor_find_intron(mainexon->prev_site),
					   /*altintron*/Acceptor_find_intron(altexon->prev_site),
					   /*knownp*/true));
      }
    }

    List_free(&altexons);
    List_free(&altprevs);
  }

  return;
}


static Donor_T
get_acceptor_main (Intron_T *mainintron, Score_T *mainscore, Acceptor_T site) {
  Donor_T mainprev = NULL, prev_site;
  Score_T score;
  Intron_T intron;
  List_T p;
  bool saw_knownp = false;

  if ((*mainintron = site->main_intron) != NULL) {
    /* Must have been set previously by known transcripts */
    prev_site = (*mainintron)->prev_site;
    *mainscore = prev_site->score + site->points + Splice_count((*mainintron)->splice);
    return prev_site;

  } else {
    *mainintron = (Intron_T) NULL;
    *mainscore = 0;

    debug(printf("Looking for main at acceptor: %u (points %d)\n",site->chrpos,site->points));
    for (p = site->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);

      debug(printf("  Evaluating donor at %u",intron->prev_site->chrpos));
      if (Splice_knownp(intron->splice) == KNOWN_MAIN) {
	debug(printf(" ... known_main"));
	prev_site = intron->prev_site;

	*mainintron = intron;
	*mainscore = prev_site->score + site->points + Splice_count(intron->splice);
	mainprev = prev_site;
	saw_knownp = true;

      } else if (saw_knownp == true) {
	/* Cannot supplant KNOWN_MAIN */
	debug(printf(" ... not known_main"));

      } else if (intron->validp == true) {
	debug(printf(" ... valid"));
      
	prev_site = intron->prev_site;
	debug(printf(" ... score %d + points %d + intron %d (crosshyb %d) = %d",
		     prev_site->score,site->points,Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),
		     prev_site->score+site->points+Splice_count(intron->splice)));
	if ((score = prev_site->score + site->points + Splice_count(intron->splice)) > *mainscore) {
	  debug(printf("  **"));
	  *mainintron = intron;
	  *mainscore = score;
	  mainprev = prev_site;
	}
      }

      debug(printf("\n"));
    }

    debug(printf("\n"));
    return mainprev;
  }
}


static List_T
get_acceptor_alts (List_T *altintrons, Acceptor_T site, Donor_T excludeprev) {
  List_T altprevs = NULL, p;
  Donor_T prev_site;
  Intron_T intron;


  altprevs = (List_T) NULL;
  *altintrons = (List_T) NULL;

  debug(printf("Looking for alts at acceptor: %u\n",site->chrpos));
  for (p = site->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);

    debug(printf("  Evaluating donor at %u",intron->prev_site->chrpos));

    if ((prev_site = intron->prev_site) == excludeprev) {
      debug(printf(" ... excluded"));

    } else if (intron->validp == true) {
      debug(printf(" ... score %d + points %d + intron %d (crosshyb %d) = %d",
		   prev_site->score,site->points,Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),
		   prev_site->score+site->points+Splice_count(intron->splice)));
      altprevs = List_push(altprevs,(void *) prev_site);
      *altintrons = List_push(*altintrons,(void *) intron);
    }

    debug(printf("\n"));
  }

  debug(printf("\n"));
  return altprevs;
}


static void
process_acceptor_struct_obs (Acceptor_T site, long int *tally_matches, long int *tally_mismatches,
			     Genomicpos_T chrlength, Genome_T genome, Genomicpos_T chroffset,
			     bool forwardp) {
  List_T altprevs, altintrons, p, q;

  Donor_T mainprev, altprev;
  Intron_T mainintron, mainintron2, mainintron3, altintron, altintron2;
  Score_T mainscore;
  Genomicpos_T genestart1, genestart2;


  /* Main route */
  if ((mainprev = get_acceptor_main(&mainintron,&mainscore,site)) != NULL) {
    mainintron->introntype = MAIN;
    site->main_intron = mainintron;
    site->score = mainscore;
    site->prev = mainprev;

    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

#if 0
      if (altprev->prev != NULL) {
	printf("Considering switch exon at %u..%u\n",altprev->chrpos,altprev->prev->chrpos);
      }
#endif

      if (mainprev->prev != NULL && mainprev->prev->prev != NULL && 
	  mainprev->prev->prev->chrpos == altprev->chrpos &&
	  Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {

	mainintron2 = Acceptor_find_intron(mainprev->prev);
	/* main_density = (Splice_density(mainintron->splice) + Splice_density(mainintron2->splice))/2.0; */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already a known skipped exon */
	} else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Skipped exon */
	  altintron->introntype = SKIPPED_EXON;
	  site->alt_structs =
	    List_push(site->alt_structs,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/SKIPPED_EXON,
						mainintron,mainintron2,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/false));

	  /* site->alt_incrscore = mainincrscore; */
	  /* altprev->terminalp = false; */
	}

      } else if (mainprev->prev != NULL && mainprev->prev->prev != NULL &&
		 mainprev->prev->prev->prev != NULL && mainprev->prev->prev->prev->prev != NULL &&
		 mainprev->prev->prev->prev->prev->chrpos == altprev->chrpos &&
		 Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {

	mainintron2 = Acceptor_find_intron(mainprev->prev);
	mainintron3 = Acceptor_find_intron(mainprev->prev->prev->prev);
	/* main_density = (Splice_density(mainintron->splice) + Splice_density(mainintron2->splice) + Splice_density(mainintron3->splice))/3.0; */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already a known skipped double exon */
	} else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Skipped double exon */
	  altintron->introntype = SKIPPED_DOUBLE_EXON;
	  site->alt_structs =
	    List_push(site->alt_structs,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/SKIPPED_DOUBLE_EXON,
						mainintron,mainintron2,mainintron,
						altintron,/*altintron2*/NULL,/*knownp*/false));

	  /* site->alt_incrscore = mainincrscore; */
	  /* altprev->terminalp = false; */
	}

      } else if (altprev->prev != NULL && altprev->prev->prev != NULL &&
		 altprev->prev->prev->chrpos == mainprev->chrpos &&
		 Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/Acceptor_find_intron(altprev->prev)) == false) {

	altintron2 = Acceptor_find_intron(altprev->prev);
	/* main_density = Splice_density(mainintron->splice); */
	/* alt_density = (Splice_density(altintron->splice) + Splice_density(altintron2->splice))/2.0; */

	if (Splice_known_count(altintron->splice) > 0 && Splice_known_count(altintron2->splice) > 0) {
	  /* Already a known extra exon */
	} else if (altintron->introntype == NONE && altintron2->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Extra exon */
	  altintron->introntype = EXTRA_EXON;
	  altintron2->introntype = EXTRA_EXON;
	  site->alt_structs =
	    List_push(site->alt_structs,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/EXTRA_EXON,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,altintron2,/*knownp*/false));

	  /* site->alt_incrscore = mainincrscore; */
	  /* altprev->terminalp = false; */
	}

      } else if (altprev->prev == NULL) {
	genestart1 = Cappaths_solve_genestart(mainprev->chrpos,tally_matches,tally_mismatches,
					      /*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);
	genestart2 = Cappaths_solve_genestart(altprev->chrpos,tally_matches,tally_mismatches,
					      /*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);

	if (mainprev->prev != NULL &&
	    Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {

	  /* main_density = Splice_density(mainintron->splice); */
	  /* alt_density = Splice_density(altintron->splice); */

	  if (Splice_known_count(altintron->splice) > 0) {
	    /* Already a known initial short */
	  } else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	    if (homologousp(genome,chroffset,/*mainend*/site->prev->chrpos,
			    /*altstart*/genestart2,/*altend*/altprev->chrpos) == true) {
	      debug4(printf("Found homology at %u..%u with %u\n",genestart2,altprev->chrpos,site->prev->chrpos));
	    } else {
	      /* Alternate initial short path */
	      altintron->introntype = INITIAL_SHORT;
	      site->alt_structs =
		List_push(site->alt_structs,
			  (void *) Acceptor_alt_new(altprev,/*introntype*/INITIAL_SHORT,
						    mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						    altintron,/*altintron2*/NULL,/*knownp*/false));
	    }
	  }
	} else if (genestart1 != genestart2 &&
		   Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {

	  /* main_density = Splice_density(mainintron->splice); */
	  /* alt_density = Splice_density(altintron->splice); */

	  if (Splice_known_count(altintron->splice) > 0) {
	    /* Already known */
	  } else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	    /* Alternate first exon */
	    altintron->introntype = INITIAL_ALT_EXON;
	    site->alt_structs =
	      List_push(site->alt_structs,
			(void *) Acceptor_alt_new(altprev,/*introntype*/INITIAL_ALT_EXON,
						  mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						  altintron,/*altintron2*/NULL,/*knownp*/false));
	  }
	}
      }
    }

    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}


/* Need tally_matches to compute genestart */
static void
process_acceptor_struct_known (Acceptor_T site, IIT_T end_exons_iit, char *chr,
			       bool forwardp) {
  List_T altprevs, altintrons, p, q;
  Score_T mainscore;

  Donor_T mainprev, altprev;
  Intron_T mainintron, altintron;
  Genomicpos_T genestart1, genestart2;


  /* Main route */
  if ((mainprev = get_acceptor_main(&mainintron,&mainscore,site)) != NULL) {
    mainintron->introntype = MAIN;
    site->main_intron = mainintron;
    site->score = mainscore;
    site->prev = mainprev;


    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

      if (mainprev->prev != NULL && mainprev->prev->prev != NULL && 
	  mainprev->prev->prev->chrpos == altprev->chrpos &&
	  Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {
	/* Skipped exon */
	site->alt_structs =
	  List_push(site->alt_structs,
		    (void *) Acceptor_alt_new(altprev,/*introntype*/SKIPPED_EXON,
					      mainintron,/*mainintron2*/Acceptor_find_intron(mainprev->prev),
					      /*mainintron3*/NULL,altintron,/*altintron2*/NULL,/*knownp*/true));
	
      } else if (mainprev->prev != NULL && mainprev->prev->prev != NULL &&
		 mainprev->prev->prev->prev != NULL && mainprev->prev->prev->prev->prev != NULL &&
		 mainprev->prev->prev->prev->prev->chrpos == altprev->chrpos &&
		 Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {
	/* Skipped double exon */
	site->alt_structs =
	  List_push(site->alt_structs,
		    (void *) Acceptor_alt_new(altprev,/*introntype*/SKIPPED_DOUBLE_EXON,
					      mainintron,/*mainintron2*/Acceptor_find_intron(mainprev->prev),
					      /*mainintron3*/Acceptor_find_intron(mainprev->prev->prev->prev),
					      altintron,/*altintron2*/NULL,/*knownp*/true));

      } else if (altprev->prev != NULL && altprev->prev->prev != NULL &&
		 altprev->prev->prev->chrpos == mainprev->chrpos &&
		 Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/Acceptor_find_intron(altprev->prev)) == false) {
	/* Extra exon */
	site->alt_structs =
	  List_push(site->alt_structs,
		    (void *) Acceptor_alt_new(altprev,/*introntype*/EXTRA_EXON,
					      mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
					      altintron,/*altintron2*/Acceptor_find_intron(altprev->prev),
					      /*knownp*/true));

      } else if (altprev->prev == NULL &&
		 Acceptor_alt_equiv(site->alt_structs,altintron,/*altintron2*/NULL) == false) {
	genestart1 = Cappaths_solve_genestart(mainprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
					      end_exons_iit,chr,/*chrlength*/0,forwardp);
	genestart2 = Cappaths_solve_genestart(altprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
					      end_exons_iit,chr,/*chrlength*/0,forwardp);

	if (mainprev->prev != NULL) {
	  /* Alternate initial short path */
	  site->alt_structs =
	    List_push(site->alt_structs,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/INITIAL_SHORT,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/true));

	} else if (genestart1 != genestart2) {
	  /* Alternate first exon */
	  site->alt_structs =
	    List_push(site->alt_structs,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/INITIAL_ALT_EXON,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/true));
	}
      }
    }

    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}




static void
process_acceptor_sites_obs (Acceptor_T site, long int *tally_matches, long int *tally_mismatches,
			    Genomicpos_T chrlength, Genome_T genome, Genomicpos_T chroffset,
			    bool forwardp) {
  List_T altprevs, altintrons, p, q, a;

  Donor_T mainprev, altprev;
  Intron_T mainintron, altintron;
  Genomicpos_T genestart1, genestart2;

  Acceptor_alt_T alt_struct;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainintron = site->main_intron;

    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  /* altprev->chrpos != mainprev->chrpos && */
	  same_acceptor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
	  Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	/* Check distance_diff to avoid calling a switch exon as an alternate donor */
	
	/* main_density = Splice_density(mainintron->splice); */
	/* alt_density = Splice_density(altintron->splice); */

	if (Splice_known_count(altintron->splice) > 0) {
	  /* Already known */
	} else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Alternate splice site (donor) */
	  altintron->introntype = ALT_DONOR;
	  site->alt_sites =
	    List_push(site->alt_sites,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/false));
	}

      } else if (altprev->prev == NULL &&
		 Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	genestart1 = Cappaths_solve_genestart(mainprev->chrpos,tally_matches,tally_mismatches,
					      /*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);
	genestart2 = Cappaths_solve_genestart(altprev->chrpos,tally_matches,tally_mismatches,
					      /*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);

	if (mainprev->prev == NULL && genestart1 == genestart2) {
	  /* main_density = Splice_density(mainintron->splice); */
	  /* alt_density = Splice_density(altintron->splice); */

	  if (Splice_known_count(altintron->splice) > 0) {
	    /* Already known */
	  } else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	    /* Alternate splice site (donor) */
	    altintron->introntype = ALT_DONOR;
	    site->alt_sites =
	      List_push(site->alt_sites,
			(void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						  mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						  altintron,/*altintron2*/NULL,/*knownp*/false));
	  }
	}
      }
    }

    for (a = site->alt_structs; a != NULL; a = List_next(a)) {
      alt_struct = (Acceptor_alt_T) List_head(a);
      mainprev = alt_struct->prev_site;
      mainintron = alt_struct->alt_intron;

      for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
	altprev = (Donor_T) List_head(p);
	altintron = (Intron_T) List_head(q);

	if (altprev->chrpos == mainprev->chrpos) {
	  /* Skip */
	} else if (altprev->prev != NULL && mainprev->prev != NULL &&
		   altprev->chrpos != mainprev->chrpos &&
		   same_acceptor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
		   Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	
	  /* main_density = Splice_density(mainintron->splice); */
	  /* alt_density = Splice_density(altintron->splice); */

	  if (Splice_known_count(altintron->splice) > 0) {
	    /* Already known */
	  } else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	    /* Alternate splice site (donor) */
	    altintron->introntype = ALT_DONOR;
	    site->alt_sites =
	      List_push(site->alt_sites,
			(void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						  mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						  altintron,/*altintron2*/NULL,/*knownp*/false));
	  }

	} else if (altprev->prev == NULL &&
		   Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	  genestart1 = Cappaths_solve_genestart(mainprev->chrpos,tally_matches,tally_mismatches,
						/*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);
	  genestart2 = Cappaths_solve_genestart(altprev->chrpos,tally_matches,tally_mismatches,
						/*end_exons_iit*/NULL,/*chr*/NULL,chrlength,forwardp);

	  if (mainprev->prev == NULL && genestart1 == genestart2) {
	    /* main_density = Splice_density(mainintron->splice); */
	    /* alt_density = Splice_density(altintron->splice); */

	    if (Splice_known_count(altintron->splice) > 0) {
	      /* Already known */
	    } else if (altintron->introntype == NONE /* alt_density > 0.05 * main_density */) {
	      /* Alternate splice site (donor) */
	      altintron->introntype = ALT_DONOR;
	      site->alt_sites =
		List_push(site->alt_sites,
			  (void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						    mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						    altintron,/*altintron2*/NULL,/*knownp*/false));
	    }
	  }
	}
      }
    }

    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}


static void
process_acceptor_sites_known (Acceptor_T site, IIT_T end_exons_iit, char *chr,
			      bool forwardp) {
  List_T altprevs, altintrons, p, q, a;

  Donor_T mainprev, altprev;
  Intron_T mainintron, altintron;
  Genomicpos_T genestart1, genestart2;

  Acceptor_alt_T alt_struct;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainintron = site->main_intron;

    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

      if (altprev->prev != NULL && mainprev->prev != NULL &&
	  /* altprev->chrpos != mainprev->chrpos && */
	  same_acceptor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
	  Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	
	debug5(printf("Found alt_donor off main, middle\n"));
	debug5(printf("mainintron: %u..%u\n",mainintron->prev_site->chrpos,mainintron->site->chrpos));
	debug5(printf("altintron: %u..%u\n",altintron->prev_site->chrpos,altintron->site->chrpos));
	site->alt_sites =
	  List_push(site->alt_sites,
		    (void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
					      mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
					      altintron,/*altintron2*/NULL,/*knownp*/true));

      } else if (altprev->prev == NULL &&
		 Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	genestart1 = Cappaths_solve_genestart(mainprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
					      end_exons_iit,chr,/*chrlength*/0,forwardp);
	genestart2 = Cappaths_solve_genestart(altprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
					      end_exons_iit,chr,/*chrlength*/0,forwardp);

	if (mainprev->prev == NULL && genestart1 == genestart2) {
	  /* Alternate splice site (donor) */
	  debug5(printf("Found alt_donor off main, first\n"));
	  debug5(printf("mainintron: %u..%u\n",mainintron->prev_site->chrpos,mainintron->site->chrpos));
	  debug5(printf("altintron: %u..%u\n",altintron->prev_site->chrpos,altintron->site->chrpos));

	  site->alt_sites =
	    List_push(site->alt_sites,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/true));
	}
      }
    }

    for (a = site->alt_structs; a != NULL; a = List_next(a)) {
      alt_struct = (Acceptor_alt_T) List_head(a);
      mainprev = alt_struct->prev_site;
      mainintron = alt_struct->alt_intron;

      for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
	altprev = (Donor_T) List_head(p);
	altintron = (Intron_T) List_head(q);

	if (altprev->chrpos == mainprev->chrpos) {
	  /* Skip */
	} else if (altprev->prev != NULL && mainprev->prev != NULL &&
		   altprev->chrpos != mainprev->chrpos &&
		   same_acceptor_chrpos(altprev,mainprev,/*nlevels*/0) == true &&
		   Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	
	  debug5(printf("Found alt_donor off alt, middle\n"));
	  debug5(printf("mainintron: %u..%u\n",mainintron->prev_site->chrpos,mainintron->site->chrpos));
	  debug5(printf("altintron: %u..%u\n",altintron->prev_site->chrpos,altintron->site->chrpos));
	  site->alt_sites =
	    List_push(site->alt_sites,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						altintron,/*altintron2*/NULL,/*knownp*/true));

	} else if (altprev->prev == NULL &&
		   Acceptor_alt_equiv(site->alt_sites,altintron,/*altintron2*/NULL) == false) {
	  genestart1 = Cappaths_solve_genestart(mainprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
						end_exons_iit,chr,/*chrlength*/0,forwardp);
	  genestart2 = Cappaths_solve_genestart(altprev->chrpos,/*tally_matches*/NULL,/*tally_mismatches*/NULL,
						end_exons_iit,chr,/*chrlength*/0,forwardp);

	  if (mainprev->prev == NULL && genestart1 == genestart2) {
	    /* Alternate splice site (donor) */
	    debug5(printf("Found alt_donor off alt, first\n"));
	    debug5(printf("mainintron: %u..%u\n",mainintron->prev_site->chrpos,mainintron->site->chrpos));
	    debug5(printf("altintron: %u..%u\n",altintron->prev_site->chrpos,altintron->site->chrpos));
	    site->alt_sites =
	      List_push(site->alt_sites,
			(void *) Acceptor_alt_new(altprev,/*introntype*/ALT_DONOR,
						  mainintron,/*mainintron2*/NULL,/*mainintron3*/NULL,
						  altintron,/*altintron2*/NULL,/*knownp*/true));
	  }
	}
      }
    }


    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}


static void
process_acceptor_switch_obs (Acceptor_T site, long int *tally_matches, long int *tally_mismatches,
			     Genomicpos_T chrlength, Genome_T genome, Genomicpos_T chroffset,
			     bool forwardp) {
  List_T altprevs, altintrons, p, q;

  Donor_T mainprev, altprev;
  Intron_T mainintron, mainintron2, altintron, altintron2;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainintron = site->main_intron;

    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

      if (altprev->prev != NULL && altprev->prev->prev != NULL &&
	  mainprev->prev != NULL && mainprev->prev->prev != NULL &&
	  altprev->chrpos != mainprev->chrpos &&
	  same_acceptor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
	  same_donor_chrpos(altprev->prev,mainprev->prev,/*nlevels*/0) == true &&
	  Acceptor_alt_equiv(site->alt_switches,altintron,/*altintron2*/Acceptor_find_intron(altprev->prev)) == false) {
	/* printf("Considering switch exon at %u..%u\n",altprev->chrpos,altprev->prev->chrpos); */

	mainintron2 = Acceptor_find_intron(mainprev->prev);
	altintron2 = Acceptor_find_intron(altprev->prev);
	/* main_density = (Splice_density(mainintron->splice) + Splice_density(mainintron2->splice))/2.0; */
	/* alt_density = (Splice_density(altintron->splice) + Splice_density(altintron2->splice))/2.0; */

	if (Splice_known_count(altintron->splice) > 0 && Splice_known_count(altintron2->splice) > 0) {
	  /* Already known */
	} else if (altintron->introntype == NONE && altintron2->introntype == NONE /* alt_density > 0.05 * main_density */) {
	  /* Alternate mutually exclusive exon */
	  altintron->introntype = SWITCH_EXON;
	  site->alt_switches =
	    List_push(site->alt_switches,
		      (void *) Acceptor_alt_new(altprev,/*introntype*/SWITCH_EXON,
						mainintron,mainintron2,/*mainintron3*/NULL,
						altintron,altintron2,/*knownp*/false));
	}
      }
    }

    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}


/* Need tally_matches to compute genestart */
static void
process_acceptor_switch_known (Acceptor_T site, IIT_T end_exons_iit, char *chr,
			       bool forwardp) {
  List_T altprevs, altintrons, p, q;

  Donor_T mainprev, altprev;
  Intron_T mainintron, altintron;


  /* Main route */
  if ((mainprev = site->prev) != NULL) {
    mainintron = site->main_intron;
    site->prev = mainprev;


    /* Alternate route for acceptor */
    altprevs = get_acceptor_alts(&altintrons,site,/*excludeprev*/mainprev);
    for (p = altprevs, q = altintrons; p != NULL; p = List_next(p), q = List_next(q)) {
      altprev = (Donor_T) List_head(p);
      altintron = (Intron_T) List_head(q);

      /* For same_donor_chrpos, allowing nlevels of 1 gets cases where we have every-other skipped exons (e.g., 1,3,5 out of 1,2,3,4,5) */
      if (altprev->prev != NULL && altprev->prev->prev != NULL &&
	  mainprev->prev != NULL && mainprev->prev->prev != NULL &&
	  altprev->chrpos != mainprev->chrpos &&
	  same_acceptor_chrpos(altprev,mainprev,/*nlevels*/1) == false &&
	  same_donor_chrpos(altprev->prev,mainprev->prev,/*nlevels*/0) == true &&
	  Acceptor_alt_equiv(site->alt_switches,altintron,/*altintron2*/Acceptor_find_intron(altprev->prev)) == false) {
#if 0
	printf("Considering switch exon:\n");
	printf("  mainprev %u, mainprev->prev %u, mainprev->prev->prev %u\n",mainprev->chrpos,mainprev->prev->chrpos,mainprev->prev->prev->chrpos);
	printf("  altprev %u, altprev->prev %u, altprev->prev->prev %u\n",altprev->chrpos,altprev->prev->chrpos,altprev->prev->prev->chrpos);
#endif

	/* Alternate mutually exclusive exon */
	site->alt_switches =
	  List_push(site->alt_switches,
		    (void *) Acceptor_alt_new(altprev,/*introntype*/SWITCH_EXON,
					      mainintron,/*mainintron2*/Acceptor_find_intron(mainprev->prev),
					      /*mainintron3*/NULL,altintron,/*altintron2*/Acceptor_find_intron(altprev->prev),
					      /*knownp*/true));
      }
    }

    List_free(&altintrons);
    List_free(&altprevs);
  }

  return;
}


static void
traverse_graph_fwd (Donor_T *sites_donor_fwd, int nsites_donor_fwd,
		    Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd,
		    long int *tally_matches_low, long int *tally_mismatches_low,
		    Genome_T genome, char *chr, Genomicpos_T chroffset,
		    Genomicpos_T chrlength) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos;
  int donor_fwd_i, acceptor_fwd_i;

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_struct_obs(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_struct_obs(sites_acceptor_fwd[acceptor_fwd_i],tally_matches_low,tally_mismatches_low,
				  chrlength,genome,chroffset,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_sites_obs(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_sites_obs(sites_acceptor_fwd[acceptor_fwd_i],tally_matches_low,tally_mismatches_low,
				 chrlength,genome,chroffset,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_switch_obs(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_switch_obs(sites_acceptor_fwd[acceptor_fwd_i],tally_matches_low,tally_mismatches_low,
				  chrlength,genome,chroffset,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  return;
}

static void
traverse_graph_fwd_known (Donor_T *sites_donor_fwd, int nsites_donor_fwd,
			  Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd,
			  IIT_T end_exons_iit, char *chr) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos;
  int donor_fwd_i = 0, acceptor_fwd_i = 0;

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_struct_known(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_struct_known(sites_acceptor_fwd[acceptor_fwd_i],end_exons_iit,chr,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_sites_known(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_sites_known(sites_acceptor_fwd[acceptor_fwd_i],end_exons_iit,chr,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  donor_fwd_i = acceptor_fwd_i = 0;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    if (donor_fwd_pos <= acceptor_fwd_pos) {

      process_donor_switch_known(sites_donor_fwd[donor_fwd_i]);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {

      process_acceptor_switch_known(sites_acceptor_fwd[acceptor_fwd_i],end_exons_iit,chr,/*forwardp*/true);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  return;
}


/************************************************************************/


static void
traverse_graph_rev (Donor_T *sites_donor_rev, int nsites_donor_rev,
		    Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev,
		    long int *tally_matches_high, long int *tally_mismatches_high,
		    Genome_T genome, char *chr, Genomicpos_T chroffset,
		    Genomicpos_T chrlength) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos;
  int donor_rev_i, acceptor_rev_i;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_struct_obs(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_struct_obs(sites_acceptor_rev[acceptor_rev_i],tally_matches_high,tally_mismatches_high,
				  chrlength,genome,chroffset,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_sites_obs(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_sites_obs(sites_acceptor_rev[acceptor_rev_i],tally_matches_high,tally_mismatches_high,
				 chrlength,genome,chroffset,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_switch_obs(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_switch_obs(sites_acceptor_rev[acceptor_rev_i],tally_matches_high,tally_mismatches_high,
				  chrlength,genome,chroffset,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  return;
}


static void
traverse_graph_rev_known (Donor_T *sites_donor_rev, int nsites_donor_rev,
			  Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev,
			  IIT_T end_exons_iit, char *chr) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos;
  int donor_rev_i, acceptor_rev_i;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_struct_known(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_struct_known(sites_acceptor_rev[acceptor_rev_i],end_exons_iit,chr,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_sites_known(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_sites_known(sites_acceptor_rev[acceptor_rev_i],end_exons_iit,chr,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  donor_rev_i = acceptor_rev_i = 0;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    if (donor_rev_pos >= acceptor_rev_pos) {

      process_donor_switch_known(sites_donor_rev[donor_rev_i]);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {

      process_acceptor_switch_known(sites_acceptor_rev[acceptor_rev_i],end_exons_iit,chr,/*forwardp*/false);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  return;
}



/************************************************************************
 *   Paths
 ************************************************************************/

struct Path_T {
  bool knownp;
  Introntype_T type;
  bool forwardp;

  double density;		/* Computed for mainpath only */
  List_T main_introns;
  List_T alt_introns;

  bool primaryp;
  bool altp;

  Genomicpos_T genestart;	/* Start of first intron */
  Genomicpos_T geneend;		/* End of last intron */
};


static void
Path_free (Path_T *old) {

#if 0
  for (p = (*old)->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    Intron_free(&intron);
  }
#endif
  List_free(&(*old)->alt_introns);
  List_free(&(*old)->main_introns);
  FREE(*old);
  return;
}

static int
Path_fwd_cmp (const void *a, const void *b) {
  Path_T x = * (Path_T *) a;
  Path_T y = * (Path_T *) b;

  if (x->knownp == true && y->knownp == false) {
    return -1;
  } else if (y->knownp == true && x->knownp == false) {
    return +1;
  } else if (x->genestart < y->genestart) {
    return -1;
  } else if (y->genestart < x->genestart) {
    return +1;
  } else if (x->geneend < y->geneend) {
    return -1;
  } else if (y->geneend < x->geneend) {
    return +1;
  } else {
    return 0;
  }
}

static int
Path_rev_cmp (const void *a, const void *b) {
  Path_T x = * (Path_T *) a;
  Path_T y = * (Path_T *) b;

  if (x->knownp == true && y->knownp == false) {
    return -1;
  } else if (y->knownp == true && x->knownp == false) {
    return +1;
  } else if (x->genestart > y->genestart) {
    return -1;
  } else if (y->genestart > x->genestart) {
    return +1;
  } else if (x->geneend > y->geneend) {
    return -1;
  } else if (y->geneend > x->geneend) {
    return +1;
  } else {
    return 0;
  }
}


static Path_T
Path_new_main (List_T introns, Introntype_T type, bool forwardp, bool altp, bool knownp) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));

  new->knownp = knownp;
  new->type = type;

  new->main_introns = introns;
  new->alt_introns = (List_T) NULL;

  new->primaryp = true;
  new->forwardp = forwardp;
  new->altp = altp;
  return new;
}
 

static Path_T
Path_donor_alt (Donor_alt_T alt, bool forwardp) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));

  new->knownp = alt->knownp;
  new->type = alt->type;

  new->main_introns = List_push(NULL,alt->main_intron);
  new->alt_introns = List_push(NULL,alt->alt_intron);

  new->primaryp = true;
  new->forwardp = forwardp;
  new->altp = true;
  return new;
}

static Path_T
Path_acceptor_alt (Acceptor_alt_T alt, bool forwardp) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));

  new->knownp = alt->knownp;
  new->type = alt->type;

  new->main_introns = List_push(NULL,alt->main_intron);
  if (alt->main_intron_2 != NULL) {
    new->main_introns = List_push(new->main_introns,alt->main_intron_2);
  }
  if (alt->main_intron_3 != NULL) {
    new->main_introns = List_push(new->main_introns,alt->main_intron_3);
  }

  new->alt_introns = List_push(NULL,alt->alt_intron);
  if (alt->alt_intron_2 != NULL) {
    new->alt_introns = List_push(new->alt_introns,alt->alt_intron_2);
  }

  new->primaryp = true;
  new->forwardp = forwardp;
  new->altp = true;
  return new;
}


static Path_T
Path_new_retained_intron (Intron_T main_intron, bool forwardp, bool knownp) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));

  new->knownp = knownp;
  new->type = RETAINED_INTRON;

  new->main_introns = List_push(NULL,main_intron);
  new->alt_introns = List_push(NULL,main_intron); /* Deleted after computing genebounds */

  new->primaryp = true;
  new->forwardp = forwardp;
  new->altp = true;
  return new;
}



static Path_T
Path_copy (Path_T orig) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));

  new->knownp = orig->knownp;
  new->type = orig->type;
  new->forwardp = orig->forwardp;

  new->main_introns = List_copy(orig->main_introns);
  new->alt_introns = List_copy(orig->alt_introns);

  new->primaryp = orig->primaryp;
  new->altp = orig->altp;

  new->genestart = orig->genestart;
  new->geneend = orig->geneend;

  return new;
}


static bool
Path_primaryp (Path_T path, char *chr, long int *primary_extents, long int *crosshyb_extents) {
  bool primaryp;
  Genomicpos_T chrstart, chrend, chrpos;

  if (crosshyb_extents == NULL) {
    return true;

  } else if (primary_extents == NULL) {
    return false;

  } else {
    primaryp = false;
    if (path->genestart < path->geneend) {
      chrstart = path->genestart;
      chrend = path->geneend;
    } else {
      chrstart = path->geneend;
      chrend = path->genestart;
    }

    for (chrpos = chrstart; chrpos <= chrend; chrpos++) {
      if (primary_extents[chrpos] >= 20 &&
	  (double) primary_extents[chrpos]/(double) (primary_extents[chrpos] + crosshyb_extents[chrpos]) > 0.80) {
	primaryp = true;
      }
    }

    return primaryp;
  }
}


#if 0
static long int
compute_sum (long int *maxcount, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend) {
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



static long int
tally_coverage (int *npositions, Genomicpos_T chrstart, Genomicpos_T chrend,
		long int *tally_matches_low, long int *tally_mismatches_low,
		long int *tally_matches_high, long int *tally_mismatches_high) {
  long int sum = 0, maxcount;
  Genomicpos_T low, high;

  if (chrstart < chrend) {
    low = chrstart + 1U;
    high = chrend - 1U;
  } else {
    low = chrend + 1U;
    high = chrstart - 1U;
  }
  
  sum = compute_sum(&maxcount,tally_matches_low,low,high);
  sum += compute_sum(&maxcount,tally_mismatches_low,low,high);
  sum += compute_sum(&maxcount,tally_matches_high,low,high);
  sum += compute_sum(&maxcount,tally_mismatches_high,low,high);

  *npositions = high - low + 1U;
  return sum;
}
#endif


static long int
tally_mincount (Genomicpos_T chrstart, Genomicpos_T chrend,
		long int *tally_matches_low, long int *tally_mismatches_low,
		long int *tally_matches_high, long int *tally_mismatches_high) {
  long int mincount, count;
  Genomicpos_T low, high, pos;

  if (chrstart < chrend) {
    low = chrstart + 1U;
    high = chrend - 1U;
  } else {
    low = chrend + 1U;
    high = chrstart - 1U;
  }
  
  mincount = tally_matches_low[low] + tally_mismatches_low[low] +
    tally_matches_high[low] + tally_mismatches_high[low];

  for (pos = low + 1U; pos <= high; pos++) {
    if ((count = tally_matches_low[pos] + tally_mismatches_low[pos] +
	 tally_matches_high[pos] + tally_mismatches_high[pos]) < mincount) {
      mincount = count;
    }
  }

  return mincount;
}


/* USE_READ_COUNTS: Found to be nearly identical to coverage */

static double
gene_density (double *exon_density_1,
#ifdef USE_READ_COUNTS
	      double *exon_density_2, Bamreader_T bamreader, char *chr,
#endif
	      List_T introns, Genomicpos_T genestart, Genomicpos_T geneend,
	      int insertlength, int readlength, int min_overhang,
	      long int *tally_matches_low, long int *tally_mismatches_low,
	      long int *tally_matches_high, long int *tally_mismatches_high) {
  double intron_sum = 0.0, exon_sum = 0.0;
  int nintrons = 0;
  List_T p;
  Intron_T intron;

  int querypos, npositions;
#ifdef USE_READ_COUNTS
  int read_counts = 0, exon_npositions;
#endif
  Genomicpos_T exonstart, exonend, pos;
  double *query_counts_low, *query_counts_high;


  /* Compute intron density */
  for (p = introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    intron_sum += Splice_density(intron->splice,genestart,geneend,insertlength,readlength,
				 min_overhang);
    nintrons++;
  }


  npositions = 0;
  if (genestart < geneend) {
    /* Count positions */
    exonstart = genestart;
    for (p = introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      exonend = intron->donorpos;
      npositions += exonend - exonstart + 1;
      exonstart = intron->acceptorpos;
    }
    exonend = geneend;
    npositions += exonend - exonstart + 1;

    if (npositions <= insertlength) {
      npositions = 0;
    } else {
      query_counts_low = (double *) CALLOC(npositions,sizeof(double));
      query_counts_high = (double *) CALLOC(npositions,sizeof(double));

      /* Fill query_counts */
      querypos = 0;
      exonstart = genestart;
      for (p = introns; p != NULL; p = List_next(p)) {
	intron = (Intron_T) List_head(p);
	exonend = intron->donorpos;
#ifdef USE_READ_COUNTS
	read_counts += Bamread_nreads(&exon_npositions,bamreader,chr,exonstart,exonend);
#endif
	for (pos = exonstart; pos <= exonend; pos++) {
	  query_counts_low[querypos] = (double) (tally_matches_low[pos] + tally_mismatches_low[pos]);
	  query_counts_high[querypos] = (double) (tally_matches_high[pos] + tally_mismatches_high[pos]);
	  querypos++;
	}
	exonstart = intron->acceptorpos;
      }

      exonend = geneend;
#ifdef USE_READ_COUNTS
      read_counts += Bamread_nreads(&exon_npositions,bamreader,chr,exonstart,exonend);
#endif
      for (pos = exonstart; pos <= exonend; pos++) {
	query_counts_low[querypos] = (double) (tally_matches_low[pos] + tally_mismatches_low[pos]);
	query_counts_high[querypos] = (double) (tally_matches_high[pos] + tally_mismatches_high[pos]);
	querypos++;
      }


#ifdef WEIGHT_ENDS
      /* Weight ends */
      for (querypos = 0, i = 1; querypos < readlength; querypos++, i++) {
	query_counts_low[querypos] *= (double) readlength/(double) i;
      }
      for (querypos = npositions - insertlength + readlength - 1, i = 1;
	   querypos >= npositions - insertlength; querypos--, i++) {
	if (querypos >= 0) {
	  query_counts_low[querypos] *= (double) readlength/(double) i;
	}
      }

      for (querypos = npositions - 1, i = 1; querypos >= npositions - readlength;
	   querypos--, i++) {
	if (querypos >= 0) {
	  query_counts_high[querypos] *= (double) readlength/(double) i;
	}
      }
      for (querypos = insertlength - readlength, i = 1;
	   querypos < insertlength; querypos++, i++) {
	query_counts_high[querypos] *= (double) readlength/(double) i;
      }
#endif

    }

  } else {
    /* Count positions */
    exonstart = genestart;
    for (p = introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      exonend = intron->donorpos;
      npositions += exonstart - exonend + 1;
      exonstart = intron->acceptorpos;
    }
    exonend = geneend;
    npositions += exonstart - exonend + 1;

    if (npositions <= insertlength) {
      npositions = 0;
    } else {
      query_counts_low = (double *) CALLOC(npositions,sizeof(double));
      query_counts_high = (double *) CALLOC(npositions,sizeof(double));


      /* Fill query counts */
      querypos = 0;
      exonstart = genestart;
      for (p = introns; p != NULL; p = List_next(p)) {
	intron = (Intron_T) List_head(p);
	exonend = intron->donorpos;
#ifdef USE_READ_COUNTS
	read_counts += Bamread_nreads(&exon_npositions,bamreader,chr,exonstart,exonend);
#endif
	for (pos = exonstart; pos >= exonend; pos--) {
	  query_counts_low[querypos] = (double) (tally_matches_low[pos] + tally_mismatches_low[pos]);
	  query_counts_high[querypos] = (double) (tally_matches_high[pos] + tally_mismatches_high[pos]);
	  querypos++;
	}
	exonstart = intron->acceptorpos;
      }

      exonend = geneend;
#ifdef USE_READ_COUNTS
      read_counts += Bamread_nreads(&exon_npositions,bamreader,chr,exonstart,exonend);
#endif
      for (pos = exonstart; pos >= exonend; pos--) {
	query_counts_low[querypos] = (double) (tally_matches_low[pos] + tally_mismatches_low[pos]);
	query_counts_high[querypos] = (double) (tally_matches_high[pos] + tally_mismatches_high[pos]);
	querypos++;
      }


#ifdef WEIGHT_ENDS
      /* Weight ends */
      for (querypos = 0, i = 1; querypos < readlength; querypos++, i++) {
	query_counts_high[querypos] *= (double) readlength/(double) i;
      }
      for (querypos = npositions - insertlength + readlength - 1, i = 1;
	   querypos >= npositions - insertlength; querypos--, i++) {
	if (querypos >= 0) {
	  query_counts_high[querypos] *= (double) readlength/(double) i;
	}
      }

      for (querypos = npositions - 1, i = 1; querypos >= npositions - readlength;
	   querypos--, i++) {
	if (querypos >= 0) {
	  query_counts_low[querypos] *= (double) readlength/(double) i;
	}
      }
      for (querypos = insertlength - readlength, i = 1;
	   querypos < insertlength; querypos++, i++) {
	query_counts_low[querypos] *= (double) readlength/(double) i;
      }
#endif

    }
  }


  if (npositions == 0) {
    *exon_density_1 = 0.0;
#ifdef USE_READ_COUNTS
    *exon_density_2 = 0.0;
#endif

  } else {
    exon_sum = 0.0;
    for (querypos = 0; querypos < npositions; querypos++) {
      exon_sum += query_counts_low[querypos] + query_counts_high[querypos];
    }
    *exon_density_1 = exon_sum/(double) (readlength + readlength)/(double) (npositions - insertlength);
#ifdef USE_READ_COUNTS
    *exon_density_2 = (double) read_counts/2.0/(double) (npositions - insertlength);
#endif

    FREE(query_counts_high);
    FREE(query_counts_low);
  }

  /* intron_density */
  return (double) intron_sum/(double) nintrons;
}


static double
avg_splice_density (List_T introns, Genomicpos_T genestart, Genomicpos_T geneend, int insertlength,
		    int readlength, int min_overhang) {
  double sum = 0.0;
  int nintrons = 0;
  List_T p;
  Intron_T intron;

  for (p = introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    sum += Splice_count(intron->splice);
    nintrons++;
  }

  if (nintrons == 0) {
    return 0.0;
  } else {
    return (double) sum/(double) nintrons;
  }
}


static void
print_path (FILE *fp, List_T introns) {
  List_T p;
  Intron_T intron;

  if ((p = introns) == NULL) {
    fprintf(fp,"NA");
  } else {
    p = introns;
    intron = (Intron_T) List_head(p);
    fprintf(fp,"%u..%u",intron->donorpos,intron->acceptorpos);

    for (p = List_next(p); p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      fprintf(fp,",%u..%u",intron->donorpos,intron->acceptorpos);
    }
  }
}


/* Best way to define path is by using introns, not exons */
static void
Path_print (FILE *fp, Path_T path, double mainpath_density, int npairs, char *chr, 
	    Genomicpos_T main_intron_pathstart, Genomicpos_T main_intron_pathend,
	    Genomicpos_T genestart, Genomicpos_T geneend, IIT_T knowngenes_iit,
	    int *genei, int *known_genei, int *novel_genei, int pathnum, bool forwardp,
	    long int *tally_matches_low, long int *tally_mismatches_low,
	    long int *tally_matches_high, long int *tally_mismatches_high,
	    int insertlength, int readlength, int min_overhang, bool primaryp) {
  List_T introns, p;
  Intron_T intron;
  double intron_density = 0.0, main_density, alt_density;
#ifdef USE_READ_COUNTS
  double exon_density_2;
#endif
  long int retained_coverage;
  bool highexpr_p = false;
  Genomicpos_T last_chrpos;
  int nreads_intron, npositions, i;
  double rpkm, rpk;
  char *genelabel;
  bool allocp;

  if (path->type == MAIN) {
    introns = path->main_introns;
  } else {
    introns = path->alt_introns;
  }
    
  if (path->type == MAIN) {
    if (primaryp == false) fprintf(fp,"#");
    /* Previously identified gene by using genei.  Now using intron path start/end. */
    if (path->knownp == true) {
      *genei = (*known_genei)++;
      fprintf(fp,">Gene.known.chr%s.%s.%u..%u %s:%u..%u\n",
	      chr,forwardp == true ? "fwd" : "rev",main_intron_pathstart,main_intron_pathend,
	      chr,genestart,geneend);
    } else {
      *genei = (*novel_genei)++;
      fprintf(fp,">Gene.novel.chr%s.%s.%u..%u %s:%u..%u\n",
	      chr,forwardp == true ? "fwd" : "rev",main_intron_pathstart,main_intron_pathend,
	      chr,genestart,geneend);
    }

    if (primaryp == false) fprintf(fp,"#");
    rpk = 1000.0*mainpath_density;
    rpkm = 1000000.0*rpk/(double) npairs;

    fprintf(fp,"type=%s rpkm=%f\n",Introntype_string(path->type),rpkm);

  } else if (path->type == RETAINED_INTRON || path->type == TERMINAL_RETAINED) {
    /* Compare against mainpath density */
    intron = (Intron_T) List_head(path->main_introns);

    retained_coverage = tally_mincount(intron->donorpos,intron->acceptorpos,
				       tally_matches_low,tally_mismatches_low,
				       tally_matches_high,tally_mismatches_high);
    alt_density = (double) retained_coverage/(double) (readlength + readlength);


    if (path->knownp == false && alt_density <= 0.05 * mainpath_density) {
      /* Shouldn't get here, since already checked */
      return;
    } else {
      if (primaryp == false) fprintf(fp,"#");
      /* Previously identified alternate by using genei.pathnum.  Now using intron path start/end and alt path */
      fprintf(fp,">Alt.%s.chr%s.%s.%u..%u.%s.",
	      path->knownp == true ? "known" : "novel",chr,
	      forwardp == true ? "fwd" : "rev",main_intron_pathstart,main_intron_pathend,
	      Introntype_string(path->type));
      print_path(fp,path->main_introns);
      fprintf(fp,"=>");
      fprintf(fp,"NA"); /* alt_introns is NULL */
      fprintf(fp," %s:%u..%u %s\n",chr,genestart,geneend,Introntype_string(path->type));

      if (primaryp == false) fprintf(fp,"#");
      if (alt_density + mainpath_density == 0) {
	fprintf(fp,"type=%s frac=NA ref=%f alt=%f\n",
		Introntype_string(path->type),mainpath_density,alt_density);
      } else {
	fprintf(fp,"type=%s frac=%f ref=%f alt=%f\n",
		Introntype_string(path->type),alt_density/(alt_density + mainpath_density),mainpath_density,alt_density);
      }
    }

  } else {
    main_density = avg_splice_density(path->main_introns,genestart,geneend,insertlength,readlength,min_overhang);
    alt_density = avg_splice_density(path->alt_introns,genestart,geneend,insertlength,readlength,min_overhang);
    if (path->knownp == false && (alt_density + main_density == 0.0 || alt_density <= 0.05 * main_density)) {
      return;
    } else {
      if (primaryp == false) fprintf(fp,"#");
      /* Previously identified alternate by using genei.pathnum.  Now using intron path start/end and alt path */
      fprintf(fp,">Alt.%s.chr%s.%s.%u..%u.%s.",
	      path->knownp == true ? "known" : "novel",chr,
	      forwardp == true ? "fwd" : "rev",main_intron_pathstart,main_intron_pathend,
	      Introntype_string(path->type));
      print_path(fp,path->main_introns);
      fprintf(fp,"=>");
      print_path(fp,introns);
      fprintf(fp," %s:%u..%u %s\n",chr,genestart,geneend,Introntype_string(path->type));

      if (primaryp == false) fprintf(fp,"#");
      if (alt_density + main_density < MIN_COUNT) {
	fprintf(fp,"type=%s frac=NA ref=%.1f alt=%.1f\n",Introntype_string(path->type),main_density,alt_density);
      } else {
	fprintf(fp,"type=%s frac=%f ref=%.1f alt=%.1f\n",
		Introntype_string(path->type),
		(alt_density + PSEUDOCOUNT)/(alt_density + main_density + PSEUDOCOUNT + PSEUDOCOUNT),
		main_density,alt_density);
      }
    }
  }

  /* Print exon and introns */
  last_chrpos = path->genestart; /* was path->pathstart */
  for (p = introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    if (primaryp == false) fprintf(fp,"#");
#if 0
    nreads_intron = Bamread_nreads(&npositions,bamreader,chr,intron->donorpos,intron->acceptorpos);
    fprintf(fp,"%u %u %d/%d %d %d %d /%d\n",
	    last_chrpos,intron->donorpos,
	    nreads_intron,npositions,
	    Splice_known_count(intron->splice),
	    Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),
	    Splice_npositions(intron->splice));
#else
    fprintf(fp,"%u %u %d %d %d /%d",
	    last_chrpos,intron->donorpos,
	    Splice_known_count(intron->splice),
	    Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),
	    Splice_npositions(intron->splice,genestart,geneend,/*insertlength*/200,/*readlength*/75,
			      /*min_overhang*/15));
    if (intron->ngenei > 0) {
      genelabel = IIT_label(knowngenes_iit,intron->genei[0],&allocp);
      fprintf(fp," %s",genelabel);
      if (allocp) FREE(genelabel);
      for (i = 1; i < intron->ngenei; i++) {
	genelabel = IIT_label(knowngenes_iit,intron->genei[i],&allocp);
	fprintf(fp,",%s",genelabel);
	if (allocp) FREE(genelabel);
      }
    }
    fprintf(fp,"\n");
#endif

    last_chrpos = intron->acceptorpos;
  }

  if (primaryp == false) fprintf(fp,"#");
  fprintf(fp,"%u %u\n",last_chrpos,geneend /* was path->pathend */);

  return;
}


#if 0
static bool
check_path (Path_T path) {
  List_T p;
  Intron_T intron;

  /* Find altstart and altend */
  for (p = path->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    /*acceptor = intron->site; */
    /* donor = intron->prev_site; */
    if (intron->introntype == MAIN) {
      /* Skip */
    } else if (intron->prev_site->switchp == true && path->type != SWITCH_EXON) {
      debug5(printf("check_path finds that %u was involved in a switch\n",intron->prev_site->chrpos));
      return false;
    } else if (intron->site->switchp == true && path->type != SWITCH_EXON) {
      debug5(printf("check_path finds that %u was involved in a switch\n",intron->site->chrpos));
      return false;
    }
  }

  return true;
}
#endif


/************************************************************************
 *   Terminals
 ************************************************************************/


static Path_T
get_path_main (Acceptor_T acceptor, bool forwardp, bool knownp) {
  Path_T mainpath;
  Intron_T main_intron;
  Donor_T donor;
  Donor_alt_T alt;
  List_T p;

  mainpath = Path_new_main(/*introns*/NULL,/*pathtype*/MAIN,forwardp,/*altp*/false,knownp);
  while (acceptor != NULL) {
    main_intron = acceptor->main_intron;
    if (main_intron->validp == true) {
      mainpath->main_introns = List_push(mainpath->main_introns,(void *) main_intron);
    }
    if ((donor = acceptor->prev) == NULL) {
      fprintf(stderr,"Singleton acceptor at %u\n",acceptor->chrpos);
      abort();
    } else {
      acceptor = donor->prev;
      if (acceptor != NULL) {
	acceptor->main_exon = donor->main_exon; /* Make forward link for main */
	for (p = donor->alt_structs; p != NULL; p = List_next(p)) {
	  alt = (Donor_alt_T) List_head(p);
	  if (alt->prev_site->main_exon == NULL) {
	    alt->prev_site->main_exon = donor->main_exon; /* Make forward link for alt structs */
	  }
	}
	for (p = donor->alt_sites; p != NULL; p = List_next(p)) {
	  alt = (Donor_alt_T) List_head(p);
	  if (alt->prev_site->main_exon == NULL) {
	    alt->prev_site->main_exon = donor->main_exon; /* Make forward link for alt sites */
	  }
	}
      }
    }
  }

  return mainpath;
}



#if 0
static List_T
get_paths_alt_old (Path_T mainpath) {
  List_T result = NULL, p, q;
  Donor_T donor;
  Donor_alt_T alt1;
  Acceptor_T acceptor;
  Acceptor_alt_T alt;
  Intron_T main_intron = NULL;
  bool forwardp;

  forwardp = mainpath->forwardp;
  for (q = mainpath->introns; q != NULL; q = List_next(q)) {
    main_intron = (Intron_T) List_head(q);

    donor = main_intron->prev_site;
    for (p = donor->alt_sites; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      result = List_push(result,(void *) Path_donor_alt(alt1,forwardp));
    }

    for (p = donor->alt_switches; p != NULL; p = List_next(p)) {
      alt1 = (Donor_alt_T) List_head(p);
      result = List_push(result,(void *) Path_donor_alt(alt1,forwardp));
    }


    acceptor = main_intron->site;
    for (p = acceptor->alt_structs; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      result = List_push(result,(void *) Path_acceptor_alt(alt,forwardp));
    }

    for (p = acceptor->alt_sites; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      result = List_push(result,(void *) Path_acceptor_alt(alt,forwardp));
    }

    for (p = acceptor->alt_switches; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      result = List_push(result,(void *) Path_acceptor_alt(alt,forwardp));
    }
  }

  return result;
}
#endif


static List_T
get_paths_alt (List_T paths, Acceptor_T acceptor, bool forwardp, bool knownp) {
  List_T result = NULL, p, q;
  Donor_T donor;
  Donor_alt_T alt1;
  Acceptor_alt_T alt;
  Exon_T exon;
  Intron_T intron;

  if (acceptor->paths_usedp == false) {
    acceptor->paths_usedp = true;

    for (p = acceptor->alt_structs; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      if (knownp == true || alt->knownp == false) {
	paths = List_push(paths,(void *) Path_acceptor_alt(alt,forwardp));
      }
    }

    for (p = acceptor->alt_sites; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      if (knownp == true || alt->knownp == false) {
	paths = List_push(paths,(void *) Path_acceptor_alt(alt,forwardp));
      }
    }

    for (p = acceptor->alt_switches; p != NULL; p = List_next(p)) {
      alt = (Acceptor_alt_T) List_head(p);
      if (knownp == true || alt->knownp == false) {
	paths = List_push(paths,(void *) Path_acceptor_alt(alt,forwardp));
      }
    }

    for (q = acceptor->introns; q != NULL; q = List_next(q)) {
      intron = (Intron_T) List_head(q);
      donor = intron->prev_site;
      for (p = donor->alt_sites; p != NULL; p = List_next(p)) {
	alt1 = (Donor_alt_T) List_head(p);
	if (knownp == true || alt1->knownp == false) {
	  paths = List_push(paths,(void *) Path_donor_alt(alt1,forwardp));
	}
      }

      for (p = donor->alt_switches; p != NULL; p = List_next(p)) {
	alt1 = (Donor_alt_T) List_head(p);
	if (knownp == true || alt1->knownp == false) {
	  paths = List_push(paths,(void *) Path_donor_alt(alt1,forwardp));
	}
      }

      for (p = donor->exons; p != NULL; p = List_next(p)) {
	exon = (Exon_T) List_head(p);
	paths = get_paths_alt(paths,/*acceptor*/exon->prev_site,forwardp,knownp);
      }
    }
  }

  return paths;
}





static void
Path_add_genebounds (Path_T path,long int *tally_matches_low, long int *tally_mismatches_low,	
		     long int *tally_matches_high, long int *tally_mismatches_high,
		     IIT_T end_exons_iit, char *chr, Genomicpos_T chrlength, bool forwardp) {
  Genomicpos_T pathstart, pathend;
  List_T introns;
  Intron_T first_intron, last_intron;
  Exon_T exon;
  
  if (path->altp == false) {
    introns = path->main_introns;
  } else {
    introns = path->alt_introns;
  }

  first_intron = (Intron_T) List_head(introns);
  pathstart = first_intron->prev_site->chrpos;

  last_intron = (Intron_T) List_last_value(introns);
  pathend = last_intron->site->chrpos;


  if (forwardp == true) {
    if (path->altp == true && (exon = first_intron->prev_site->main_exon) != NULL &&
	exon->prev_site != NULL) {
      path->genestart = exon->prev_site->chrpos;
    } else {
      path->genestart = Cappaths_solve_genestart(pathstart,tally_matches_low,tally_mismatches_low,
						 end_exons_iit,chr,chrlength,forwardp);
    }

    if (path->altp == true && (exon = last_intron->site->main_exon) != NULL &&
	exon->site != NULL) {
      path->geneend = exon->site->chrpos;
    } else {
      path->geneend = Cappaths_solve_geneend(pathend,tally_matches_high,tally_mismatches_high,
					     end_exons_iit,chr,chrlength,forwardp);
    }

  } else {
    if (path->altp == true && (exon = first_intron->prev_site->main_exon) != NULL &&
	exon->prev_site != NULL) {
      path->genestart = exon->prev_site->chrpos;
    } else {
      path->genestart = Cappaths_solve_genestart(pathstart,tally_matches_high,tally_mismatches_high,
						 end_exons_iit,chr,chrlength,forwardp);
    }

    if (path->altp == true && (exon = last_intron->site->main_exon) != NULL &&
	exon->site != NULL) {
      path->geneend = exon->site->chrpos;
    } else {
      path->geneend = Cappaths_solve_geneend(pathend,tally_matches_low,tally_mismatches_low,
					     end_exons_iit,chr,chrlength,forwardp);
    }
  }

  return;
}


static List_T
Site_get_paths (Path_T *mainpath, Acceptor_T acceptor, char *chr, int mincount_alt,
		long int *tally_matches_low, long int *tally_mismatches_low,
		long int *tally_matches_high, long int *tally_mismatches_high,
		long int *primary_extents, long int *crosshyb_extents,
		IIT_T end_exons_iit, Genomicpos_T chrlength, bool forwardp,
		bool altpaths_p, bool knownp) {
  List_T altpaths, orig_altpaths, p;
  Path_T path;

  if (acceptor->known_mainpath != NULL) {
    *mainpath = Path_copy(acceptor->known_mainpath);

    orig_altpaths = NULL;
    for (p = acceptor->known_altpaths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      orig_altpaths = List_push(orig_altpaths,(void *) Path_copy(path));
    }

    altpaths = get_paths_alt(orig_altpaths,acceptor,forwardp,knownp);
    for (p = altpaths; p != orig_altpaths; p = List_next(p)) {
      path = (Path_T) List_head(p);
      Path_add_genebounds(path,tally_matches_low,tally_mismatches_low,
			  tally_matches_high,tally_mismatches_high,
			  end_exons_iit,chr,chrlength,forwardp);
    }
    return altpaths;

  } else if (acceptor->known_termpath != NULL) {
    *mainpath = Path_copy(acceptor->known_termpath);
    return (List_T) NULL;

  } else {
    *mainpath = get_path_main(acceptor,forwardp,knownp);
    (*mainpath)->main_introns = Intron_path_trim((*mainpath)->main_introns);
    if ((*mainpath)->main_introns == NULL) {
      Path_free(&(*mainpath));
      *mainpath = (Path_T) NULL;
      return (List_T) NULL;

    } else {
      Path_add_genebounds(*mainpath,tally_matches_low,tally_mismatches_low,
			  tally_matches_high,tally_mismatches_high,
			  end_exons_iit,chr,chrlength,forwardp);
      if (Path_primaryp(*mainpath,chr,primary_extents,crosshyb_extents) == false) {
	(*mainpath)->primaryp = false;
	return (List_T) NULL;
	
      } else if (altpaths_p == false) {
	(*mainpath)->primaryp = true;
	return (List_T) NULL;

      } else {
	(*mainpath)->primaryp = true;
	altpaths = get_paths_alt(/*altpaths*/NULL,acceptor,forwardp,knownp);
	for (p = altpaths; p != NULL; p = List_next(p)) {
	  path = (Path_T) List_head(p);
	  Path_add_genebounds(path,tally_matches_low,tally_mismatches_low,
			      tally_matches_high,tally_mismatches_high,
			      end_exons_iit,chr,chrlength,forwardp);
	}
	
	return altpaths;
      }
    }
  }
}


#if 0
static Path_T
Path_new_termpath_old (List_T term_introns, int n_alt_introns,
		       List_T main_introns, int n_main_introns,
		       Genomicpos_T geneend, Introntype_T terminal_type, bool knownp) {
  Path_T new;
  List_T p;
  Donor_T donor;
  Intron_T intron;
  bool usedp;
  Genomicpos_T last_chrpos, altstart, altend;
  int n_main_introns, n_alt_introns;

  /* Find altstart and altend */
  usedp = true;
  p = introns;
  while (p != NULL && usedp == true) {
    intron = (Intron_T) List_head(p);
    if (intron->site->usedp == false || intron->prev_site->usedp == false) {
      usedp = false;
    } else {
      p = List_next(p);
    }
  }

  donor = intron->prev_site;
  if (donor->prev == NULL) {
    return (Path_T) NULL;
  } else {
    new = (Path_T) MALLOC(sizeof(*new));
    new->knownp = knownp;
    new->type = terminal_type;
    new->main_introns = (List_T) NULL;
    new->alt_introns = (List_T) NULL;
    new->primaryp = true;
    new->altp = true;

    altstart = donor->prev->chrpos;
    altend = geneend;

    new->genestart = altstart;
    new->geneend = geneend;

    p = introns;
    usedp = true;
    while (p != NULL && usedp == true) {
      intron = (Intron_T) List_head(p);
      if (intron->site->usedp == false || intron->prev_site->usedp == false) {
	usedp = false;
      } else {
	p = List_next(p);
      }
    }

    last_chrpos = altstart;
    while (p != NULL) {
      intron = (Intron_T) List_head(p);
      new->alt_introns = List_push(new->alt_introns,(void *) intron);
      p = List_next(p);
    }

    new->alt_introns = List_reverse(new->alt_introns);
    return new;
  }
}
#endif


static Path_T
Path_new_termpath (Acceptor_T terminal, int n_alt_introns,
		   Acceptor_T main_terminal, int n_main_introns,
		   Genomicpos_T geneend, Introntype_T terminal_type, bool knownp,
		   bool retainedp) {
  Path_T new;
  Intron_T intron;
  Genomicpos_T donorpos;
  int i;

  donorpos = terminal->prev->chrpos;

  new = (Path_T) MALLOC(sizeof(*new));
  new->knownp = knownp;
  new->type = terminal_type;

  new->primaryp = true;
  new->altp = true;

  new->alt_introns = (List_T) NULL;
  for (i = 0; i < n_alt_introns; i++) {
    intron = terminal->main_intron;
    new->alt_introns = List_push(new->alt_introns,(void *) intron);
    terminal = terminal->prev->prev;
  }
  new->alt_introns = List_reverse(new->alt_introns);

  if (terminal == NULL) {
    List_free(&new->alt_introns);
    FREE(new);
    return (Path_T) NULL;
  } else if (terminal->prev->prev == NULL) {
    List_free(&new->alt_introns);
    FREE(new);
    return (Path_T) NULL;
  } else {
    new->genestart = terminal->prev->prev->chrpos;
    new->geneend = geneend;

    new->main_introns = (List_T) NULL;
    if (retainedp == true) {
      for (i = 0; i < n_main_introns; i++) {
	intron = main_terminal->main_intron;
	new->main_introns = List_push(new->main_introns,(void *) intron);
	main_terminal = main_terminal->prev->prev;
      }
    } else {
      /* Don't consider full mainpath, just the intron that initially branches away */
      i = 0;
      while (i < n_main_introns && ((Intron_T) main_terminal->main_intron)->donorpos != donorpos) {
	intron = main_terminal->main_intron;
	main_terminal = main_terminal->prev->prev;
	i++;
      }
      if (i < n_main_introns) {
	intron = main_terminal->main_intron;
	new->main_introns = List_push(new->main_introns,(void *) intron);
      }
    }
    new->main_introns = List_reverse(new->main_introns);
    
    return new;
  }
}




static Genomicpos_T
donor_distalbound (Donor_T this);
static Genomicpos_T
acceptor_distalbound (Acceptor_T this);

static Genomicpos_T
donor_distalbound (Donor_T this) {
  if (this->prev == NULL) {
    return this->chrpos;
  } else {
    return acceptor_distalbound(this->prev);
  }
}

static Genomicpos_T
acceptor_distalbound (Acceptor_T this) {
  if (this->prev == NULL) {
    return this->chrpos;
  } else {
    return donor_distalbound(this->prev);
  }
}


/************************************************************************/

struct Gene_T {
  bool forwardp;
  Path_T mainpath;
  List_T altpaths;
  List_T alt_termpaths;
  List_T ret_termpaths;
};


void
Gene_free (Gene_T *old) {
  List_T p;
  Path_T path;

  Path_free(&(*old)->mainpath);

  for (p = (*old)->altpaths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    Path_free(&path);
  }
  List_free(&(*old)->altpaths);

  for (p = (*old)->alt_termpaths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    Path_free(&path);
  }
  List_free(&(*old)->alt_termpaths);

  for (p = (*old)->ret_termpaths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    Path_free(&path);
  }
  List_free(&(*old)->ret_termpaths);

  FREE(*old);
  return;
}


void
Splicegraph_genes_gc (List_T *genes) {
  List_T p;
  Gene_T gene;

  for (p = *genes; p != NULL; p = List_next(p)) {
    gene = (Gene_T) List_head(p);
    Gene_free(&gene);
  }
  List_free(&(*genes));
  return;
}



Gene_T
Gene_new (Path_T mainpath, List_T altpaths, List_T alt_termpaths, List_T ret_termpaths,
	  bool forwardp) {
  Gene_T new = (Gene_T) MALLOC(sizeof(*new));

  new->forwardp = forwardp;
  new->mainpath = mainpath;
  new->altpaths = altpaths;
  new->alt_termpaths = alt_termpaths;
  new->ret_termpaths = ret_termpaths;
  return new;
}


/************************************************************************
 *   Splicegraph
 ************************************************************************/

#define T Splicegraph_T
struct T {
  Uinttable_T donor_fwd_sitestable;
  Uinttable_T acceptor_fwd_sitestable;
  Uinttable_T donor_rev_sitestable;
  Uinttable_T acceptor_rev_sitestable;

  Donor_T *sites_donor_fwd;
  Acceptor_T *sites_acceptor_fwd;
  Donor_T *sites_donor_rev;
  Acceptor_T *sites_acceptor_rev;

  int nsites_donor_fwd;
  int nsites_acceptor_fwd;
  int nsites_donor_rev;
  int nsites_acceptor_rev;



};


void
Splicegraph_free (T *old) {
  Genomicpos_T *keys;
  Donor_T donor;
  Acceptor_T acceptor;
  int n, i;

  if ((n = Uinttable_length((*old)->donor_fwd_sitestable)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys((*old)->donor_fwd_sitestable,/*sortp*/false);
    for (i = 0; i < n; i++) {
      donor = (Donor_T) Uinttable_get((*old)->donor_fwd_sitestable,keys[i]);
      Donor_free(&donor);
    }
    FREE(keys);
  }
  Uinttable_free(&(*old)->donor_fwd_sitestable);

  if ((n = Uinttable_length((*old)->acceptor_fwd_sitestable)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys((*old)->acceptor_fwd_sitestable,/*sortp*/false);
    for (i = 0; i < n; i++) {
      acceptor = (Acceptor_T) Uinttable_get((*old)->acceptor_fwd_sitestable,keys[i]);
      Acceptor_free(&acceptor);
    }
    FREE(keys);
  }
  Uinttable_free(&(*old)->acceptor_fwd_sitestable);

  if ((n = Uinttable_length((*old)->donor_rev_sitestable)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys((*old)->donor_rev_sitestable,/*sortp*/false);
    for (i = 0; i < n; i++) {
      donor = (Donor_T) Uinttable_get((*old)->donor_rev_sitestable,keys[i]);
      Donor_free(&donor);
    }
    FREE(keys);
  }
  Uinttable_free(&(*old)->donor_rev_sitestable);

  if ((n = Uinttable_length((*old)->acceptor_rev_sitestable)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys((*old)->acceptor_rev_sitestable,/*sortp*/false);
    for (i = 0; i < n; i++) {
      acceptor = (Acceptor_T) Uinttable_get((*old)->acceptor_rev_sitestable,keys[i]);
      Acceptor_free(&acceptor);
    }
    FREE(keys);
  }
  Uinttable_free(&(*old)->acceptor_rev_sitestable);

  FREE(*old);
  return;
}



T
Splicegraph_new () {
  T new = (T) MALLOC(sizeof(*new));
  int nsplices = 10000;

  new->donor_fwd_sitestable = Uinttable_new(nsplices);
  new->acceptor_fwd_sitestable = Uinttable_new(nsplices);
  new->donor_rev_sitestable = Uinttable_new(nsplices);
  new->acceptor_rev_sitestable = Uinttable_new(nsplices);

  return new;
}



void
Splicegraph_print_genes (FILE *fp, List_T genes, int npairs, IIT_T knowngenes_iit, char *chr,
			 long int *tally_matches_low, long int *tally_mismatches_low,
			 long int *tally_matches_high, long int *tally_mismatches_high,
			 int insertlength, int readlength, int min_overhang) {

  List_T g, p;
  Gene_T gene;
  Path_T mainpath, path, *paths;
  Intron_T first_intron, last_intron;
  Genomicpos_T main_intron_pathstart, main_intron_pathend; /* Used to define the main gene */
  int known_genei = 0, novel_genei = 0, genei, pathnum, npaths, i;
  bool primaryp;

  for (g = genes; g != NULL; g = List_next(g)) {
    gene = (Gene_T) List_head(g);
    pathnum = 1;

    mainpath = gene->mainpath;
    if (mainpath->knownp == true || mainpath->density > 0.0) {
      primaryp = mainpath->primaryp;
      first_intron = (Intron_T) List_head(mainpath->main_introns);
      main_intron_pathstart = first_intron->donorpos;

      last_intron = (Intron_T) List_last_value(mainpath->main_introns);
      main_intron_pathend = last_intron->acceptorpos;

      Path_print(fp,mainpath,mainpath->density,npairs,
		 chr,main_intron_pathstart,main_intron_pathend,
		 mainpath->genestart,mainpath->geneend,knowngenes_iit,
		 &genei,&known_genei,&novel_genei,pathnum,gene->forwardp,
		 tally_matches_low,tally_mismatches_low,
		 tally_matches_high,tally_mismatches_high,
		 insertlength,readlength,min_overhang,primaryp);

      if ((npaths = List_length(gene->altpaths) + List_length(gene->alt_termpaths) + List_length(gene->ret_termpaths)) > 0) {
	paths = (Path_T *) CALLOC(npaths,sizeof(Path_T));

	npaths = 0;
	for (p = gene->altpaths; p != NULL; p = List_next(p)) {
	  paths[npaths++] = (Path_T) List_head(p);
	}
	for (p = gene->alt_termpaths; p != NULL; p = List_next(p)) {
	  paths[npaths++] = (Path_T) List_head(p);
	}
	for (p = gene->ret_termpaths; p != NULL; p = List_next(p)) {
	  paths[npaths++] = (Path_T) List_head(p);
	}
	if (gene->forwardp == true) {
	  qsort(paths,npaths,sizeof(Path_T),Path_fwd_cmp);
	} else {
	  qsort(paths,npaths,sizeof(Path_T),Path_rev_cmp);
	}

	for (i = 0; i < npaths; i++) {
	  path = paths[i];
	  Path_print(fp,path,mainpath->density,npairs,
		     chr,main_intron_pathstart,main_intron_pathend,
		     path->genestart,path->geneend,knowngenes_iit,
		     &genei,&known_genei,&novel_genei,pathnum++,gene->forwardp,
		     tally_matches_low,tally_mismatches_low,
		     tally_matches_high,tally_mismatches_high,
		     insertlength,readlength,min_overhang,primaryp);
	}
	FREE(paths);
      }
    }
  }

  return;
}


List_T
Splicegraph_splices (List_T genes) {
  List_T splices = NULL, g, p, q;
  Gene_T gene;
  Path_T path;
  Intron_T intron;
  int sign;

  for (g = genes; g != NULL; g = List_next(g)) {
    gene = (Gene_T) List_head(g);

    path = gene->mainpath;
    if (path->forwardp == true) {
      sign = +1;
    } else {
      sign = -1;
    }

    for (q = path->main_introns; q != NULL; q = List_next(q)) {
      intron = (Intron_T) List_head(q);
      splices = List_push(splices,(void *) Splice_new_known(intron->donorpos,intron->acceptorpos,sign,/*mainpath_p*/true));
    }

    for (p = gene->altpaths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      for (q = path->alt_introns; q != NULL; q = List_next(q)) {
	intron = (Intron_T) List_head(q);
	splices = List_push(splices,(void *) Splice_new_known(intron->donorpos,intron->acceptorpos,sign,/*mainpath_p*/false));
      }
    }

    for (p = gene->alt_termpaths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      for (q = path->alt_introns; q != NULL; q = List_next(q)) {
	intron = (Intron_T) List_head(q);
	splices = List_push(splices,(void *) Splice_new_known(intron->donorpos,intron->acceptorpos,sign,/*mainpath_p*/false));
      }
    }

    /* No need to process gene->ret_termpaths, which should have no introns */

  }

  return splices;
}


/************************************************************************/


static Genomicpos_T
donor_used_p (Donor_T this);
static Genomicpos_T
acceptor_used_p (Acceptor_T this);

static Genomicpos_T
donor_used_p (Donor_T this) {
  if (this == NULL) {
    return 0U;
  } else if (this->usedp > 0U) {
    return this->usedp;
  } else {
    return acceptor_used_p(this->prev);
  }
}

static Genomicpos_T
acceptor_used_p (Acceptor_T this) {
  if (this == NULL) {
    return 0U;
  } else if (this->usedp > 0U) {
    return this->usedp;
  } else {
    return donor_used_p(this->prev);
  }
}

static Intron_T
acceptor_used_point (Acceptor_T this) {
  if (this == NULL) {
    return (Intron_T) NULL;
  } else if (this->prev == NULL) {
    return (Intron_T) NULL;
  } else if (this->usedp > 0 && this->prev->usedp > 0) {
    return this->main_intron;
  } else {
    return acceptor_used_point(this->prev->prev);
  }
}


static void
donor_mark_used (Donor_T this, Genomicpos_T terminal_chrpos);
static void
acceptor_mark_used (Acceptor_T this, Genomicpos_T terminal_chrpos);


static void
donor_mark_used (Donor_T this, Genomicpos_T terminal_chrpos) {
  List_T p;
  Exon_T exon;

  if (this == NULL) {
    return;

  } else if (this->usedp != 0) {
    /* Already marked */
    return;

  } else {
    this->usedp = terminal_chrpos;
#if 0
    acceptor_mark_used(this->prev,terminal_chrpos);
#else
    for (p = this->exons; p != NULL; p = List_next(p)) {
      exon = (Exon_T) List_head(p);
      acceptor_mark_used(exon->prev_site,terminal_chrpos);
    }
#endif
    return;
  }
}


static void
acceptor_mark_used (Acceptor_T this, Genomicpos_T terminal_chrpos) {
  List_T p;
  Intron_T intron;

  if (this == NULL) {
    return;

  } else if (this->usedp != 0) {
    /* Already marked */
    return;

  } else {
    this->usedp = terminal_chrpos;
#if 0
    donor_mark_used(this->prev,terminal_chrpos);
#else
    for (p = this->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      donor_mark_used(intron->prev_site,terminal_chrpos);
    }
#endif
    return;
  }
}

static int
acceptor_nintrons (Acceptor_T this, Genomicpos_T donor_chrpos, Genomicpos_T acceptor_chrpos) {
  if (this == NULL) {
    return 0;
  } else if (this->prev == NULL) {
    return 0;
  } else if (this->chrpos == acceptor_chrpos && this->prev->chrpos == donor_chrpos) {
    return 0;
  } else {
    return 1 + acceptor_nintrons(this->prev->prev,donor_chrpos,acceptor_chrpos);
  }
}


static Genomicpos_T
flat_chrpos_low (double *window_diff, Genomicpos_T low, Genomicpos_T high,
		 int allowed_zeroes) {
  Genomicpos_T chrpos;
  bool zerop = false;
  int nconsecutive_zeroes = 0;

  debug6(printf("starting flat_chrpos_low from high %u downto low %u\n",high,low));
  for (chrpos = high; chrpos >= low; chrpos--) {
    debug6(printf("  chrpos %u, window_diff %f\n",chrpos,window_diff[chrpos]));
    if (zerop == true) {
      if (window_diff[chrpos] == 0.0) {
	/* Continuing zero region */
	nconsecutive_zeroes++;
	if (nconsecutive_zeroes >= allowed_zeroes) {
	  debug6(printf("\n"));
	  return chrpos;
	}
      } else {
	/* Ending zero region */
	zerop = false;
      }

    } else {
      if (window_diff[chrpos] == 0.0) {
	/* Starting zero region */
	nconsecutive_zeroes = 1;
	zerop = true;
      } else {
	/* Continuing non-zero region */
      }
    }
  }

  debug6(printf("\n"));
  return low;
}


static Genomicpos_T
flat_chrpos_high (double *window_diff, Genomicpos_T low, Genomicpos_T high,
		  int allowed_zeroes) {
  Genomicpos_T chrpos;
  bool zerop = false;
  int nconsecutive_zeroes = 0;

  debug6(printf("starting flat_chrpos_high from low %u to high %u\n",low,high));
  for (chrpos = low; chrpos <= high; chrpos++) {
    debug6(printf("  chrpos %u, window_diff %f\n",chrpos,window_diff[chrpos]));
    if (zerop == true) {
      if (window_diff[chrpos] == 0.0) {
	/* Continuing zero region */
	nconsecutive_zeroes++;
	if (nconsecutive_zeroes >= allowed_zeroes) {
	  debug6(printf("\n"));
	  return chrpos;
	}
      } else {
	/* Ending zero region */
	zerop = false;
      }

    } else {
      if (window_diff[chrpos] == 0.0) {
	/* Starting zero region */
	nconsecutive_zeroes = 1;
	zerop = true;
      } else {
	/* Continuing non-zero region */
      }
    }
  }

  debug6(printf("\n"));
  return high;
}


static double
terminal_density (Acceptor_T terminal, int n, Genomicpos_T genestart, Genomicpos_T geneend,
		  int insertlength, int readlength, int min_overhang) {
  double sum = 0.0;
  Intron_T intron;
  int i;

  for (i = 0; i < n; i++) {
    intron = terminal->main_intron;
    sum += Splice_count(intron->splice);
    terminal = terminal->prev->prev;
  }

  return sum/(double) n;
}



static List_T
get_paths_term_alt_known (List_T alt_termpaths, List_T terminals, Acceptor_T main_terminal, char *chr,
			  long int *tally_matches_low, long int *tally_mismatches_low,
			  long int *tally_matches_high, long int *tally_mismatches_high,
			  long int *primary_extents, long int *crosshyb_extents,
			  IIT_T end_exons_iit, Genomicpos_T chrlength, int mincount_alt,
			  bool forwardp) {
  List_T p;
  Path_T termpath, path;
  Acceptor_T terminal;
  Intron_T join_point;
  Genomicpos_T join_donor_chrpos, join_acceptor_chrpos;
  Genomicpos_T geneend1, geneend2;
  int n_main_introns, n_alt_introns;


  for (p = terminals; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    join_point = acceptor_used_point(terminal);
    if (join_point == NULL) {
      join_donor_chrpos = 0U;
      join_acceptor_chrpos = 0U;
    } else {
      join_donor_chrpos = join_point->prev_site->chrpos;
      join_acceptor_chrpos = join_point->site->chrpos;
    }

    n_main_introns = acceptor_nintrons(main_terminal,join_donor_chrpos,join_acceptor_chrpos);
    n_alt_introns = acceptor_nintrons(terminal,join_donor_chrpos,join_acceptor_chrpos);
    
    if (n_main_introns == 1 && n_alt_introns == 1) {
      if (forwardp == true) {
	geneend1 = Cappaths_solve_geneend(main_terminal->chrpos,tally_matches_high,tally_mismatches_high,
					  end_exons_iit,chr,chrlength,forwardp);
	geneend2 = Cappaths_solve_geneend(terminal->chrpos,tally_matches_high,tally_mismatches_high,
					  end_exons_iit,chr,chrlength,forwardp);
      } else {
	geneend1 = Cappaths_solve_geneend(main_terminal->chrpos,tally_matches_low,tally_mismatches_low,
					  end_exons_iit,chr,chrlength,forwardp);
	geneend2 = Cappaths_solve_geneend(terminal->chrpos,tally_matches_low,tally_mismatches_low,
					  end_exons_iit,chr,chrlength,forwardp);
      }
      if (geneend1 == geneend2) {
	/* Alternate splice site */
	Site_get_paths(&path,terminal,chr,mincount_alt,
		       tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		       primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		       /*altpaths_p*/false,/*knownp*/true);
	if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					  main_terminal,n_main_introns,
					  geneend2,/*terminal_type*/TERMINAL_ALT_SITE,
					  /*knownp*/true,/*retainedp*/false)) != NULL) {
	  alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	  terminal->known_termpath = termpath;
	}
	Path_free(&path);

      } else {
	/* Alternate last exon */
	Site_get_paths(&path,terminal,chr,mincount_alt,
		       tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		       primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		       /*altpaths_p*/false,/*knownp*/true);
	if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					  main_terminal,n_main_introns,
					  geneend2,/*terminal_type*/TERMINAL_ALT_EXON,
					  /*knownp*/true,/*retainedp*/false)) != NULL) {
	  alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	  terminal->known_termpath = termpath;
	}
	Path_free(&path);
      }

    } else if ((forwardp == true && terminal->chrpos < main_terminal->chrpos) ||
	       (forwardp == false && terminal->chrpos > main_terminal->chrpos)) {
      /* Short form */
      Site_get_paths(&path,terminal,chr,mincount_alt,
		     tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		     primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		     /*altpaths_p*/false,/*knownp*/true);
      if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					main_terminal,n_main_introns,
					path->geneend,/*terminal_type*/TERMINAL_SHORT,
					/*knownp*/true,/*retainedp*/false)) != NULL) {
	alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	terminal->known_termpath = termpath;
      }
      Path_free(&path);

    } else {
      /* Long form */
      Site_get_paths(&path,terminal,chr,mincount_alt,
		     tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		     primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		     /*altpaths_p*/false,/*knownp*/true);
      if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					main_terminal,n_main_introns,
					path->geneend,/*terminal_type*/TERMINAL_LONG,
					/*knownp*/true,/*retainedp*/false)) != NULL) {
	alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	terminal->known_termpath = termpath;
      }
      Path_free(&path);
    }
  }

  return alt_termpaths;
}


static List_T
get_paths_term_alt_obs (List_T alt_termpaths, List_T terminals, Acceptor_T main_terminal, char *chr,
			long int *tally_matches_low, long int *tally_mismatches_low,
			long int *tally_matches_high, long int *tally_mismatches_high,
			long int *primary_extents, long int *crosshyb_extents,
			IIT_T end_exons_iit, Genomicpos_T chrlength, int mincount_alt,
			Genomicpos_T main_genestart, Genomicpos_T main_geneend,
			Genomicpos_T alt_geneend, int insertlength,
			int readlength, int min_overhang, bool forwardp) {
  List_T p;
  Path_T termpath, path;
  Acceptor_T terminal;
  Intron_T join_point;
  Genomicpos_T join_donor_chrpos, join_acceptor_chrpos;
  Genomicpos_T geneend1, geneend2;

  int n_main_introns, n_alt_introns;
  double main_density, alt_density;

  for (p = terminals; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    join_point = acceptor_used_point(terminal);
    if (join_point == NULL) {
      join_donor_chrpos = 0U;
      join_acceptor_chrpos = 0U;
    } else {
      join_donor_chrpos = join_point->prev_site->chrpos;
      join_acceptor_chrpos = join_point->site->chrpos;
    }

    n_main_introns = acceptor_nintrons(main_terminal,join_donor_chrpos,join_acceptor_chrpos);
    n_alt_introns = acceptor_nintrons(terminal,join_donor_chrpos,join_acceptor_chrpos);
    
    main_density = terminal_density(main_terminal,n_main_introns,
				    main_genestart,main_geneend,insertlength,readlength,min_overhang);
    alt_density = terminal_density(terminal,n_alt_introns,
				   main_genestart,alt_geneend,insertlength,readlength,min_overhang);

    if (alt_density <= 0.05 * main_density) {
      /* Skip */

    } else if (n_main_introns == 1 && n_alt_introns == 1) {
      if (forwardp == true) {
	geneend1 = Cappaths_solve_geneend(main_terminal->chrpos,tally_matches_high,tally_mismatches_high,
					  end_exons_iit,chr,chrlength,forwardp);
	geneend2 = Cappaths_solve_geneend(terminal->chrpos,tally_matches_high,tally_mismatches_high,
					  end_exons_iit,chr,chrlength,forwardp);
      } else {
	geneend1 = Cappaths_solve_geneend(main_terminal->chrpos,tally_matches_low,tally_mismatches_low,
					  end_exons_iit,chr,chrlength,forwardp);
	geneend2 = Cappaths_solve_geneend(terminal->chrpos,tally_matches_low,tally_mismatches_low,
					  end_exons_iit,chr,chrlength,forwardp);
      }
      if (geneend1 == geneend2) {
	/* Alternate splice site */
	Site_get_paths(&path,terminal,chr,mincount_alt,
		       tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		       primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		       /*altpaths_p*/false,/*knownp*/false);
	if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					  main_terminal,n_main_introns,
					  geneend2,/*terminal_type*/TERMINAL_ALT_SITE,
					  /*knownp*/false,/*retainedp*/false)) != NULL) {
	  alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	}
	Path_free(&path);

      } else {
	/* Alternate last exon */
	Site_get_paths(&path,terminal,chr,mincount_alt,
		       tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		       primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		       /*altpaths_p*/false,/*knownp*/false);
	if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					  main_terminal,n_main_introns,
					  geneend2,/*terminal_type*/TERMINAL_ALT_EXON,
					  /*knownp*/false,/*retainedp*/false)) != NULL) {
	  alt_termpaths = List_push(alt_termpaths,(void *) termpath);
	}
	Path_free(&path);
      }

    } else if ((forwardp == true && terminal->chrpos < main_terminal->chrpos) ||
	       (forwardp == false && terminal->chrpos > main_terminal->chrpos)) {
      /* Short form */
      Site_get_paths(&path,terminal,chr,mincount_alt,
		     tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		     primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		     /*altpaths_p*/false,/*knownp*/false);
      if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					main_terminal,n_main_introns,
					path->geneend,/*terminal_type*/TERMINAL_SHORT,
					/*knownp*/false,/*retainedp*/false)) != NULL) {
	alt_termpaths = List_push(alt_termpaths,(void *) termpath);
      }
      Path_free(&path);

    } else {
      /* Long form */
      Site_get_paths(&path,terminal,chr,mincount_alt,
		     tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		     primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		     /*altpaths_p*/false,/*knownp*/false);
      if ((termpath = Path_new_termpath(terminal,n_alt_introns,
					main_terminal,n_main_introns,
					path->geneend,/*terminal_type*/TERMINAL_LONG,
					/*knownp*/false,/*retainedp*/false)) != NULL) {
	alt_termpaths = List_push(alt_termpaths,(void *) termpath);
      }
      Path_free(&path);
    }
  }

  return alt_termpaths;
}


static List_T
get_paths_term_retained (List_T ret_termpaths, List_T terminals, Acceptor_T main_terminal, char *chr,
			 long int *tally_matches_low, long int *tally_mismatches_low,
			 long int *tally_matches_high, long int *tally_mismatches_high,
			 long int *primary_extents, long int *crosshyb_extents,
			 IIT_T end_exons_iit, Genomicpos_T chrlength, int mincount_alt,
			 bool forwardp, bool knownp) {
  List_T p;
  Path_T termpath, path;
  Acceptor_T terminal;
  Intron_T join_point;
  Genomicpos_T join_donor_chrpos, join_acceptor_chrpos;
  Genomicpos_T geneend1, geneend2;
  int n_main_introns, n_alt_introns;

  for (p = terminals; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);

    join_point = terminal->main_intron;
    debug5(printf("Retained terminal is at %s:%u.  Retained terminal intron is at %s:%u..%u\n",
		  chr,terminal->chrpos,chr,join_point->prev_site->chrpos,join_point->site->chrpos));

    /* Previously checked for knownp == true || (terminal->score - basescore) >= short_form_support */
    Site_get_paths(&path,terminal,chr,mincount_alt,
		   tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
		   primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
		   /*altpaths_p*/false,knownp);
    if ((termpath = Path_new_termpath(terminal,/*n_alt_introns*/0,
				      main_terminal,/*n_main_introns*/1,
				      terminal->retained_geneend,/*terminal_type*/TERMINAL_RETAINED,
				      knownp,/*retainedp*/true)) != NULL) {
      ret_termpaths = List_push(ret_termpaths,(void *) termpath);
      if (knownp == true) {
	terminal->known_termpath = termpath;
      }
    }
    Path_free(&path);
  }

  return ret_termpaths;
}



static List_T
Terminals_find_genes_obs (List_T genes, List_T terminals, char *chr, int genei,
			  long int *tally_matches_low, long int *tally_mismatches_low,
			  long int *tally_matches_high, long int *tally_mismatches_high,
			  double *window_diff_low, double *window_diff_high,
			  long int *primary_extents, long int *crosshyb_extents,
			  IIT_T end_exons_iit, int insertlength, int readlength,
			  int min_overhang, Genomicpos_T chrlength,
			  int mincount_alt, bool altpaths_p, bool forwardp) {
  List_T altpaths, alt_termpaths, ret_termpaths, retained_terminals, p;
  Acceptor_T main_terminal;
  Intron_T intron;
  Path_T mainpath, path;
  Genomicpos_T main_genestart, main_geneend, alt_geneend;
  double main_density, alt_density;
  long int retained_coverage;


  /* Main terminal */
  main_terminal = (Acceptor_T) List_head(terminals);
  altpaths = Site_get_paths(&mainpath,main_terminal,chr,mincount_alt,
			    tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
			    primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,altpaths_p,
			    /*knownp*/false);

  if (mainpath == NULL) {
    fprintf(stderr,"#Gene %d had no main path\n",genei);
    return genes;

  } else {
    debug5(printf("Main path\n"));
    debug5(Intron_path_print(mainpath->main_introns));
    main_genestart = mainpath->genestart;
    main_geneend = mainpath->geneend;

    gene_density(&mainpath->density,
#ifdef USE_READ_COUNTS
		 &exon_density_2,bamreader,chr,
#endif
		 mainpath->main_introns,main_genestart,main_geneend,insertlength,readlength,
		 min_overhang,tally_matches_low,tally_mismatches_low,
		 tally_matches_high,tally_mismatches_high);


    /* Retained introns and terminals */
    retained_terminals = (List_T) NULL;
    for (p = mainpath->main_introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->known_retained_p == true) {
	/* Already found to be retained */
      } else if (forwardp == true) {
	if ((alt_geneend = flat_chrpos_high(window_diff_high,/*low*/intron->donorpos,/*high*/intron->acceptorpos,
					    /*allowed_zeroes*/20)) >= intron->acceptorpos) {
	  debug5(printf("Fwd intron at %s:%u..%u is retained novel because geneend is %u\n",
			chr,intron->donorpos,intron->acceptorpos,alt_geneend));
	  /* main_density = Splice_count(intron->splice); */

	  retained_coverage = tally_mincount(intron->donorpos,intron->acceptorpos,
					     tally_matches_low,tally_mismatches_low,
					     tally_matches_high,tally_mismatches_high);
	  alt_density = (double) retained_coverage/(double) (readlength + readlength);

	  if (alt_density > 0.05 * mainpath->density) {
	    path = Path_new_retained_intron(intron,forwardp,/*knownp*/false);
	    Path_add_genebounds(path,tally_matches_low,tally_mismatches_low,
				tally_matches_high,tally_mismatches_high,
				end_exons_iit,chr,chrlength,forwardp);
	    List_free(&path->alt_introns);
	    path->alt_introns = (List_T) NULL;
	    altpaths = List_push(altpaths,(void *) path);
	  }

	} else if (intron->site->retained_geneend == 0U &&
		   alt_geneend > intron->donorpos + RETAINED_BUFFER &&
		   alt_geneend < intron->acceptorpos - RETAINED_BUFFER) {
	  debug5(printf("Intron at %s:%u..%u has donor coverage novel from geneend %u\n",
			chr,intron->donorpos,intron->acceptorpos,alt_geneend));
	  intron->site->retained_geneend = alt_geneend;
	  retained_terminals = List_push(retained_terminals,(void *) intron->site);
	}

      } else {
	if ((alt_geneend = flat_chrpos_low(window_diff_low,/*low*/intron->acceptorpos,/*high*/intron->donorpos,
				       /*allowed_zeroes*/20)) <= intron->donorpos) {
	  debug5(printf("Rev intron at %s:%u..%u is retained novel because geneend is %u\n",
			chr,intron->donorpos,intron->acceptorpos,alt_geneend));
	  /* main_density = Splice_count(intron->splice); */

	  retained_coverage = tally_mincount(intron->donorpos,intron->acceptorpos,
					     tally_matches_low,tally_mismatches_low,
					     tally_matches_high,tally_mismatches_high);
	  alt_density = (double) retained_coverage/(double) (readlength + readlength);

	  if (alt_density > 0.05 * mainpath->density) {
	    path = Path_new_retained_intron(intron,forwardp,/*knownp*/false);
	    Path_add_genebounds(path,tally_matches_low,tally_mismatches_low,
				tally_matches_high,tally_mismatches_high,
				end_exons_iit,chr,chrlength,forwardp);
	    List_free(&path->alt_introns);
	    path->alt_introns = (List_T) NULL;
	    altpaths = List_push(altpaths,(void *) path);
	  }

	} else if (intron->site->retained_geneend == 0U &&
		   alt_geneend < intron->donorpos - RETAINED_BUFFER &&
		   alt_geneend > intron->acceptorpos + RETAINED_BUFFER) {
	  debug5(printf("Intron at %s:%u..%u has donor coverage novel from geneend %u\n",
			chr,intron->donorpos,intron->acceptorpos,alt_geneend));
	  intron->site->retained_geneend = alt_geneend;
	  retained_terminals = List_push(retained_terminals,(void *) intron->site);
	}
      }
    }

    /* Alternate terminals */
    alt_termpaths = get_paths_term_alt_obs(/*alt_termpaths*/NULL,List_next(terminals),main_terminal,chr,
					   tally_matches_low,tally_mismatches_low,
					   tally_matches_high,tally_mismatches_high,
					   primary_extents,crosshyb_extents,end_exons_iit,
					   chrlength,mincount_alt,
					   main_genestart,main_geneend,alt_geneend,
					   insertlength,readlength,min_overhang,forwardp);

    ret_termpaths = NULL;
    for (p = main_terminal->known_ret_termpaths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      ret_termpaths = List_push(ret_termpaths,(void *) Path_copy(path));
    }
    ret_termpaths = get_paths_term_retained(ret_termpaths,retained_terminals,main_terminal,chr,
					    tally_matches_low,tally_mismatches_low,
					    tally_matches_high,tally_mismatches_high,
					    primary_extents,crosshyb_extents,end_exons_iit,
					    chrlength,mincount_alt,forwardp,/*knownp*/false);
    List_free(&retained_terminals);

    return List_push(genes,(void *) Gene_new(mainpath,altpaths,alt_termpaths,ret_termpaths,forwardp));
  }
}


static Genomicpos_T
donor_coverage_known (IIT_T middle_exons_iit, IIT_T end_exons_iit, char *chr,
		      Genomicpos_T donorpos, Genomicpos_T acceptorpos, bool forwardp) {
  int *matches, nmatches, i;
  Interval_T interval;
  Genomicpos_T farthestpos;

  if (forwardp == true) {
    farthestpos = 0U;
    matches = IIT_get(&nmatches,middle_exons_iit,chr,donorpos+1U,donorpos+1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      interval = IIT_interval(middle_exons_iit,matches[i]);
      if (Interval_sign(interval) > 0 && Interval_high(interval) > farthestpos) {
	farthestpos = Interval_high(interval);
      }
    }
    FREE(matches);

    matches = IIT_get(&nmatches,end_exons_iit,chr,donorpos+1U,donorpos+1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      interval = IIT_interval(end_exons_iit,matches[i]);
      if (Interval_sign(interval) > 0 && Interval_high(interval) > farthestpos) {
	farthestpos = Interval_high(interval);
      }
    }
    FREE(matches);

    if (farthestpos > donorpos + RETAINED_BUFFER && farthestpos < acceptorpos - RETAINED_BUFFER) {
      return farthestpos;
    } else {
      return 0U;
    }

  } else {
    farthestpos = -1U;
    matches = IIT_get(&nmatches,middle_exons_iit,chr,donorpos-1U,donorpos-1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      interval = IIT_interval(middle_exons_iit,matches[i]);
      if (Interval_sign(interval) < 0 && Interval_low(interval) < farthestpos) {
	farthestpos = Interval_low(interval);
      }
    }
    FREE(matches);

    matches = IIT_get(&nmatches,end_exons_iit,chr,donorpos-1U,donorpos-1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      interval = IIT_interval(end_exons_iit,matches[i]);
      if (Interval_sign(interval) < 0 && Interval_low(interval) < farthestpos) {
	farthestpos = Interval_low(interval);
      }
    }
    FREE(matches);

    if (farthestpos < donorpos - RETAINED_BUFFER && farthestpos > acceptorpos + RETAINED_BUFFER) {
      return farthestpos;
    } else {
      return 0U;
    }
  }
}



static List_T
Terminals_find_genes_known (List_T genes, List_T terminals, char *chr, int genei,
			    long int *tally_matches_low, long int *tally_mismatches_low,
			    long int *tally_matches_high, long int *tally_mismatches_high,
			    long int *primary_extents, long int *crosshyb_extents,
			    IIT_T middle_exons_iit, IIT_T end_exons_iit, Genomicpos_T chrlength,
			    int mincount_alt, bool altpaths_p, bool forwardp) {
  List_T altpaths, alt_termpaths, ret_termpaths, retained_terminals = NULL, p;
  Acceptor_T main_terminal;
  Intron_T intron;
  Path_T mainpath, termpath, path;
  Genomicpos_T main_genestart, main_geneend, alt_geneend;

  /* Main terminal */
  main_terminal = (Acceptor_T) List_head(terminals);
  altpaths = Site_get_paths(&mainpath,main_terminal,chr,mincount_alt,
			    tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
			    primary_extents,crosshyb_extents,end_exons_iit,chrlength,forwardp,
			    altpaths_p,/*knownp*/true);

  if (mainpath == NULL) {
    fprintf(stderr,"#Gene %d had no main path\n",genei);
    return genes;

  } else {
    debug5(printf("Main path\n"));
    debug5(Intron_path_print(mainpath->main_introns));

    main_genestart = mainpath->genestart;
    main_geneend = mainpath->geneend;

    /* Retained introns and terminals */
    for (p = mainpath->main_introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (IIT_contained(middle_exons_iit,chr,intron->low+1U,intron->high-1U) == true ||
	  IIT_contained(end_exons_iit,chr,intron->low+1U,intron->high-1U) == true) {
	debug5(printf("Intron at %s:%u..%u is retained known\n",chr,intron->donorpos,intron->acceptorpos));
	path = Path_new_retained_intron(intron,forwardp,/*knownp*/true);
	Path_add_genebounds(path,tally_matches_low,tally_mismatches_low,
			    tally_matches_high,tally_mismatches_high,
			    end_exons_iit,chr,chrlength,forwardp);
	List_free(&path->alt_introns);
	path->alt_introns = (List_T) NULL;
	altpaths = List_push(altpaths,(void *) path);
	intron->known_retained_p = true;
      } else if ((alt_geneend = donor_coverage_known(middle_exons_iit,end_exons_iit,chr,intron->donorpos,intron->acceptorpos,forwardp)) > 0U) {
	debug5(printf("Intron at %s:%u..%u has donor coverage known\n",chr,intron->donorpos,intron->acceptorpos));
	intron->site->retained_geneend = alt_geneend;
	retained_terminals = List_push(retained_terminals,(void *) intron->site);
      }
    }

    /* Alternate terminals */
    alt_termpaths = get_paths_term_alt_known(/*alt_termpaths*/NULL,List_next(terminals),main_terminal,chr,
					     tally_matches_low,tally_mismatches_low,
					     tally_matches_high,tally_mismatches_high,
					     primary_extents,crosshyb_extents,end_exons_iit,
					     chrlength,mincount_alt,forwardp);
    ret_termpaths = get_paths_term_retained(/*ret_termpaths*/NULL,retained_terminals,main_terminal,chr,
					    tally_matches_low,tally_mismatches_low,
					    tally_matches_high,tally_mismatches_high,
					    primary_extents,crosshyb_extents,end_exons_iit,
					    chrlength,mincount_alt,forwardp,/*knownp*/true);
    List_free(&retained_terminals);

    main_terminal->known_mainpath = mainpath;
    main_terminal->known_altpaths = altpaths;
    main_terminal->known_alt_termpaths = alt_termpaths;
    main_terminal->known_ret_termpaths = ret_termpaths;

    return List_push(genes,(void *) Gene_new(mainpath,altpaths,alt_termpaths,ret_termpaths,forwardp));
  }
}




static List_T
Terminals_gather_fwd_list (Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd) {
  List_T survivors = NULL, list = NULL;
  Acceptor_T *terminals, terminal, site;
  int nterminals, i;

#if 0
  /* Fix cases where a donor is terminal, meaning it didn't complete */
  for (i = 0; i < nsites_donor_fwd; i++) {
    donor = sites_donor_fwd[i];
    if (donor->terminalp == true) {
      donor->terminalp = false;
      if (donor->prev != NULL) {
	donor->prev->terminalp = true;
      }
    }
  }
#endif

  for (i = 0; i < nsites_acceptor_fwd; i++) {
    site = sites_acceptor_fwd[i];
    if (site->terminalp == true && site->prev != NULL) {
      list = List_push(list,site);
    }
  }

  nterminals = List_length(list);
  terminals = (Acceptor_T *) List_to_array(list,NULL);
  qsort(terminals,nterminals,sizeof(Acceptor_T),Acceptor_score_cmp);
  List_free(&list);

  for (i = 0; i < nterminals; i++) {
    terminal = terminals[i];
    if (acceptor_used_p(terminal) == 0U) {
      /* printf("Terminal %u is a survivor\n",terminal->chrpos); */
      survivors = List_push(survivors,terminal);
      acceptor_mark_used(terminal,terminal->chrpos);
    } else {
      /* printf("Terminal %u is not a survivor\n",terminal->chrpos); */
      /* Might speed up process */
      /* path_mark_used(terminal); */
    }
  }
    
  FREE(terminals);
  return List_reverse(survivors);
}

static List_T
Terminals_gather_rev_list (Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev) {
  List_T survivors = NULL, list = NULL;
  Acceptor_T *terminals, terminal, site;
  int nterminals, i;

#if 0
  /* Fix cases where a donor is terminal, meaning it didn't complete */
  for (i = 0; i < nsites_donor_rev; i++) {
    donor = sites_donor_rev[i];
    if (donor->terminalp == true) {
      donor->terminalp = false;
      if (donor->prev != NULL) {
	donor->prev->terminalp = true;
      }
    }
  }
#endif

  for (i = 0; i < nsites_acceptor_rev; i++) {
    site = sites_acceptor_rev[i];
    if (site->terminalp == true && site->prev != NULL) {
      list = List_push(list,site);
    }
  }

  nterminals = List_length(list);
  terminals = (Acceptor_T *) List_to_array(list,NULL);
  qsort(terminals,nterminals,sizeof(Acceptor_T),Acceptor_score_cmp);
  List_free(&list);
  
  for (i = 0; i < nterminals; i++) {
    terminal = terminals[i];
    if (acceptor_used_p(terminal) == 0U) {
      survivors = List_push(survivors,terminal);
      acceptor_mark_used(terminal,terminal->chrpos);
    } else {
      /* Might speed up process */
      /* path_mark_used(terminal); */
    }
  }

  FREE(terminals);
  return List_reverse(survivors);
}



static Uinttable_T
Terminals_gather_table (Donor_T *sites_donor, int nsites_donor,
			Acceptor_T *sites_acceptor, int nsites_acceptor) {
  Uinttable_T terminal_table;
  List_T list = NULL;
  Acceptor_T *terminals, terminal, site;
  Donor_T donor;
  Genomicpos_T terminal_chrpos;
  int nterminals, i;

#if 0
  /* Fix cases where a donor is terminal, meaning it didn't complete */
  for (i = 0; i < nsites_donor; i++) {
    donor = sites_donor[i];
    if (donor->terminalp == true) {
      donor->terminalp = false;
      if (donor->prev != NULL) {
	donor->prev->terminalp = true;
      }
    }
  }
#endif

  for (i = 0; i < nsites_donor; i++) {
    donor = sites_donor[i];
    donor->usedp = 0U;
  }

  for (i = 0; i < nsites_acceptor; i++) {
    site = sites_acceptor[i];
    site->usedp = 0U;
    if (site->terminalp == false) {
      debug3(printf("Acceptor %u is not a terminal\n",site->chrpos));
    } else if (site->prev == NULL) {
      debug3(printf("Acceptor %u is a terminal, but has no prev\n",site->chrpos));
    } else if (site->main_intron->validp == false) {
      /* Shouldn't reach here */
      debug3(printf("Acceptor %u is a terminal, but main intron is not valid\n",site->chrpos));
    } else if (Splice_valid_end_p(site->main_intron->splice) == false) {
      debug3(printf("Acceptor %u is a terminal, but main intron is not valid for end\n",site->chrpos));
    } else {
      debug3(printf("Acceptor %u is a valid terminal\n",site->chrpos));
      list = List_push(list,site);
    }
  }

  nterminals = List_length(list);
  terminals = (Acceptor_T *) List_to_array(list,NULL);
  qsort(terminals,nterminals,sizeof(Acceptor_T),Acceptor_score_cmp);
  List_free(&list);
  
  terminal_table = Uinttable_new(nterminals);

  for (i = 0; i < nterminals; i++) {
    terminal = terminals[i];
    if ((terminal_chrpos = acceptor_used_p(terminal)) == 0U) {
      debug3(printf("Terminal %u is a survivor, so marking it as usedp\n",
		    terminal->chrpos));
      Uinttable_put(terminal_table,terminal->chrpos,(void *) List_push(NULL,terminal));
      acceptor_mark_used(terminal,terminal->chrpos);

    } else if (terminal->known_terminalp == true) {
      debug3(printf("Terminal %u is a known terminal, so marking it as usedp\n",
		    terminal->chrpos));
      Uinttable_put(terminal_table,terminal->chrpos,(void *) List_push(NULL,terminal));
      acceptor_mark_used(terminal,terminal->chrpos);

    } else {
      debug3(printf("Terminal %u is not a survivor, so grouping it with survivor at %u\n",
		    terminal->chrpos,terminal_chrpos));
      Uinttable_put(terminal_table,terminal_chrpos,
		    (void *) List_push((List_T) Uinttable_get(terminal_table,terminal_chrpos),
				       terminal));
      /* Might speed up process */
      /* path_mark_used(terminal); */
    }
  }
    
  FREE(terminals);
  return terminal_table;
}



#if 0
static List_T
Terminals_sort (List_T list) {
  List_T sorted = NULL;
  Acceptor_T *terminals;
  int nterminals, i;

  nterminals = List_length(list);
  terminals = (Acceptor_T *) List_to_array(list,NULL);
  qsort(terminals,nterminals,sizeof(Acceptor_T),Acceptor_position_cmp);
  List_free(&list);

  sorted = (List_T) NULL;
  for (i = nterminals-1; i >= 0; i--) {
    sorted = List_push(sorted,terminals[i]);
  }
  FREE(terminals);

  return sorted;
}
#endif


/* Returns the introns in high-to-low order, i.e., the head of the
   list has the highest position */
static List_T
Site_get_introns (Acceptor_T acceptor) {
  Donor_T donor;
  Intron_T intron;
  List_T p;

  if (acceptor == NULL) {
    return (List_T) NULL;

  } else if ((donor = acceptor->prev) == NULL) {
    fprintf(stderr,"Singleton down at %u\n",acceptor->chrpos);
    abort();

  } else {
    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->prev_site == donor) {
	if (intron->validp == false) {
	  /* Recursive call.  Intron can be invalid due to trimming. */
	  return Site_get_introns(donor->prev);
	} else {
	  /* Recursive call */
	  return List_push(Site_get_introns(donor->prev),intron);
	}
      }
    }
    fprintf(stderr,"At acceptor %u, intron not found for donor at %u\n",acceptor->chrpos,donor->chrpos);
    abort();
  }
}


/************************************************************************
 *   Conflict resolution
 ************************************************************************/

static List_T
Terminals_get_introns (List_T terminals) {
  List_T p;
  List_T introns = NULL;
  Acceptor_T terminal;

  for (p = terminals; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    introns = List_append(Site_get_introns(terminal),introns);
  }

  return introns;
}


static int
find_conflicting_introns (List_T terminals_fwd, List_T terminals_rev, char *chr) {
  int nconflicts = 0;
  List_T intronlist, intronlist_fwd, intronlist_rev;
  Intron_T *introns;
  int nintrons, i, j;

  intronlist_fwd = Terminals_get_introns(terminals_fwd);
  intronlist_rev = Terminals_get_introns(terminals_rev);
  intronlist = List_append(intronlist_fwd,intronlist_rev);

  if ((nintrons = List_length(intronlist)) == 0) {
    return 0;

  } else {
    introns = (Intron_T *) List_to_array(intronlist,NULL);
    List_free(&intronlist);

    qsort(introns,nintrons,sizeof(Intron_T),Intron_low_cmp);
    for (i = 0; i < nintrons; i++) {
      /* Previously, was introns[j]->low < introns[i]->low + 20, but now checking all overlaps */
      for (j = i+1; j < nintrons && introns[j]->low < introns[i]->high + 20; j++) {

	if (introns[j]->low < introns[i]->low + 20 && introns[j]->low < introns[i]->high - 20) {
	  /* Endpoints do not conflict */
	} else if (introns[j]->forwardp == introns[i]->forwardp) {
	  /* No conflict */
	} else {
	  fprintf(stderr,"Low position conflict between %s:%u..%u (count %d, pathscore %d) and %s:%u..%u (count %d, pathscore %d)",
		  chr,introns[i]->low,introns[i]->high,Splice_count(introns[i]->splice),introns[i]->pathscore,
		  chr,introns[j]->low,introns[j]->high,Splice_count(introns[j]->splice),introns[j]->pathscore);
	  if (introns[i]->pathscore < introns[j]->pathscore) {
	    fprintf(stderr," => eliminate first\n");
	    introns[i]->validp = false;
	    nconflicts++;
	  } else if (introns[j]->pathscore < introns[i]->pathscore) {
	    fprintf(stderr," => eliminate second\n");
	    introns[j]->validp = false;
	    nconflicts++;
	  } else if (Splice_count(introns[i]->splice) < Splice_count(introns[j]->splice)) {
	    fprintf(stderr," => eliminate first\n");
	    introns[i]->validp = false;
	    nconflicts++;
	  } else if (Splice_count(introns[j]->splice) < Splice_count(introns[i]->splice)) {
	    fprintf(stderr," => eliminate second\n");
	    introns[j]->validp = false;
	    nconflicts++;
	  } else {
	    fprintf(stderr," => eliminate neither\n");
	  }
	}
      }
    }

    qsort(introns,nintrons,sizeof(Intron_T),Intron_high_cmp);
    for (i = 0; i < nintrons; i++) {
      /* Previously, was introns[j]->high > introns[i]->high - 20, but now checking all overlaps */
      for (j = i+1; j < nintrons && introns[j]->high > introns[i]->low - 20; j++) {

	if (introns[j]->high < introns[i]->high - 20 && introns[j]->high > introns[i]->low + 20) {
	  /* Endpoints do not conflict */
	} else if (introns[j]->forwardp == introns[i]->forwardp) {
	  /* No conflict */
	} else {
	  fprintf(stderr,"High position conflict between %s:%u..%u (count %d, pathscore %d) and %s:%u..%u (count %d, pathscore %d)",
		  chr,introns[i]->low,introns[i]->high,Splice_count(introns[i]->splice),introns[i]->pathscore,
		  chr,introns[j]->low,introns[j]->high,Splice_count(introns[j]->splice),introns[j]->pathscore);
	  if (introns[i]->pathscore < introns[j]->pathscore) {
	    fprintf(stderr," => eliminate first\n");
	    introns[i]->validp = false;
	    nconflicts++;
	  } else if (introns[j]->pathscore < introns[i]->pathscore) {
	    fprintf(stderr," => eliminate second\n");
	    introns[j]->validp = false;
	    nconflicts++;
	  } else if (Splice_count(introns[i]->splice) < Splice_count(introns[j]->splice)) {
	    fprintf(stderr," => eliminate first\n");
	    introns[i]->validp = false;
	    nconflicts++;
	  } else if (Splice_count(introns[j]->splice) < Splice_count(introns[i]->splice)) {
	    fprintf(stderr," => eliminate second\n");
	    introns[j]->validp = false;
	    nconflicts++;
	  } else {
	    fprintf(stderr," => eliminate neither\n");
	  }
	}
      }
    }

    FREE(introns);
    return nconflicts;
  }
}


#if 0
/* Should really compute this, rather than just read down->score,
   because path may have been trimmed.  However, this adds up only the
   intron scores. */
static Score_T
Acceptor_get_score (Acceptor_T acceptor) {
  Donor_T donor;
  Intron_T intron;
  List_T p;

  if (acceptor == NULL) {
    return 0;

  } else if ((donor = acceptor->prev) == NULL) {
    fprintf(stderr,"Singleton acceptor at %u\n",acceptor->chrpos);
    abort();

  } else {
    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->prev_site == donor) {
	if (intron->validp == false) {
	  /* Intron can be invalid due to trimming. */
	  return Site_get_score(donor->prev);
	} else {
	  return Splice_count(intron->splice) + Site_get_score(donor->prev);
	}
      }
    }
    fprintf(stderr,"At acceptor %u, intron not found for donor at %u\n",acceptor->chrpos,donor->chrpos);
    abort();
  }
}
#endif



static void
assign_pathscore_to_introns (Acceptor_T acceptor, Score_T pathscore) {
  Donor_T donor;
  Intron_T intron;
  List_T p;

  if (acceptor == NULL) {
    return;

  } else if ((donor = acceptor->prev) == NULL) {
    fprintf(stderr,"Singleton acceptor at %u\n",acceptor->chrpos);
    abort();

  } else {
    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->prev_site == donor) {
	if (intron->validp == true) {
	  /* Intron can be invalid due to trimming. */
	  intron->pathscore = pathscore;
	}
      }
    }
    assign_pathscore_to_introns(donor->prev,pathscore);
  }
}



/************************************************************************
 *   Output
 *   For output, we can assume that paths contain only valid introns
 ************************************************************************/


#if 0
static void
Site_print_fwd (FILE *fp, Acceptor_T acceptor, char *chr) {
  Donor_T donor;
  Intron_T intron;
  List_T p;

  if (acceptor == NULL) {
    return;

  } else if ((donor = acceptor->prev) == NULL) {
    fprintf(stderr,"Singleton acceptor at %u\n",acceptor->chrpos);
    abort();

  } else {
    Site_print_fwd(fp,donor->prev,chr);

    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->prev_site == donor) {
	if (intron->validp == false) {
	  /* Intron can be invalid due to trimming, although should have been resolved before printing */
	  fprintf(stderr,"Got an invalid intron at %s:%u..%u at print time\n",chr,donor->chrpos,acceptor->chrpos);
	  abort();
	} else {
	  fprintf(fp,">%d %s:%u..%u\n",Splice_count(intron->splice),chr,donor->chrpos,acceptor->chrpos);
	}
      }
    }

    return;
  }
}


static void
Site_print_rev (FILE *fp, Acceptor_T acceptor, char *chr) {
  Donor_T donor;
  Intron_T intron;
  List_T p;

  if (acceptor == NULL) {
    return;

  } else if ((donor = acceptor->prev) == NULL) {
    fprintf(stderr,"Singleton acceptor at %u\n",acceptor->chrpos);
    abort();

  } else {
    Site_print_rev(fp,donor->prev,chr);

    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      if (intron->prev_site == donor) {
	if (intron->validp == false) {
	  /* Intron can be invalid due to trimming, although should have been resolved before printing */
	  fprintf(stderr,"Got an invalid intron at %s:%u..%u at print time\n",chr,donor->chrpos,acceptor->chrpos);
	  abort();
	} else {
	  /* Intron can be invalid due to trimming, although should have been resolved before printing */
	  fprintf(fp,">%d %s:%u..%u\n",Splice_count(intron->splice),chr,donor->chrpos,acceptor->chrpos);
	}
      }
    }

    return;
  }
}
#endif

    
static void
compute_cum_int (long int *cum, long int *x, long int *y, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  cum[0] = 0.0;
  for (chrpos = 1; chrpos <= chrlength; chrpos++) {
    cum[chrpos] = x[chrpos] + y[chrpos] + cum[chrpos - 1U];
  }
  return;
}

static void
compute_log_tally (double *log_tally, long int *x, long int *y, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  /* For non-cum results, want < and not <= */
  for (chrpos = 0; chrpos < chrlength; chrpos++) {
    log_tally[chrpos] = log((double) (x[chrpos]+y[chrpos]+1));
  }
  
  return;
}

static void
compute_cum_double (double *cum, double *x, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  cum[0] = 0.0;
  for (chrpos = 1; chrpos <= chrlength; chrpos++) {
    cum[chrpos] = x[chrpos] + cum[chrpos - 1U];
  }
  return;
}

static void
compute_diff (double *window_diff, double *cumlog_tally, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;
  double sum_right, sum_left;

  for (chrpos = 30U + 1; chrpos <= chrlength - 30U - 1; chrpos++) {
    sum_right = cumlog_tally[chrpos + 30U - 1] - cumlog_tally[chrpos - 1];
    sum_left = cumlog_tally[chrpos - 1] - cumlog_tally[chrpos - 30U - 1];
    window_diff[(int) chrpos] = sum_right - sum_left;
    /* window_sumx[(int) chrpos] - window_sumx[(int) (chrpos - 30U)]; */
    /* debug3(printf("Putting diff %.1f into chrpos %u\n",window_diff[(int) chrpos],chrpos)); */
  }

  return;
}


/************************************************************************/


/************************************************************************/


#if 0
static Genomicpos_T
distance_diff (Genomicpos_T x, Genomicpos_T y) {
  if (x < y) {
    return y - x;
  } else {
    return x - y;
  }
}
#endif


static bool
identicalp (Genomicpos_T low1, Genomicpos_T high1, Genomicpos_T low2, Genomicpos_T high2) {
  if (low1 == low2 && high1 == high2) {
    return true;
  } else {
    return false;
  }
}


static int
splice_match (IIT_T splices_iit, char *chr, Genomicpos_T intron_low, Genomicpos_T intron_high, int sign) {
  int nsplices = 0;
  int *matches, nmatches, i;
  Interval_T interval;
  char *label;
  bool allocp;

  matches = IIT_get(&nmatches,splices_iit,/*divstring*/chr,
		    /*coordstart*/intron_low,/*coordend*/intron_high,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    interval = IIT_interval(splices_iit,matches[i]);
    if (Interval_sign(interval) == sign) {
      if (identicalp(intron_low,intron_high,Interval_low(interval),Interval_high(interval)) == true) {
	label = IIT_label(splices_iit,matches[i],&allocp);
	if (sscanf(label,"%d",&nsplices) != 1) {
	  fprintf(stderr,"Could not parse label %s\n",label);
	}
	if (allocp == true) {
	  FREE(label);
	}
	FREE(matches);
	return nsplices;
      }
    }
  }

  FREE(matches);
  return 0;
}



static void
put_splice_sites (Uinttable_T donor_sitestable, Uinttable_T acceptor_sitestable,
		  List_T obs_splices, List_T knowngenes, bool forwardp) {
  List_T g, p, q;
  Gene_T gene;
  Path_T path;
  Intron_T intron;
  Donor_T donor;
  Acceptor_T acceptor;
  Splice_T splice, knownsplice;
  Genomicpos_T donorpos, acceptorpos, known_prev_acceptorpos;
  int nsplices, sign;

  nsplices = List_length(obs_splices);
  if (forwardp == true) {
    sign = +1;
  } else {
    sign = -1;
  }


  /* Handle known genes first */
  for (g = knowngenes; g != NULL; g = List_next(g)) {
    gene = (Gene_T) List_head(g);
    if (gene->forwardp == forwardp) {
      path = gene->mainpath;

      known_prev_acceptorpos = 0U;
      for (p = path->main_introns; p != NULL; p = List_next(p)) {
	intron = (Intron_T) List_head(p);
	donorpos = intron->donorpos;
	acceptorpos = intron->acceptorpos;

	if ((donor = (Donor_T) Uinttable_get(donor_sitestable,donorpos)) == NULL) {
	  Uinttable_put(donor_sitestable,donorpos,(void *) Donor_new(donorpos,known_prev_acceptorpos));
	}

	if ((acceptor = (Acceptor_T) Uinttable_get(acceptor_sitestable,acceptorpos)) == NULL) {
	  acceptor = Acceptor_new(acceptorpos);
	  Uinttable_put(acceptor_sitestable,acceptorpos,(void *) acceptor);
	}
	if (List_next(p) == NULL) {
	  /* For last intron, set known_terminalp to be true */
	  acceptor->known_terminalp = true;
	}

	if (Uinttable_get(acceptor->splices_table,donorpos) == NULL) {
	  splice = Splice_new_known(donorpos,acceptorpos,sign,/*mainpath_p*/true);
	  Uinttable_put(acceptor->splices_table,donorpos,(void *) splice);
	}
	known_prev_acceptorpos = acceptorpos;
      }

      for (q = gene->altpaths; q != NULL; q = List_next(q)) {
	path = (Path_T) List_head(q);
	for (p = path->alt_introns; p != NULL; p = List_next(p)) {
	  intron = (Intron_T) List_head(p);
	  donorpos = intron->donorpos;
	  acceptorpos = intron->acceptorpos;

	  if ((donor = (Donor_T) Uinttable_get(donor_sitestable,donorpos)) == NULL) {
	    Uinttable_put(donor_sitestable,donorpos,(void *) Donor_new(donorpos,/*known_prev_acceptorpos*/0U));
	  }

	  if ((acceptor = (Acceptor_T) Uinttable_get(acceptor_sitestable,acceptorpos)) == NULL) {
	    acceptor = Acceptor_new(acceptorpos);
	    Uinttable_put(acceptor_sitestable,acceptorpos,(void *) acceptor);
	  }
	  if (Uinttable_get(acceptor->splices_table,donorpos) == NULL) {
	    splice = Splice_new_known(donorpos,acceptorpos,sign,/*mainpath_p*/false);
	    Uinttable_put(acceptor->splices_table,donorpos,(void *) splice);
	  }
	}
      }

      for (q = gene->alt_termpaths; q != NULL; q = List_next(q)) {
	path = (Path_T) List_head(q);
	for (p = path->alt_introns; p != NULL; p = List_next(p)) {
	  intron = (Intron_T) List_head(p);
	  donorpos = intron->donorpos;
	  acceptorpos = intron->acceptorpos;

	  if ((donor = (Donor_T) Uinttable_get(donor_sitestable,donorpos)) == NULL) {
	    Uinttable_put(donor_sitestable,donorpos,(void *) Donor_new(donorpos,/*known_prev_acceptorpos*/0U));
	  }

	  if ((acceptor = (Acceptor_T) Uinttable_get(acceptor_sitestable,acceptorpos)) == NULL) {
	    acceptor = Acceptor_new(acceptorpos);
	    Uinttable_put(acceptor_sitestable,acceptorpos,(void *) acceptor);
	  }
	  if (Uinttable_get(acceptor->splices_table,donorpos) == NULL) {
	    splice = Splice_new_known(donorpos,acceptorpos,sign,/*mainpath_p*/false);
	    Uinttable_put(acceptor->splices_table,donorpos,(void *) splice);
	  }
	}
      }
    }
  }


  for (p = obs_splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (Splice_sign(splice) == sign) {
      donorpos = Splice_donorpos(splice);
      acceptorpos = Splice_acceptorpos(splice);

      if ((donor = (Donor_T) Uinttable_get(donor_sitestable,donorpos)) == NULL) {
	Uinttable_put(donor_sitestable,donorpos,(void *) Donor_new(donorpos,/*known_prev_acceptorpos*/0U));
      }

      if ((acceptor = (Acceptor_T) Uinttable_get(acceptor_sitestable,acceptorpos)) == NULL) {
	acceptor = Acceptor_new(acceptorpos);
	Uinttable_put(acceptor_sitestable,acceptorpos,(void *) acceptor);
      }
      if ((knownsplice = Uinttable_get(acceptor->splices_table,donorpos)) != NULL) {
	debug(printf("Saw known fwd splice at %u..%u\n",donorpos,acceptorpos));
	Splice_transfer_info(knownsplice,splice);
      } else {
	debug(printf("Saw novel fwd splice at %u..%u\n",donorpos,acceptorpos));
	Uinttable_put(acceptor->splices_table,donorpos,(void *) splice);
      }
    }
  }

  return;
}


static void
free_donor_sites (Donor_T *sites, int nsites) {
  int i;

  for (i = 0; i < nsites; i++) {
    Donor_free(&(sites[i]));
  }
  FREE(sites);
  return;
}

static void
free_acceptor_sites (Acceptor_T *sites, int nsites) {
  int i;

  for (i = 0; i < nsites; i++) {
    Acceptor_free(&(sites[i]));
  }
  FREE(sites);
  return;
}



/************************************************************************
 *   Debugging procedures
 ************************************************************************/


static void
dump_sites (Donor_T *sites_donor, int nsites_donor,
	    Acceptor_T *sites_acceptor, int nsites_acceptor,
	    char *chr) {
  Donor_T donor;
  Acceptor_T acceptor;
  int i;

#if 0
  for (i = 0; i < nsites_iblock; i++) {
    site = sites_iblock[i];
    printf(">I %s:%u\n",chr,site->chrpos);
  }
#endif

  for (i = 0; i < nsites_donor; i++) {
    donor = sites_donor[i];
    printf(">I %s:%u up\n",chr,donor->chrpos);
  }

  for (i = 0; i < nsites_acceptor; i++) {
    acceptor = sites_acceptor[i];
    printf(">I %s:%u down\n",chr,acceptor->chrpos);
  }

  return;
}


static void
dump_exons (Donor_T *sites_donor, int nsites_donor, char *chr) {
  List_T p;
  Donor_T donor;
  Exon_T exon;
  int i;

  for (i = 0; i < nsites_donor; i++) {
    donor = sites_donor[i];
    for (p = donor->exons; p != NULL; p = List_next(p)) {
      exon = (Exon_T) List_head(p);
      printf(">%d %s:%u..%u\n",exon->incrscore,chr,exon->prev_site->chrpos,donor->chrpos);
    }
  }

  return;
}


static void
dump_introns (Acceptor_T *sites_acceptor, int nsites_acceptor, char *chr) {
  List_T p;
  Acceptor_T acceptor;
  Intron_T intron;
  int i;

  for (i = 0; i < nsites_acceptor; i++) {
    acceptor = sites_acceptor[i];
    for (p = acceptor->introns; p != NULL; p = List_next(p)) {
      intron = (Intron_T) List_head(p);
      printf(">%d/%d %s:%u..%u\n",
	     Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),
	     chr,intron->prev_site->chrpos,acceptor->chrpos);
    }
  }

  return;
}



static void
dump_donor_fwd (Donor_T site) {
  Exon_T exon;
  List_T p;

  printf("Fwd donor at %u with score %d, and exons:\n",site->chrpos,site->score);
  for (p = site->exons; p != NULL; p = List_next(p)) {
    exon = (Exon_T) List_head(p);
    printf("  to %u, incrscore %d",exon->prev_site->chrpos,exon->incrscore);
    if (exon->prev_site == site->prev) {
      printf(" **");
    }
    printf("\n");
  }
  printf("\n");
  return;
}

static void
dump_donor_rev (Donor_T site) {
  Exon_T exon;
  List_T p;

  printf("Rev donor at %u with score %d, and exons:\n",site->chrpos,site->score);
  for (p = site->exons; p != NULL; p = List_next(p)) {
    exon = (Exon_T) List_head(p);
    printf("  to %u, incrscore %d",exon->prev_site->chrpos,exon->incrscore);
    if (exon->prev_site == site->prev) {
      printf(" **");
    }
    printf("\n");
  }
  printf("\n");
  return;
}

static void
dump_acceptor_fwd (Acceptor_T site) {
  Intron_T intron;
  Acceptor_alt_T alt;
  List_T p, q;

  printf("Fwd acceptor at %u with score %d, terminalp %d/%d and introns:\n",
	 site->chrpos,site->score,site->known_terminalp,site->terminalp);
  for (p = site->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    printf("  to %u, count %d/%d, validp %d",
	   intron->prev_site->chrpos,Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),intron->validp);
    if (intron->prev_site == site->prev) {
      printf(" main");
    }
    for (q = site->alt_structs; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    for (q = site->alt_sites; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    for (q = site->alt_switches; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    printf("\n");
  }
  printf("\n");
  return;
}

static void
dump_acceptor_rev (Acceptor_T site) {
  Intron_T intron;
  Acceptor_alt_T alt;
  List_T p, q;

  printf("Rev acceptor at %u with score %d, terminalp %d/%d and introns:\n",
	 site->chrpos,site->score,site->known_terminalp,site->terminalp);
  for (p = site->introns; p != NULL; p = List_next(p)) {
    intron = (Intron_T) List_head(p);
    printf("  to %u, count %d/%d, validp %d",
	   intron->prev_site->chrpos,Splice_primary_count(intron->splice),Splice_crosshyb_count(intron->splice),intron->validp);
    if (intron->prev_site == site->prev) {
      printf(" main");
    }
    for (q = site->alt_structs; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    for (q = site->alt_sites; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    for (q = site->alt_switches; q != NULL; q = List_next(q)) {
      alt = (Acceptor_alt_T) List_head(q);
      if (intron->prev_site == alt->prev_site) {
	printf(" %s",Introntype_string(alt->type));
      }
    }
    printf("\n");
  }
  printf("\n");
  return;
}



static void
dump_graph_fwd (Genomicpos_T start, Genomicpos_T end,
		Donor_T *sites_donor_fwd, int nsites_donor_fwd,
		Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos;
  int donor_fwd_i = 0, acceptor_fwd_i = 0;

  printf("Entered dump_graph_fwd\n");
  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {

    if (donor_fwd_pos <= acceptor_fwd_pos) {
      if (donor_fwd_pos >= start && donor_fwd_pos <= end) {
	dump_donor_fwd(sites_donor_fwd[donor_fwd_i]);
      }
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else {
      if (acceptor_fwd_pos >= start && acceptor_fwd_pos <= end) {
	dump_acceptor_fwd(sites_acceptor_fwd[acceptor_fwd_i]);
      }
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  return;
}

static void
dump_graph_rev (Genomicpos_T start, Genomicpos_T end,
		Donor_T *sites_donor_rev, int nsites_donor_rev,
		Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos;
  int donor_rev_i = 0, acceptor_rev_i = 0;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {

    if (donor_rev_pos >= acceptor_rev_pos) {
      if (donor_rev_pos >= start && donor_rev_pos <= end) {
	dump_donor_rev(sites_donor_rev[donor_rev_i]);
      }
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else {
      if (acceptor_rev_pos >= start && acceptor_rev_pos <= end) {
	dump_acceptor_rev(sites_acceptor_rev[acceptor_rev_i]);
      }
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }
  }

  return;
}


/************************************************************************
 *   Construction of graphs
 ************************************************************************/

#if 0
/* Want to maximize number of exons, and not sure how to compare
   intron scores with exon scores */
static int
intron_score_iit (IIT_T splices_iit, Genomicpos_T low, Genomicpos_T high, char *chr, bool forwardp) {
  int sign;

  sign = (forwardp == true) ? +1 : -1;

  return (int) splice_match(splices_iit,chr,low,high,sign);
}
#endif


/************************************************************************
 *   Construction step 2. Adding exons
 ************************************************************************/

#if 0
static Genomicpos_T
edge_up (Genomicpos_T low, double *window_diff) {
  Genomicpos_T chrpos;
  double largest_diff = 0.0;
  Genomicpos_T bestpos;

  for (chrpos = low - 20; chrpos <= low + 20; chrpos++) {
    if (window_diff[chrpos] > largest_diff) {
      largest_diff = window_diff[chrpos];
      bestpos = chrpos;
    }
  }

  return bestpos;
}


static Genomicpos_T
edge_down (Genomicpos_T high, double *window_diff) {
  Genomicpos_T chrpos;
  long int largest_diff = 0.0;
  Genomicpos_T bestpos;

  for (chrpos = high + 20; chrpos >= high - 20; chrpos--) {
    if (-window_diff[chrpos] > largest_diff) {
      largest_diff = -window_diff[chrpos];
      bestpos = chrpos;
    }
  }

  return bestpos;
}



static bool
exon_evaluate (double *window_diff, Genomicpos_T low, Genomicpos_T high) {
  Genomicpos_T edge;

  edge = edge_up(low,window_diff);
  if (edge > low + 5 || edge < low - 5) {
    return false;
  }

  edge = edge_down(high,window_diff);
  if (edge > high + 5 || edge < high - 5) {
    return false;
  }
  
  return true;
}
#endif



static void
add_exons_fwd_aux (Donor_T site,
		   Donor_T *sites_donor_fwd, Acceptor_T *sites_acceptor_fwd,
#if 0
		   Site_T *sites_iblock_fwd,
#endif
		   int donor_fwd_i, int acceptor_fwd_i,
#if 0
		   , int iblock_fwd_i
#endif
		   int min_exonlength, int max_exonlength,
		   double *window_diff_low) {
  Acceptor_T prev_site;
  Genomicpos_T chrpos = site->chrpos, flat_bound, donor_bound;
  Score_T incrscore;
  int j;
  
  debug1(printf("Finding exons for position %u\n",chrpos));
  /* site->exons = (List_T) NULL; -- Don't remove information from previous run */

#ifdef COMPUTE_DONOR_BOUND
  /* Find prev donor that has exons */
  i = donor_fwd_i - 1;
  while (i >= 0 &&
	 (sites_donor_fwd[i]->chrpos + ALTSITE_DISTANCE >= chrpos ||
	  sites_donor_fwd[i]->exons == NULL)) {
    i--;
  }
  if (i < 0) {
    donor_bound = 0U;
  } else {
    donor_bound = sites_donor_fwd[i]->chrpos;
  }
#else
  donor_bound = 0U;
#endif

  flat_bound = flat_chrpos_low(window_diff_low,donor_bound,chrpos,/*allowed_zeroes*/20);
  debug1(printf("flat_bound from %u downto %u is at %u\n",chrpos,donor_bound,flat_bound));

#if 0
  if (iblock_fwd_i == 0) {
    iblock_bound = 0U;
  } else {
    iblock_bound = sites_iblock_fwd[iblock_fwd_i - 1]->chrpos;
  }
#endif

  /* This loop is necessary, even if min_exonlength is small */
  for (j = acceptor_fwd_i - 1; j >= 0 && sites_acceptor_fwd[j]->chrpos + min_exonlength > chrpos; j--) {
    debug1(printf("  Skipped acceptor at %u because of min_exonlength + %d > %u\n",
		  sites_acceptor_fwd[j]->chrpos,min_exonlength,chrpos));
  }

  /* If intergenic goes to intron positions, must be strict > and not >= */
  while (j >= 0 && 
	 sites_acceptor_fwd[j]->chrpos + max_exonlength >= chrpos && 
	 sites_acceptor_fwd[j]->chrpos > flat_bound) {
    prev_site = (Acceptor_T) sites_acceptor_fwd[j];
    prev_site->utr_p = false;

    /* Previously had exon_evaluate(window_diff_low,prev_site->chrpos,chrpos) == true) */
    if (Uinttable_get(site->exons_usedp,prev_site->chrpos) == NULL) {
      site->exons = List_push(site->exons,Exon_new(prev_site,site,/*incrscore*/1));
      Uinttable_put(site->exons_usedp,prev_site->chrpos,(void *) prev_site);
      debug1(printf("  Making exon from %u to %u\n",prev_site->chrpos,site->chrpos));
    }

    j--;
  }


  if (site->exons == NULL) {
    site->utr_p = true;
  } else {
    site->exons = List_reverse(site->exons);
  }

  return;
}


static void
add_exons_rev_aux (Donor_T site,
		   Donor_T *sites_donor_rev, Acceptor_T *sites_acceptor_rev,
#if 0
		   Site_T *sites_iblock_rev,
#endif
		   int donor_rev_i, int acceptor_rev_i,
#if 0
		   , int iblock_rev_i
#endif
		   int min_exonlength, int max_exonlength,
		   double *window_diff_high) {
  Acceptor_T prev_site;
  Genomicpos_T chrpos = site->chrpos, flat_bound, donor_bound;
  int j;
  
  debug1(printf("Finding exons for position %u\n",chrpos));
  /* site->exons = (List_T) NULL; -- Don't remove information from previous run */

#ifdef COMPUTE_DONOR_BOUND
  i = donor_rev_i - 1;
  while (i >= 0 && 
	 (sites_donor_rev[i]->chrpos <= chrpos + ALTSITE_DISTANCE ||
	  sites_donor_rev[i]->exons == NULL)) {
    i--;
  }
  if (i < 0) {
    donor_bound = /*chrlength*/ -1U;
  } else {
    donor_bound = sites_donor_rev[i]->chrpos;
  }
#else
  donor_bound = -1U;
#endif

  flat_bound = flat_chrpos_high(window_diff_high,chrpos,donor_bound,/*allowed_zeroes*/20);
  debug1(printf("flat_bound from %u upto %u is at %u\n",chrpos,donor_bound,flat_bound));

#if 0
  if (iblock_rev_i == 0) {
    iblock_bound = chrlength;
  } else {
    iblock_bound = sites_iblock_rev[iblock_rev_i - 1]->chrpos;
  }
#endif

  /* This loop is necessary, even if min_exonlength is small */
  for (j = acceptor_rev_i - 1; j >= 0 && sites_acceptor_rev[j]->chrpos < chrpos + min_exonlength; j--) {
    debug1(printf("  Skipped acceptor at %u because of min_exonlength - %u < %u\n",
		  sites_acceptor_rev[j]->chrpos,min_exonlength,chrpos));
  }

 /* If intergenic goes to intron positions, must be strict < and not <= */
  while (j >= 0 && 
	 sites_acceptor_rev[j]->chrpos <= chrpos + max_exonlength && 
	 sites_acceptor_rev[j]->chrpos < flat_bound) {
    
    prev_site = (Acceptor_T) sites_acceptor_rev[j];
    prev_site->utr_p = false;

    /* Previously had exon_evaluate(window_diff_high,chrpos,prev_site->chrpos) == true */
    if (Uinttable_get(site->exons_usedp,prev_site->chrpos) == NULL) {
      site->exons = List_push(site->exons,Exon_new(prev_site,site,/*incrscore*/1));
      Uinttable_put(site->exons_usedp,prev_site->chrpos,(void *) prev_site);
      debug1(printf("  Making exon from %u to %u\n",prev_site->chrpos,site->chrpos));
    }

    j--;
  }

  if (site->exons == NULL) {
    site->utr_p = true;
  } else {
    site->exons = List_reverse(site->exons);
  }

  return;
}


static void
add_exons_fwd (Donor_T *sites_donor_fwd, int nsites_donor_fwd,
	       Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd,
#if 0
	       , Site_T *sites_iblock_fwd, int nsites_iblock_fwd
#endif
	       int min_exonlength, int max_exonlength,
	       double *window_diff_low) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos, smallest;
  int donor_fwd_i = 0, acceptor_fwd_i = 0;

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  /* iblock_fwd_pos = (nsites_iblock_fwd == 0) ? (Genomicpos_T) -1U : sites_iblock_fwd[0]->chrpos; */

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd
#if 0
	 || iblock_fwd_i < nsites_iblock_fwd
#endif
	 ) {
    smallest = (donor_fwd_pos < acceptor_fwd_pos) ? donor_fwd_pos : acceptor_fwd_pos;
    /* smallest = (iblock_fwd_pos < smallest) ? iblock_fwd_pos : smallest; */

    if (donor_fwd_pos == smallest) {
      add_exons_fwd_aux(sites_donor_fwd[donor_fwd_i],sites_donor_fwd,sites_acceptor_fwd,
#if 0
			sites_iblock_fwd,
#endif
			donor_fwd_i,acceptor_fwd_i,
#if 0
			,iblock_fwd_i
#endif
			min_exonlength,max_exonlength,window_diff_low);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else if (acceptor_fwd_pos == smallest) {
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

#if 0
    } else {
      iblock_fwd_pos = (++iblock_fwd_i >= nsites_iblock_fwd) ? (Genomicpos_T) -1U : sites_iblock_fwd[iblock_fwd_i]->chrpos;
#endif
    }

  }

  return;
}


static void
add_exons_rev (Donor_T *sites_donor_rev, int nsites_donor_rev,
	       Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev,
#if 0
	       , Site_T *sites_iblock_rev, int nsites_iblock_rev
#endif
	       int min_exonlength, int max_exonlength,
	       double *window_diff_high) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos, largest;
  int donor_rev_i = 0, acceptor_rev_i = 0;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  /* iblock_rev_pos = (nsites_iblock_rev == 0) ? (Genomicpos_T) 0U : sites_iblock_rev[0]->chrpos; */

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev
#if 0
	 || iblock_rev_i < nsites_iblock_rev
#endif
	 ) {
    
    largest = (donor_rev_pos > acceptor_rev_pos) ? donor_rev_pos : acceptor_rev_pos;
    /* largest = (iblock_rev_pos > largest) ? iblock_rev_pos : largest; */

    if (donor_rev_pos == largest) {
      add_exons_rev_aux(sites_donor_rev[donor_rev_i],sites_donor_rev,sites_acceptor_rev,
			donor_rev_i,acceptor_rev_i,min_exonlength,max_exonlength,
			window_diff_high);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else if (acceptor_rev_pos == largest) {
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

#if 0
    } else {
      iblock_rev_pos = (++iblock_rev_i >= nsites_iblock_rev) ? (Genomicpos_T) 0U : sites_iblock_rev[iblock_rev_i]->chrpos;
#endif
    }

  }

  return;
}


static void
add_exons_fwd_known (Donor_T *sites_donor_fwd, int nsites_donor_fwd,
		     Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd,
		     IIT_T knownexons_iit, char *chr) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos, smallest, acceptorpos;
  int donor_fwd_i = 0, acceptor_fwd_i = 0;
  Donor_T site;
  Acceptor_T prev_site;
  Interval_T interval;
  int *matches, nmatches, i, j;

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;

  while (donor_fwd_i < nsites_donor_fwd || acceptor_fwd_i < nsites_acceptor_fwd) {
    smallest = (donor_fwd_pos < acceptor_fwd_pos) ? donor_fwd_pos : acceptor_fwd_pos;

    if (donor_fwd_pos == smallest) {
      site = sites_donor_fwd[donor_fwd_i];
      site->exons = (List_T) NULL; /* Okay to set to NULL for known sites */
      debug1(printf("Finding exons for position %u\n",site->chrpos));

      matches = IIT_get(&nmatches,knownexons_iit,chr,/*x*/smallest,/*y*/smallest,/*sortp*/false);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(knownexons_iit,matches[i]);
	if (Interval_sign(interval) > 0 && Interval_high(interval) == smallest) {
	  acceptorpos = Interval_low(interval);
	  for (j = acceptor_fwd_i - 1; j >= 0 && sites_acceptor_fwd[j]->chrpos != acceptorpos; j--) ;
	  if (j < 0) {
	    debug1(printf("  ??? No exon found for acceptorpos %u\n",acceptorpos));
	  } else {
	    prev_site = sites_acceptor_fwd[j];
	    site->exons = List_push(site->exons,Exon_new(prev_site,site,/*incrscore*/1));
	    Uinttable_put(site->exons_usedp,prev_site->chrpos,(void *) prev_site);
	    debug1(printf("  Making exon from %u to %u\n",prev_site->chrpos,site->chrpos));
	  }
	}
      }
      FREE(matches);
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else if (acceptor_fwd_pos == smallest) {
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

    }
  }

  return;
}


static void
add_exons_rev_known (Donor_T *sites_donor_rev, int nsites_donor_rev,
		     Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev,
		     IIT_T knownexons_iit, char *chr) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos, largest, acceptorpos;
  int donor_rev_i = 0, acceptor_rev_i = 0;
  Donor_T site;
  Acceptor_T prev_site;
  Interval_T interval;
  int *matches, nmatches, i, j;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;

  while (donor_rev_i < nsites_donor_rev || acceptor_rev_i < nsites_acceptor_rev) {
    largest = (donor_rev_pos > acceptor_rev_pos) ? donor_rev_pos : acceptor_rev_pos;

    if (donor_rev_pos == largest) {
      site = sites_donor_rev[donor_rev_i];
      site->exons = (List_T) NULL; /* Okay to set to NULL for known sites */
      debug1(printf("Finding exons for position %u\n",site->chrpos));

      matches = IIT_get(&nmatches,knownexons_iit,chr,/*x*/largest,/*y*/largest,/*sortp*/false);
      for (i = 0; i < nmatches; i++) {
	interval = IIT_interval(knownexons_iit,matches[i]);
	if (Interval_sign(interval) < 0 && Interval_low(interval) == largest) {
	  acceptorpos = Interval_high(interval);
	  for (j = acceptor_rev_i - 1; j >= 0 && sites_acceptor_rev[j]->chrpos != acceptorpos; j--) ;
	  if (j < 0) {
	    debug1(printf("  ??? No exon found for acceptorpos %u\n",acceptorpos));
	  } else {
	    prev_site = sites_acceptor_rev[j];
	    site->exons = List_push(site->exons,Exon_new(prev_site,site,/*incrscore*/1));
	    Uinttable_put(site->exons_usedp,prev_site->chrpos,(void *) prev_site);
	    debug1(printf("  Making exon from %u to %u\n",prev_site->chrpos,site->chrpos));
	  }
	}
      }
      FREE(matches);
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else if (acceptor_rev_pos == largest) {
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

    }

  }

  return;
}




/************************************************************************
 *   Construction step 5. Adding introns
 ************************************************************************/

static void
add_introns_fwd_aux (Acceptor_T site, Donor_T *sites_donor_fwd, int donor_fwd_i,
#if 0
		     , Site_T *sites_iblock_fwd, , int iblock_fwd_i
#endif
		     IIT_T introns_iit, char *chr,
		     int min_intronlength, int max_intronlength, bool knownp) {
  Donor_T prev_site;
  Genomicpos_T chrpos = site->chrpos, iblock_bound;
  Splice_T splice;
  int count;
  int j;

  debug2(printf("Finding introns for position %u\n",chrpos));
  /* site->introns = (List_T) NULL; -- Don't remove information from previous run */

#if 0
  if (iblock_fwd_i == 0) {
    iblock_bound = 0U;
  } else {
    iblock_bound = sites_iblock_fwd[iblock_fwd_i - 1]->chrpos;
  }
#else
  if (chrpos < (Genomicpos_T) max_intronlength) {
    iblock_bound = 0U;
  } else {
    iblock_bound = chrpos - max_intronlength;
  }
#endif
  debug2(printf("iblock_bound from %u at %u\n",chrpos,iblock_bound));

  for (j = donor_fwd_i - 1; j >= 0 && sites_donor_fwd[j]->chrpos + min_intronlength >= chrpos; j--) {
    debug2(printf("  Skipped donor at %u because of min_intronlength + %u >= %u\n",
		  sites_donor_fwd[j]->chrpos,min_intronlength,chrpos));
  }

 /* If intergenic block goes to intron positions, must be >= and not > */
  while (j >= 0 && sites_donor_fwd[j]->chrpos >= iblock_bound) {
    prev_site = (Donor_T) sites_donor_fwd[j];

    if ((splice = (Splice_T) Uinttable_get(site->splices_table,/*donorpos*/prev_site->chrpos)) != NULL &&
	Splice_intron_usedp(splice) == false) {
      site->introns = List_push(site->introns,Intron_new(prev_site,site,splice,introns_iit,chr,
							 /*forwardp*/true));

      if (knownp == true) {
	count = Splice_known_count(splice);
      } else {
	count = Splice_primary_count(splice);
	/* crosshyb_count = Splice_crosshyb_count(splice); */
      }

      site->points += count;
      prev_site->points += count;
      site->score += count; /* for start of chain */
      prev_site->score += count;
      Splice_set_intron_usedp(splice);
	
      debug2(printf("  Making intron from %u to %u, count %d\n",
		    prev_site->chrpos,site->chrpos,count));
    }

    j--;
  }

  site->introns = List_reverse(site->introns);

  return;
}


static void
add_introns_rev_aux (Acceptor_T site, Donor_T *sites_donor_rev, int donor_rev_i,
#if 0
 		     , Site_T *sites_iblock_rev, int iblock_rev_i
#endif
		     IIT_T introns_iit, char *chr,
		     int min_intronlength, int max_intronlength, bool knownp) {
  Donor_T prev_site;
  Genomicpos_T chrpos = site->chrpos, iblock_bound;
  Splice_T splice;
  int count;
  int j;

  debug2(printf("Finding introns for position %u\n",chrpos));
  /* site->introns = (List_T) NULL; -- Don't remove information from previous run */

#if 0
  if (iblock_rev_i == 0) {
    iblock_bound = chrlength;
  } else {
    iblock_bound = sites_iblock_rev[iblock_rev_i - 1]->chrpos;
  }
#else
  iblock_bound = chrpos + max_intronlength;
#endif
  debug2(printf("iblock_bound from %u at %u\n",chrpos,iblock_bound));

  for (j = donor_rev_i - 1; j >= 0 && sites_donor_rev[j]->chrpos <= chrpos + min_intronlength; j--) {
    debug2(printf("  Skipped donor at %u because of min_intronlength - %u <= %u\n",
		  sites_donor_rev[j]->chrpos,min_intronlength,chrpos));
  }

 /* If intergenic block goes to intron positions, must be >= and not > */
  while (j >= 0 && sites_donor_rev[j]->chrpos <= iblock_bound) {
    prev_site = (Donor_T) sites_donor_rev[j];

    if ((splice = (Splice_T) Uinttable_get(site->splices_table,/*donorpos*/prev_site->chrpos)) != NULL &&
	Splice_intron_usedp(splice) == false) {
      site->introns = List_push(site->introns,Intron_new(prev_site,site,splice,introns_iit,chr,
							 /*forwardp*/false));

      if (knownp == true) {
	count = Splice_known_count(splice);
      } else {
	count = Splice_primary_count(splice);
      }

      site->points += count;
      prev_site->points += count;
      site->score += count; /* for start of chain */
      prev_site->score += count;
      Splice_set_intron_usedp(splice);

      debug2(printf("  Making intron from %u to %u, count %d\n",
		    prev_site->chrpos,site->chrpos,count));
    }

    j--;
  }

  site->introns = List_reverse(site->introns);

  return;
}




static void
add_introns_fwd (Donor_T *sites_donor_fwd, int nsites_donor_fwd,
		 Acceptor_T *sites_acceptor_fwd, int nsites_acceptor_fwd,
#if 0
		 , Site_T *sites_iblock_fwd, int nsites_iblock_fwd
#endif
		 IIT_T introns_iit, char *chr,
		 int min_intronlength, int max_intronlength, bool knownp) {
  Genomicpos_T donor_fwd_pos, acceptor_fwd_pos, smallest;
  int donor_fwd_i = 0, acceptor_fwd_i = 0;

  donor_fwd_pos = (nsites_donor_fwd == 0) ? (Genomicpos_T) -1U : sites_donor_fwd[0]->chrpos;
  acceptor_fwd_pos = (nsites_acceptor_fwd == 0) ? (Genomicpos_T) -1U : sites_acceptor_fwd[0]->chrpos;
  /* iblock_fwd_pos = (nsites_iblock_fwd == 0) ? (Genomicpos_T) -1U : sites_iblock_fwd[0]->chrpos; */

  while (acceptor_fwd_i < nsites_acceptor_fwd) {
    smallest = (donor_fwd_pos < acceptor_fwd_pos) ? donor_fwd_pos : acceptor_fwd_pos;
    /* smallest = (iblock_fwd_pos < smallest) ? iblock_fwd_pos : smallest; */

    if (donor_fwd_pos == smallest) {
      donor_fwd_pos = (++donor_fwd_i >= nsites_donor_fwd) ? (Genomicpos_T) -1U : sites_donor_fwd[donor_fwd_i]->chrpos;

    } else if (acceptor_fwd_pos == smallest) {
      add_introns_fwd_aux(sites_acceptor_fwd[acceptor_fwd_i],
			  sites_donor_fwd,donor_fwd_i,introns_iit,chr,
			  min_intronlength,max_intronlength,knownp);
      acceptor_fwd_pos = (++acceptor_fwd_i >= nsites_acceptor_fwd) ? (Genomicpos_T) -1U : sites_acceptor_fwd[acceptor_fwd_i]->chrpos;

#if 0
    } else {
      iblock_fwd_pos = (++iblock_fwd_i >= nsites_iblock_fwd) ? (Genomicpos_T) -1U : sites_iblock_fwd[iblock_fwd_i]->chrpos;
#endif

    }
  }

  return;
}


static void
add_introns_rev (Donor_T *sites_donor_rev, int nsites_donor_rev,
		 Acceptor_T *sites_acceptor_rev, int nsites_acceptor_rev,
#if 0
		 , Site_T *sites_iblock_rev, int nsites_iblock_rev
#endif
		 IIT_T introns_iit, char *chr,
		 int min_intronlength, int max_intronlength, bool knownp) {
  Genomicpos_T donor_rev_pos, acceptor_rev_pos, largest;
  int donor_rev_i = 0, acceptor_rev_i = 0;

  donor_rev_pos = (nsites_donor_rev == 0) ? (Genomicpos_T) 0U : sites_donor_rev[0]->chrpos;
  acceptor_rev_pos = (nsites_acceptor_rev == 0) ? (Genomicpos_T) 0U : sites_acceptor_rev[0]->chrpos;
  /* iblock_rev_pos = (nsites_iblock_rev == 0) ? (Genomicpos_T) 0U : sites_iblock_rev[0]->chrpos; */

  while (acceptor_rev_i < nsites_acceptor_rev) {
    largest = (donor_rev_pos > acceptor_rev_pos) ? donor_rev_pos : acceptor_rev_pos;
    /* largest = (iblock_rev_pos > largest) ? iblock_rev_pos : largest; */

    if (donor_rev_pos == largest) {
      donor_rev_pos = (++donor_rev_i >= nsites_donor_rev) ? (Genomicpos_T) 0U : sites_donor_rev[donor_rev_i]->chrpos;

    } else if (acceptor_rev_pos == largest) {
      add_introns_rev_aux(sites_acceptor_rev[acceptor_rev_i],
			  sites_donor_rev,donor_rev_i,introns_iit,chr,
			  min_intronlength,max_intronlength,knownp);
      acceptor_rev_pos = (++acceptor_rev_i >= nsites_acceptor_rev) ? (Genomicpos_T) 0U : sites_acceptor_rev[acceptor_rev_i]->chrpos;

#if 0
    } else {
      iblock_rev_pos = (++iblock_rev_i >= nsites_iblock_rev) ? (Genomicpos_T) 0U : sites_iblock_rev[iblock_rev_i]->chrpos;
#endif

    }
  }

  return;
}




/************************************************************************/

#if 0
static void
post_sites (Site_T **sites, int *nsites, IIT_T genebounds_iit,
	    char *chr, Genomicpos_T chrstart, Genomicpos_T chrend, int sign) {
  List_T p;
  int *matches, nmatches, i, k = 0;
  Interval_T interval;
  int divno;

  divno = IIT_divint(genebounds_iit,chr);
  matches = IIT_get_signed_with_divno(&nmatches,genebounds_iit,divno,chrstart,chrend,/*sortp*/false,sign);

  if (nmatches == 0) {
    *nsites = 0;
    *sites = (Site_T *) NULL;
  } else {
    *nsites = 2*nmatches;
    *sites = (Site_T *) CALLOC(*nsites,sizeof(Site_T));
    for (i = 0; i < nmatches; i++) {
      interval = IIT_interval(genebounds_iit,matches[i]);
      (*sites)[k++] = Site_new(Interval_low(interval));
      (*sites)[k++] = Site_new(Interval_high(interval));
    }
    FREE(matches);

    if (sign > 0) {
      qsort(*sites,*nsites,sizeof(Site_T),Site_ascending_cmp);
    } else if (sign < 0) {
      qsort(*sites,*nsites,sizeof(Site_T),Site_descending_cmp);
    }
  }

  return;
}
#endif


static Donor_T *
collect_donors (int *nsites, Uinttable_T sitestable, int sign) {
  Donor_T *sites;
  Genomicpos_T *keys;
  int i, k;

  if ((*nsites = Uinttable_length(sitestable)) == 0) {
    return (Donor_T *) NULL;

  } else if (sign > 0) {
    sites = (Donor_T *) CALLOC(*nsites,sizeof(Donor_T));

    keys = (Genomicpos_T *) Uinttable_keys(sitestable,/*sortp*/true);
    k = 0;
    for (i = 0; i < *nsites; i++) {
      sites[k++] = (Donor_T) Uinttable_get(sitestable,keys[i]);
    }
    FREE(keys);
    return sites;

  } else {
    sites = (Donor_T *) CALLOC(*nsites,sizeof(Donor_T));

    keys = (Genomicpos_T *) Uinttable_keys(sitestable,/*sortp*/true);
    k = 0;
    for (i = (*nsites) - 1; i >= 0; i--) {
      sites[k++] = (Donor_T) Uinttable_get(sitestable,keys[i]);
    }
    FREE(keys);
    return sites;
  }
}


static Acceptor_T *
collect_acceptors (int *nsites, Uinttable_T sitestable, int sign) {
  Acceptor_T *sites, acceptor;
  Genomicpos_T *keys;
  int i, k;

  if ((*nsites = Uinttable_length(sitestable)) == 0) {
    return (Acceptor_T *) NULL;

  } else if (sign > 0) {
    sites = (Acceptor_T *) CALLOC(*nsites,sizeof(Acceptor_T));

    keys = (Genomicpos_T *) Uinttable_keys(sitestable,/*sortp*/true);
    k = 0;
    for (i = 0; i < *nsites; i++) {
      sites[k++] = acceptor = (Acceptor_T) Uinttable_get(sitestable,keys[i]);
      acceptor->paths_usedp = false;
    }
    FREE(keys);
    return sites;

  } else {
    sites = (Acceptor_T *) CALLOC(*nsites,sizeof(Acceptor_T));

    keys = (Genomicpos_T *) Uinttable_keys(sitestable,/*sortp*/true);
    k = 0;
    for (i = (*nsites) - 1; i >= 0; i--) {
      sites[k++] = acceptor = (Acceptor_T) Uinttable_get(sitestable,keys[i]);
      acceptor->paths_usedp = false;
    }
    FREE(keys);
    return sites;
  }
}



/* Written to minimize memory usage */
static void
build_graph_fwd_obs (T this, char *chr, List_T obs_splices, List_T knowngenes,
		     bool forwardp, int min_exonlength, int max_exonlength,
		     int min_intronlength, int max_intronlength,
		     double *window_diff_low) {

  /* Step 1. Put and collect splice sites */
  put_splice_sites(this->donor_fwd_sitestable,this->acceptor_fwd_sitestable,
		   obs_splices,knowngenes,forwardp);
  fprintf(stderr,"Chr %s has %d donor and %d donor and acceptor %s sites from splices\n",
	  chr,Uinttable_length(this->donor_fwd_sitestable),Uinttable_length(this->acceptor_fwd_sitestable),
	  (forwardp == true) ? "forward" : "reverse");

  this->sites_donor_fwd =
    collect_donors(&this->nsites_donor_fwd,this->donor_fwd_sitestable,/*sign*/+1);
  this->sites_acceptor_fwd =
    collect_acceptors(&this->nsites_acceptor_fwd,this->acceptor_fwd_sitestable,/*sign*/+1);


#if 0
  /* Step 3. Post intergenic blocks */
  post_sites(&sites_iblock_fwd,&nsites_iblock_fwd,genebounds_iit,chr,chrstart,chrend,/*sign*/+1);
  fprintf(stderr,"Found %d %s intergenic blocks\n",
	  nsites_iblock_fwd,(forwardp == true) ? "forward" : "reverse");
#endif

  /* Step 2.  Add exons and introns, which requires prior donor/acceptor/iblocks */
  /* Removed parameters for sites_iblock_fwd,nsites_iblock_fwd */
  add_exons_fwd(this->sites_donor_fwd,this->nsites_donor_fwd,
		this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
		min_exonlength,max_exonlength,window_diff_low);
  add_introns_fwd(this->sites_donor_fwd,this->nsites_donor_fwd,
		  this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
		  /*introns_iit*/NULL,/*chr*/NULL,
		  min_intronlength,max_intronlength,/*knownp*/false);

  return;
}


static void
build_graph_rev_obs (T this, char *chr, List_T obs_splices, List_T knowngenes,
		     bool forwardp, int min_exonlength, int max_exonlength,
		     int min_intronlength, int max_intronlength,
		     double *window_diff_high) {

  /* Step 1. Put and collect splice sites */
  put_splice_sites(this->donor_rev_sitestable,this->acceptor_rev_sitestable,
		   obs_splices,knowngenes,forwardp);
  fprintf(stderr,"Chr %s has %d donor and %d donor and acceptor %s sites from splices\n",
	  chr,Uinttable_length(this->donor_rev_sitestable),Uinttable_length(this->acceptor_rev_sitestable),
	  (forwardp == true) ? "forward" : "reverse");

  this->sites_donor_rev =
    collect_donors(&this->nsites_donor_rev,this->donor_rev_sitestable,/*sign*/-1);
  this->sites_acceptor_rev =
    collect_acceptors(&this->nsites_acceptor_rev,this->acceptor_rev_sitestable,/*sign*/-1);


#if 0
  /* Step 3. Post intergenic blocks */
  post_sites(&sites_iblock_rev,&nsites_iblock_rev,genebounds_iit,chr,chrstart,chrend,/*sign*/-1);
  fprintf(stderr,"Found %d %s intergenic blocks\n",
	  nsites_iblock_rev,(forwardp == true) ? "forward" : "reverse");
#endif


  /* Step 2.  Add exons and introns, which requires prior donor/acceptor/iblocks */
  /* Removed parameters for sites_iblock_rev,nsites_iblock_rev */
  add_exons_rev(this->sites_donor_rev,this->nsites_donor_rev,
		this->sites_acceptor_rev,this->nsites_acceptor_rev,
		min_exonlength,max_exonlength,window_diff_high);
  add_introns_rev(this->sites_donor_rev,this->nsites_donor_rev,
		  this->sites_acceptor_rev,this->nsites_acceptor_rev,
		  /*introns_iit*/NULL,/*chr*/NULL,
		  min_intronlength,max_intronlength,/*knownp*/false);

  return;
}



static void
build_graph_fwd_known (T this, char *chr, List_T splices, IIT_T introns_iit,
		       IIT_T knownexons_iit, bool forwardp) {

  /* Step 1. Put and collect splice sites */
  put_splice_sites(this->donor_fwd_sitestable,this->acceptor_fwd_sitestable,
		   /*obs_splices*/splices,/*knowngenes*/NULL,forwardp);
  fprintf(stderr,"Chr %s has %d donor and %d donor and acceptor %s sites from splices\n",
	  chr,Uinttable_length(this->donor_fwd_sitestable),Uinttable_length(this->acceptor_fwd_sitestable),
	  (forwardp == true) ? "forward" : "reverse");

  this->sites_donor_fwd =
    collect_donors(&this->nsites_donor_fwd,this->donor_fwd_sitestable,/*sign*/+1);
  this->sites_acceptor_fwd =
    collect_acceptors(&this->nsites_acceptor_fwd,this->acceptor_fwd_sitestable,/*sign*/+1);

  /* Step 2.  Add exons and introns, which requires prior donor/acceptor/iblocks */
  /* Removed parameters for sites_iblock_fwd,nsites_iblock_fwd */
  add_exons_fwd_known(this->sites_donor_fwd,this->nsites_donor_fwd,
		      this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
		      knownexons_iit,chr);
  add_introns_fwd(this->sites_donor_fwd,this->nsites_donor_fwd,
		  this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
		  introns_iit,chr,/*min_intronlength*/0,
		  /*max_intronlength*/1000000,/*knownp*/true);

  return;
}


static void
build_graph_rev_known (T this, char *chr, List_T splices, IIT_T introns_iit,
		       IIT_T knownexons_iit, bool forwardp) {

  /* Step 1. Put and collect splice sites */
  put_splice_sites(this->donor_rev_sitestable,this->acceptor_rev_sitestable,
		   /*obs_splices*/splices,/*knowngenes*/NULL,forwardp);
  fprintf(stderr,"Chr %s has %d donor and %d donor and acceptor %s sites from splices\n",
	  chr,Uinttable_length(this->donor_rev_sitestable),Uinttable_length(this->acceptor_rev_sitestable),
	  (forwardp == true) ? "forward" : "reverse");

  this->sites_donor_rev =
    collect_donors(&this->nsites_donor_rev,this->donor_rev_sitestable,/*sign*/-1);
  this->sites_acceptor_rev =
    collect_acceptors(&this->nsites_acceptor_rev,this->acceptor_rev_sitestable,/*sign*/-1);

  /* Step 2.  Add exons and introns, which requires prior donor/acceptor/iblocks */
  /* Removed parameters for sites_iblock_rev,nsites_iblock_rev */
  add_exons_rev_known(this->sites_donor_rev,this->nsites_donor_rev,
		      this->sites_acceptor_rev,this->nsites_acceptor_rev,
		      knownexons_iit,chr);
  add_introns_rev(this->sites_donor_rev,this->nsites_donor_rev,
		  this->sites_acceptor_rev,this->nsites_acceptor_rev,
		  introns_iit,chr,/*min_intronlength*/0,
		  /*max_intronlength*/1000000,/*knownp*/true);

  return;
}



List_T
Splicegraph_solve_obs (T this, List_T obs_splices, List_T knowngenes, char *chr,
		       long int *tally_matches_low, long int *tally_mismatches_low,
		       long int *tally_matches_high, long int *tally_mismatches_high,
		       long int *primary_extents, long int *crosshyb_extents,
		       Genome_T genome, Genomicpos_T chroffset, Genomicpos_T chrlength,
		       int insertlength, int readlength, int min_overhang,
		       int mincount_alt, int minsupport, bool need_canonical_p,
		       int min_exonlength, int max_exonlength, int auto_exonlength,
		       int min_intronlength, int max_intronlength, bool altpaths_p) {
  List_T genes;

  double *window_diff_low, *window_diff_high, *cumlog_tally, *log_tally;

  int nconflicts;
  List_T terminals_fwd = NULL, terminals_rev = NULL, p;
  List_T terminals;
  Acceptor_T terminal;
  Score_T pathscore;

  Genomicpos_T *keys;
  Uinttable_T terminal_table_fwd, terminal_table_rev;
  int n, i;

  fprintf(stderr,"Solving for %s, length %u\n",chr,chrlength);

  /* Low */
  log_tally = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_log_tally(log_tally,tally_matches_low,tally_mismatches_low,chrlength);
  
  cumlog_tally = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_cum_double(cumlog_tally,log_tally,chrlength);
  FREE(log_tally);

  window_diff_low = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_diff(window_diff_low,cumlog_tally,chrlength);
  FREE(cumlog_tally);

  /* High */
  log_tally = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_log_tally(log_tally,tally_matches_high,tally_mismatches_high,chrlength);

  cumlog_tally = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_cum_double(cumlog_tally,log_tally,chrlength);
  FREE(log_tally);

  window_diff_high = (double *) CALLOC(chrlength+1,sizeof(double));
  compute_diff(window_diff_high,cumlog_tally,chrlength);
  FREE(cumlog_tally);


  build_graph_fwd_obs(this,chr,obs_splices,knowngenes,/*forwardp*/true,
		      min_exonlength,max_exonlength,min_intronlength,max_intronlength,
		      window_diff_low);
  build_graph_rev_obs(this,chr,obs_splices,knowngenes,/*forwardp*/false,
		      min_exonlength,max_exonlength,min_intronlength,max_intronlength,
		      window_diff_high);

#if 0
  /* Used to cycle back to here based on nconflicts */
  reset_donor_sites(sites_donor_fwd,nsites_donor_fwd);
  reset_acceptor_sites(sites_acceptor_fwd,nsites_acceptor_fwd);
  reset_donor_sites(sites_donor_rev,nsites_donor_rev);
  reset_acceptor_sites(sites_acceptor_rev,nsites_acceptor_rev);
  List_free(&terminals_fwd);
  List_free(&terminals_rev);
#endif

  fprintf(stderr,"Traversing forward graph...");
  traverse_graph_fwd(this->sites_donor_fwd,this->nsites_donor_fwd,
		     this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
		     tally_matches_low,tally_mismatches_low,
		     genome,chr,chroffset,chrlength);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Traversing reverse graph...");
  traverse_graph_rev(this->sites_donor_rev,this->nsites_donor_rev,
		     this->sites_acceptor_rev,this->nsites_acceptor_rev,
		     tally_matches_high,tally_mismatches_high,
		     genome,chr,chroffset,chrlength);
  fprintf(stderr,"done\n");
  
  terminals_fwd = Terminals_gather_fwd_list(this->sites_acceptor_fwd,this->nsites_acceptor_fwd);
  terminals_rev = Terminals_gather_rev_list(this->sites_acceptor_rev,this->nsites_acceptor_rev);

  for (p = terminals_fwd; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    /* Might be better to recompute */
    pathscore = terminal->score;
    assign_pathscore_to_introns(terminal,pathscore);
  }

  for (p = terminals_rev; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    /* Might be better to recompute */
    pathscore = terminal->score;
    assign_pathscore_to_introns(terminal,pathscore);
  }

#if 0
  /* Causes misses of reasonable genes */
  nconflicts = find_conflicting_introns(terminals_fwd,terminals_rev,chr);
  fprintf(stderr,"%d conflicts found\n",nconflicts);
#else
  nconflicts = 0;
#endif
  /* Used to cycle back based on nconflicts */


  debug(printf("*** Fwd graph ***\n"));
  debug(dump_graph_fwd(/*start*/0,/*end*/-1U,
		       this->sites_donor_fwd,this->nsites_donor_fwd,
		       this->sites_acceptor_fwd,this->nsites_acceptor_fwd));
  debug(printf("*** Rev graph ***\n"));
  debug(dump_graph_rev(/*start*/0,/*end*/-1U,
		       this->sites_donor_rev,this->nsites_donor_rev,
		       this->sites_acceptor_rev,this->nsites_acceptor_rev));

#if 0
  Terminals_set_pathbounds(terminals_fwd);
  Terminals_set_pathbounds(terminals_rev);
#endif
  List_free(&terminals_rev);
  List_free(&terminals_fwd);

  terminal_table_fwd = Terminals_gather_table(this->sites_donor_fwd,this->nsites_donor_fwd,
					      this->sites_acceptor_fwd,this->nsites_acceptor_fwd);
  terminal_table_rev = Terminals_gather_table(this->sites_donor_rev,this->nsites_donor_rev,
					      this->sites_acceptor_rev,this->nsites_acceptor_rev);

  genes = (List_T) NULL;
  if ((n = Uinttable_length(terminal_table_fwd)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(terminal_table_fwd,/*sortp*/true);
    for (i = 0; i < n; i++) {
      terminals = List_reverse(Uinttable_get(terminal_table_fwd,keys[i]));
      genes = Terminals_find_genes_obs(genes,terminals,chr,/*genei*/i+1,tally_matches_low,tally_mismatches_low,
				       tally_matches_high,tally_mismatches_high,
				       window_diff_low,window_diff_high,
				       primary_extents,crosshyb_extents,/*end_exons_iit*/NULL,
				       insertlength,readlength,min_overhang,
				       chrlength,mincount_alt,altpaths_p,/*forwardp*/true);
      List_free(&terminals);
    }
    FREE(keys);
  }

  if ((n = Uinttable_length(terminal_table_rev)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(terminal_table_rev,/*sortp*/true);
    for (i = 0; i < n; i++) {
      terminals = List_reverse(Uinttable_get(terminal_table_rev,keys[i]));
      genes = Terminals_find_genes_obs(genes,terminals,chr,/*genei*/i+1,tally_matches_low,tally_mismatches_low,
				       tally_matches_high,tally_mismatches_high,
				       window_diff_low,window_diff_high,
				       primary_extents,crosshyb_extents,/*end_exons_iit*/NULL,
				       insertlength,readlength,min_overhang,
				       chrlength,mincount_alt,altpaths_p,/*forwardp*/false);
      List_free(&terminals);
    }
    FREE(keys);
  }

  FREE(window_diff_high);
  FREE(window_diff_low);

  Uinttable_free(&terminal_table_rev);
  Uinttable_free(&terminal_table_fwd);

  FREE(this->sites_donor_fwd);
  FREE(this->sites_acceptor_fwd);
  FREE(this->sites_donor_rev);
  FREE(this->sites_acceptor_rev);

  /* Reverse list so they go in chromosomal order */
  return List_reverse(genes);
}



List_T
Splicegraph_solve_known (T this, List_T splices, IIT_T introns_iit,
			 IIT_T middle_exons_iit, IIT_T end_exons_iit,
			 char *chr, Genomicpos_T chrlength) {
  List_T genes;

  List_T terminals_fwd = NULL, terminals_rev = NULL, p;
  List_T terminals;
  Acceptor_T terminal;
  Score_T pathscore;

  Genomicpos_T *keys;
  Uinttable_T terminal_table_fwd, terminal_table_rev;
  int n, i;


  build_graph_fwd_known(this,chr,splices,introns_iit,middle_exons_iit,/*forwardp*/true);
  build_graph_rev_known(this,chr,splices,introns_iit,middle_exons_iit,/*forwardp*/false);

#if 0
  reset_donor_sites(sites_donor_fwd,nsites_donor_fwd);
  reset_acceptor_sites(sites_acceptor_fwd,nsites_acceptor_fwd);
  reset_donor_sites(sites_donor_rev,nsites_donor_rev);
  reset_acceptor_sites(sites_acceptor_rev,nsites_acceptor_rev);
  List_free(&terminals_fwd);
  List_free(&terminals_rev);
#endif


  fprintf(stderr,"Traversing forward graph of known transcripts...");
  traverse_graph_fwd_known(this->sites_donor_fwd,this->nsites_donor_fwd,
			   this->sites_acceptor_fwd,this->nsites_acceptor_fwd,
			   end_exons_iit,chr);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Traversing reverse graph of known transcripts...");
  traverse_graph_rev_known(this->sites_donor_rev,this->nsites_donor_rev,
			   this->sites_acceptor_rev,this->nsites_acceptor_rev,
			   end_exons_iit,chr);
  fprintf(stderr,"done\n");
  
  terminals_fwd = Terminals_gather_fwd_list(this->sites_acceptor_fwd,this->nsites_acceptor_fwd);
  terminals_rev = Terminals_gather_rev_list(this->sites_acceptor_rev,this->nsites_acceptor_rev);

  for (p = terminals_fwd; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    /* Might be better to recompute */
    pathscore = terminal->score;
    assign_pathscore_to_introns(terminal,pathscore);
  }

  for (p = terminals_rev; p != NULL; p = List_next(p)) {
    terminal = (Acceptor_T) List_head(p);
    /* Might be better to recompute */
    pathscore = terminal->score;
    assign_pathscore_to_introns(terminal,pathscore);
  }


  debug(printf("*** Fwd graph ***\n"));
  debug(dump_graph_fwd(/*start*/0,/*end*/-1U,
		       this->sites_donor_fwd,this->nsites_donor_fwd,
		       this->sites_acceptor_fwd,this->nsites_acceptor_fwd));
  debug(printf("*** Rev graph ***\n"));
  debug(dump_graph_rev(/*start*/0,/*end*/-1U,
		       this->sites_donor_rev,this->nsites_donor_rev,
		       this->sites_acceptor_rev,this->nsites_acceptor_rev));

  List_free(&terminals_rev);
  List_free(&terminals_fwd);

  terminal_table_fwd = Terminals_gather_table(this->sites_donor_fwd,this->nsites_donor_fwd,
					      this->sites_acceptor_fwd,this->nsites_acceptor_fwd);
  terminal_table_rev = Terminals_gather_table(this->sites_donor_rev,this->nsites_donor_rev,
					      this->sites_acceptor_rev,this->nsites_acceptor_rev);

  genes = (List_T) NULL;
  if ((n = Uinttable_length(terminal_table_fwd)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(terminal_table_fwd,/*sortp*/true);
    for (i = 0; i < n; i++) {
      terminals = List_reverse(Uinttable_get(terminal_table_fwd,keys[i]));
      genes = Terminals_find_genes_known(genes,terminals,chr,/*genei*/i+1,
					 /*tally_matches_low*/NULL,/*tally_mismatches_low*/NULL,
					 /*tally_matches_high*/NULL,/*tally_mismatches_high*/NULL,
					 /*primary_extents*/NULL,/*crosshyb_extents*/NULL,
					 middle_exons_iit,end_exons_iit,chrlength,
					 /*mincount_alt*/0,/*altpaths_p*/true,/*forwardp*/true);
      List_free(&terminals);
    }
    FREE(keys);
  }

  if ((n = Uinttable_length(terminal_table_rev)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(terminal_table_rev,/*sortp*/true);
    for (i = 0; i < n; i++) {
      terminals = List_reverse(Uinttable_get(terminal_table_rev,keys[i]));
      genes = Terminals_find_genes_known(genes,terminals,chr,/*genei*/i+1,
					 /*tally_matches_low*/NULL,/*tally_mismatches_low*/NULL,
					 /*tally_matches_high*/NULL,/*tally_mismatches_high*/NULL,
					 /*primary_extents*/NULL,/*crosshyb_extents*/NULL,
					 middle_exons_iit,end_exons_iit,chrlength,
					 /*mincount_alt*/0,/*altpaths_p*/true,/*forwardp*/false);
      List_free(&terminals);
    }
    FREE(keys);
  }

  Uinttable_free(&terminal_table_rev);
  Uinttable_free(&terminal_table_fwd);

  FREE(this->sites_donor_fwd);
  FREE(this->sites_acceptor_fwd);
  FREE(this->sites_donor_rev);
  FREE(this->sites_acceptor_rev);

  return genes;
}


#if 0
int
main (int argc, char *argv[]) {
  IIT_T splices_iit;

  char *iitfile;
  char *chr;
  int divno;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"A:Me:E:i:V?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'A':
      if (!strcmp(optarg,"paths")) {
	output_type = OUTPUT_PATHS;
      } else if (!strcmp(optarg,"genes")) {
	output_type = OUTPUT_GENES;
      } else if (!strcmp(optarg,"diffs")) {
	output_type = OUTPUT_DIFFS;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
      }
      break;

    case 'M': altpaths_p = false; break;

    case 'e': min_exonlength = atoi(optarg); break;
    case 'E': max_exonlength = atoi(optarg); break;
    case 'i': min_intronlength = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);
      
  if (argc < 2) {
    /* argv[0] is program name */
    printf("Usage: splicegene <splices iit>\n");
    exit(9);
  }

  if ((splices_iit = IIT_read(argv[1],/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
    fprintf(stderr,"Could not open IIT file %s\n",argv[1]);
    exit(9);
  }

  for (divno = 1; divno < IIT_ndivs(splices_iit); divno++) {
    chr = IIT_divstring(splices_iit,divno);
    solve_chromosome(chr,splices_iit);
  }
  
  IIT_free(&splices_iit);

  return 0;
}
#endif

