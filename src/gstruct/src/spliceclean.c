static char rcsid[] = "$Id: spliceclean.c 47265 2011-09-14 20:50:59Z twu $";
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
#include <ctype.h>
#include <math.h>		/* For qsort */


#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "list.h"
#include "chrom.h"
#include "table.h"
#include "uinttable.h"

#include "iit-read.h"
#include "tally.h"

#include "getopt.h"



#define MAXLOOKBACK 30
#define MAX_WARNINGS 20

/* parse_input */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* compute_scores */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* traceback_ceiling */
#ifdef DEBUG3C
#define debug3c(x) x
#else
#define debug3c(x)
#endif

/* traceback_floor */
#ifdef DEBUG3F
#define debug3f(x) x
#else
#define debug3f(x)
#endif


/************************************************************************
 *   Program options
 ************************************************************************/

static IIT_T extents_fwd_iit = NULL;
static IIT_T extents_rev_iit = NULL;
static bool resolve_directions_p = false;

static IIT_T genebounds_iit = NULL;

static bool dump_p = false;
static bool noclean_p = false;
static bool print_included_p = true;	/* If false, prints excluded splices */
static bool runlength_p = false;
static bool print_sites_p = false;
static bool endpoints_p = false;
static bool midpoints_p = false;

static bool filterp = false;	/* Filters based on nunique, nconcordant, or maxminsupport */
static int required_nunique = -1;
static int required_nconcordant = -1;
static int required_maxminsupport = -1;

static int default_count = 1;


static struct option long_options[] = {
  /* Input options */
  {"resolvedir", no_argument, 0, 'D'},	      /* resolve_directions_p */
  {"extents_fwd", required_argument, 0, 'f'}, /* extents_fwd_iit */
  {"extents_rev", required_argument, 0, 'r'}, /* extents_rev_iit */

  {"genebounds", required_argument, 0, 'g'}, /* genebounds_iit */

  {"nunique", required_argument, 0, 'u'},     /* required_nunique, filterp */
  {"nconcordant", required_argument, 0, 'c'}, /* required_nconcordant, filterp */
  {"support", required_argument, 0, 's'}, /* required_maxminsupport, filterp */
  {"noclean", no_argument, 0, 'N'}, /* noclean_p */

  {"excluded", no_argument, 0, 'X'},  /* print_included_p */
  {"runlength", no_argument, 0, 'R'},	/* runlength_p */
  {"endpoints", no_argument, 0, 'E'},	/* endpoints_p */
  {"midpoints", no_argument, 0, 'M'},	/* midpoints_p */
  {"sites", no_argument, 0, 'S'},	/* print_sites_p */

  {"default-count", required_argument, 0, 0}, /* default_count */

  {"dump",no_argument, 0, '9'},	/* dump_p */

  /* Help options */
  {"version", no_argument, 0, 'v'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"spliceclean\n");
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
    fprintf(stdout,"\
Usage: cat <splices> | spliceclean [OPTIONS...]\n\
\n\
Options\n\
  -u, --nunique=INT            Required number of unique reads for splice\n\
  -c, --nconcordant=INT        Required number of concordant reads for splice\n\
  -s, --support=INT            Required maxminsupport for splice\n\
  -N, --noclean                Does no cleaning\n\
  --default-count=INT          Value to use when splice does not provide a count\n\
\n\
Calculation options\n\
  --D, --resolvedir            Resolve directions of conflicting splices\n\
\n\
Output options\n\
  -X, --excluded               Output only excluded splices\n\
  -R, --runlength              Output in runlength format\n\
  -E, --endpoints              Output as endpoints\n\
  -M, --midpoints              Output as midpoints\n\
  -S, --sites                  Output as splice sites\n\
\n\
");
    return;
}


/************************************************************************
 *   Splices
 ************************************************************************/

typedef struct Splice_T *Splice_T;
struct Splice_T {
  Genomicpos_T low;
  Genomicpos_T high;
  int sign;
  long int count;

  char *label;
  char *restofheader;
  char *annotation;

  Genomicpos_T bestprev_ceiling;
  Genomicpos_T bestprev_floor;

  long int total_up;
  long int total_down;

  long int score_ceiling;
  long int score_floor;

  bool boundedp;

  long int mincount;
  long int maxcount;
  bool significantp;

  bool validp;
};


static void
Splice_free (Splice_T *old) {
  if ((*old)->label != NULL) {
    FREE((*old)->label);
  }
  if ((*old)->restofheader != NULL) {
    FREE((*old)->restofheader);
  }
  if ((*old)->annotation != NULL) {
    FREE((*old)->annotation);
  }
  FREE(*old);
  return;
}

static Splice_T
Splice_new (Genomicpos_T low, Genomicpos_T high, int sign, long int count,
	    char *label, char *restofheader, char *annotation) {
  Splice_T new = (Splice_T) MALLOC(sizeof(*new));

  new->low = low;
  new->high = high;
  new->sign = sign;
  new->count = count;
  new->label = label;
  new->restofheader = restofheader;
  new->annotation = annotation;

  new->bestprev_ceiling = 0;
  new->bestprev_floor = (unsigned int) -1;
  new->score_ceiling = 0;
  new->score_floor = 0;

  if (noclean_p == true) {
    new->boundedp = true;
  } else {
    new->boundedp = false;
  }
  new->significantp = false;
  new->validp = false;

  new->mincount = 0;
  new->maxcount = 0;

  return new;
}


static int
Splice_order_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high > y->high) {
    return -1;
  } else if (y->high > x->high) {
    return +1;
  } else {
    return 0;
  }
}


static int
Splice_low_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1; 
  } else {
    return 0;
  }
}

static int
Splice_high_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->high > y->high) {
    return -1;
  } else if (y->high > x->high) {
    return +1;
  } else {
    return 0;
  }
}



static int
Splice_low_descending_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->low > y->low) {
    return -1;
  } else if (y->low > x->low) {
    return +1;
  } else {
    return 0;
  }
}

static int
Splice_high_ascending_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;
  } else {
    return 0;
  }
}

static int
Splice_count_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (y->count >  x->count) {
    return +1;
  } else {
    return 0;
  }
}


static List_T
Splice_sort_order (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_order_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}


static List_T
Splice_sort_low (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_low_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}

static List_T
Splice_sort_high (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_high_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}

static List_T
Splice_sort_low_descending (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_low_descending_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}

static List_T
Splice_sort_high_ascending (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_high_ascending_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}


static List_T
Splice_count_sort (List_T splices) {
  List_T sorted = NULL;
  Splice_T *array;
  int n, i;

  if ((n = List_length(splices)) == 0) {
    return (List_T) NULL;
  } else {
    array = (Splice_T *) List_to_array(splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_count_cmp);
    
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }
    FREE(array);
    return sorted;
  }
}

static void
compute_total_up (List_T splices) {
  long int total_up;
  List_T p;
  Splice_T splice;

  total_up = 0;
  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    total_up += splice->count;
    splice->total_up = total_up;
  }

  return;
}

static void
compute_total_down (List_T splices) {
  long int total_down;
  List_T p;
  Splice_T splice;

  total_down = 0;
  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    total_down += splice->count;
    splice->total_down = total_down;
  }

  return;
}



typedef struct Base_T *Base_T;
struct Base_T {
  Genomicpos_T prevpos;
  Genomicpos_T minextent;
  Genomicpos_T maxextent;
  List_T splices;

  bool usedp;
};

static void
Base_free (Base_T *old) {
  if ((*old)->splices != NULL) {
    List_free(&(*old)->splices);
  }
  FREE(*old);
  return;
}


static Base_T
Base_new () {
  Base_T new = (Base_T) MALLOC(sizeof(*new));

  new->minextent = -1U;
  new->maxextent = -1U;
  new->splices = (List_T) NULL;
  new->usedp = false;

  return new;
}


static void
parse_restofheader (int *nunique, int *nconcordant, int *maxminsupport, char *restofheader) {
  char *p;

  p = restofheader;

  *nunique = *nconcordant = *maxminsupport = -1;

  while (*p != '\0') {
    /* Skip spaces */
    while (*p != '\0' && isspace(*p)) p++;

    if (!strncmp(p,"nunique",strlen("nunique"))) {
      /* nunique */
      while (*p != '\0' && *p != ':') p++;
      p++;
      if (sscanf(p,"%d",&(*nunique)) < 1) {
	fprintf(stderr,"Can't parse value after nunique in %s\n",restofheader);
	abort();
      }

    } else if (!strncmp(p,"nconcordant",strlen("nconcordant"))) {
      /* nconcordant */
  
      while (*p != '\0' && *p != ':') p++;
      p++;
      if (sscanf(p,"%d",&(*nconcordant)) < 1) {
	fprintf(stderr,"Can't parse value after nconcordant in %s\n",restofheader);
	abort();
      }

    } else if (!strncmp(p,"maxminsupport",strlen("maxminsupport"))) {
      /* maxminsupport */
  
      while (*p != '\0' && *p != ':') p++;
      p++;
      if (sscanf(p,"%d",&(*maxminsupport)) < 1) {
	fprintf(stderr,"Can't parse value after maxminsupport in %s\n",restofheader);
	abort();
      }

    }
      
    /* Skip remaining text */
    while (*p != '\0' && !isspace(*p)) p++;

  }

  return;
}


static void
Splice_add (Splice_T splice, Table_T table, char *chr, Genomicpos_T chrpos) {
  Chrom_T chrom;
  Uinttable_T basetable;
  Base_T base;

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((basetable = (Uinttable_T) Table_get(table,(void *) chrom)) == NULL) {
    basetable = Uinttable_new(65522);
    Table_put(table,(void *) chrom,(void *) basetable);
  } else {
    Chrom_free(&chrom);
  }

  if ((base = (Base_T) Uinttable_get(basetable,chrpos)) == NULL) {
    base = Base_new();
    Uinttable_put(basetable,chrpos,(void *) base);
  }

  base->splices = List_push(base->splices,(void *) splice);

  return;
}


/* Empties contents of lines */
static char *
concatenate_lines (List_T lines, int content_size) {
  char *string, *temp;
  List_T l;

  string = (char *) CALLOC(content_size+1,sizeof(char));
  for (l = lines; l; l = List_next(l)) {
    temp = (char *) List_head(l);
    strcat(string,temp);
    FREE(temp);
  }
  
  /* Keep last return
  if (string[content_size-1] == '\n') {
    string[content_size-1] = '\0';
  }
  */

  return string;
}



static void
parse_splice (int *nwarnings, char *spliceline, List_T lines, int content_size,
	      Table_T fwd_low_table, Table_T fwd_high_table,
	      Table_T rev_low_table, Table_T rev_high_table) {
  Genomicpos_T low, high, start, end;
  Splice_T splice;
  char chr[1024], *label, *restofheader, *annotation, *p, *q;
  long int count;
  int len;
  int nunique, nconcordant, maxminsupport;
  bool okayp;

  if (spliceline[0] != '>') {
    fprintf(stderr,"Expected line to start with '>': %s\n",spliceline);
    abort();
  } else {
    p = &(spliceline[1]);      
    if (!isdigit(*p)) {
      q = p;
      if (*nwarnings == MAX_WARNINGS) {
	fprintf(stderr,"More than %d warnings\n",MAX_WARNINGS);
      } else if (*nwarnings < MAX_WARNINGS) {
	fprintf(stderr,"Label for %s does not begin with digit.  Assuming count of %d.\n",
		spliceline,default_count);
      }
      (*nwarnings)++;
      count = default_count;

    } else {
      /* Label must represent a count */
      q = (char *) NULL;
      count = atoi(p);
    }

    /* Skip rest of label */
    len = 0;
    while (*p != '\0' && *p != ' ') {
      len++;
      p++;
    }
    if (*p == ' ') p++;

    if (q == NULL) {
      label = (char *) NULL;
    } else {
      label = (char *) CALLOC(len+1,sizeof(char));
      strncpy(label,q,len);
    }


    /* Get chr part */
    q = p;
    len = 0;
    while (*q != '\0' && *q != ':') {
      q++;
      len++;
    }
    if (*q == '\0') {
      fprintf(stderr,"Can't parse chr part of %s\n",spliceline);
      abort();
    } else {
      strncpy(chr,p,len);
      chr[len] = '\0';
    }

    p = ++q;
    /* p = 16344406..16344388 */
    if (sscanf(p,"%u",&start) != 1) {
      fprintf(stderr,"Can't parse first chrpos in %s\n",spliceline);
      abort();
    }

    /* Advance past first chrpos */
    while (*p != '\0' && isdigit(*p)) p++;
    while (*p != '\0' && !isdigit(*p)) p++;
    /* p = 16344388 */
    if (sscanf(p,"%u",&end) != 1) {
      fprintf(stderr,"Can't parse second chrpos in %s\n",spliceline);
      abort();
    }

    /* Advance past second chrpos */
    while (*p != '\0' && isdigit(*p)) p++;

    /* Look for more content on header */
    q = p;
    while (*q != '\0' && isspace(*q)) q++;
    if (*q == '\0') {
      restofheader = (char *) NULL;
    } else {
      restofheader = (char *) CALLOC(strlen(p)+1,sizeof(char));
      strcpy(restofheader,p);
    }

    /* Make annotation */
    if (lines == NULL) {
      annotation = (char *) NULL;
    } else {
      annotation = concatenate_lines(lines,content_size);
      debug(printf("annotation is %s\n",annotation));
    }

    okayp = true;
    if (filterp == true) {
      if (restofheader == (char *) NULL) {
	fprintf(stderr,"Unable to find filtering information in %s\n",restofheader);
	exit(9);
      } else {
	parse_restofheader(&nunique,&nconcordant,&maxminsupport,restofheader);
	if (required_nunique >= 0) {
	  if (nunique < 0) {
	    fprintf(stderr,"Unable to find filtering information for nunique in %s\n",restofheader);
	    exit(9);
	  } else if (nunique < required_nunique) {
	    okayp = false;
	  }
	}
	if (required_nconcordant >= 0) {
	  if (nconcordant < 0) {
	    fprintf(stderr,"Unable to find filtering information for nconcordant in %s\n",restofheader);
	    exit(9);
	  } else if (nconcordant < required_nconcordant) {
	    okayp = false;
	  }
	}
	if (required_maxminsupport >= 0) {
	  if (maxminsupport < 0) {
	    fprintf(stderr,"Unable to find filtering information for maxminsupport in %s\n",restofheader);
	    exit(9);
	  } else if (maxminsupport < required_maxminsupport) {
	    okayp = false;
	  }
	}
      }
    }

    if (okayp == false) {
      /* Skip */

    } else if (start <= end) {
      low = start;
      high = end;

      splice = Splice_new(low,high,/*sign*/+1,count,label,restofheader,annotation);

      Splice_add(splice,fwd_low_table,chr,low);
      Splice_add(splice,fwd_high_table,chr,high);

    } else {
      low = end;
      high = start;

      splice = Splice_new(low,high,/*sign*/-1,count,label,restofheader,annotation);

      Splice_add(splice,rev_low_table,chr,low);
      Splice_add(splice,rev_high_table,chr,high);

    }

    /* FREE(annotation); -- freed by Splice_free() */
  }

  return;
}


static void
parse_input (Table_T fwd_low_table, Table_T fwd_high_table,
	     Table_T rev_low_table, Table_T rev_high_table) {
  List_T lines = NULL;
  char line[1024000], *copy, *spliceline;
  int content_size;
  int nwarnings = 0;

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '>') {
      if (lines != NULL) {
	debug(printf("Running parse_splice on %d lines\n",List_length(lines)));
	lines = List_reverse(lines);
	lines = List_pop(lines,(void **) &spliceline);
	debug(printf("spliceline is %s\n",spliceline));
	parse_splice(&nwarnings,spliceline,lines,content_size,
		     fwd_low_table,fwd_high_table,rev_low_table,rev_high_table);
	FREE(spliceline);
	List_free(&lines);
      }
      lines = (List_T) NULL;
      content_size = 0;
    } else {
      content_size += strlen(line);
    }
    copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
    strcpy(copy,line);
    debug(printf("Pushing %s",copy));
    lines = List_push(lines,(void *) copy);
  }

  if (lines != NULL) {
    debug(printf("Running parse_splice on %d lines\n",List_length(lines)));
    lines = List_reverse(lines);
    lines = List_pop(lines,(void **) &spliceline);
    debug(printf("spliceline is %s\n",spliceline));
    parse_splice(&nwarnings,spliceline,lines,content_size,
		 fwd_low_table,fwd_high_table,rev_low_table,rev_high_table);
    FREE(spliceline);
    List_free(&lines);
  }

  return;
}



static void
dump_cum (long int *cum, Genomicpos_T chrlength) {
  Genomicpos_T pos;

  for (pos = 0; pos < chrlength; pos++) {
    if (cum[pos] != 0) {
      printf("%u %ld\n",pos,cum[pos]);
    }
  }

  return;

}


/************************************************************************
 *   Dynamic programming
 ************************************************************************/

/* Assumes splices are arranged low to high */
static Splice_T
apply_ceiling (List_T splices, Genomicpos_T ceiling) {
  List_T p;
  Splice_T prevsplice = NULL, splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (splice->high > ceiling) {
      return prevsplice;
    }
    prevsplice = splice;
  }

  return prevsplice;
}

/* Assumes splices are arranged high to low */
static Splice_T
apply_floor (List_T splices, Genomicpos_T floor) {
  List_T p;
  Splice_T prevsplice = NULL, splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (splice->low < floor) {
      return prevsplice;
    }
    prevsplice = splice;
  }

  return prevsplice;
}




static Genomicpos_T
compute_ceilings (Uinttable_T low_basetable) {
  Genomicpos_T ceiling, bestprevpos, prevpos;
  long int bestscore, score;
  int nlookback;
  Genomicpos_T *keys;
  Base_T base, prevbase;
  Splice_T splice, prevsplice;
#ifdef DEBUG2
  Splice_T bestprevsplice;
#endif
  List_T p;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  /* Initialize minbaselow */
  base = (Base_T) Uinttable_get(low_basetable,keys[0]);
  for (p = base->splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    splice->score_ceiling = splice->total_up;
    splice->bestprev_ceiling = 0U;
  }

  for (i = 1; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);

    debug2(printf("At base %u, have %d splices\n",keys[i],List_length(base->splices)));
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      
      bestscore = 0;
      bestprevpos = keys[0];
      debug2(bestprevsplice = NULL);

      prevpos = base->prevpos;
      nlookback = 0;

      ceiling = splice->high /*-1*/;
      debug2(printf("  Splice %u",splice->high));
      while (prevpos >= keys[0] && nlookback < MAXLOOKBACK) {
	prevbase = (Base_T) Uinttable_get(low_basetable,prevpos);
	if ((prevsplice = apply_ceiling(prevbase->splices,ceiling)) != NULL) {
	  debug2(printf(" ... prev:%u, score:%ld+%ld = %ld",
			prevpos,prevsplice->score_ceiling,splice->total_up,
			prevsplice->score_ceiling+splice->total_up));
	    
	  if ((score = prevsplice->score_ceiling + splice->total_up) > bestscore) {
	    debug2(printf("*"));
	    bestscore = score;
	    bestprevpos = prevpos;
#ifdef DEBUG2
	    bestprevsplice = prevsplice;
#endif
	  }
	}

	prevpos = prevbase->prevpos;
	nlookback++;
      }
      debug2(printf("\n"));

      splice->score_ceiling = bestscore;
      splice->bestprev_ceiling = bestprevpos;

#ifdef DEBUG2
      if (bestprevsplice == NULL) {
	printf(" no prevsplice, bestprev is %u, score is %ld\n",
	       bestprevpos,splice->score_ceiling);
      } else {
	printf(" bestprev is pos %u, splice %u, score is %ld\n",
	       bestprevpos,bestprevsplice->high,splice->score_ceiling);
      }
#endif
    }
  }

  /* Get best overall splice */
  bestscore = 0;
  bestprevpos = 0U;

  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      
      if (splice->score_ceiling > bestscore) {
	bestscore = splice->score_ceiling;
	bestprevpos = keys[i];
      }
    }
  }

  FREE(keys);
  return bestprevpos;
}


static Genomicpos_T
compute_floors (Uinttable_T high_basetable) {
  Genomicpos_T floor, bestprevpos, prevpos;
  long int bestscore, score;
  int nlookback;
  Genomicpos_T *keys;
  Base_T base, prevbase;
  Splice_T splice, prevsplice;
#ifdef DEBUG2
  Splice_T bestprevsplice;
#endif
  List_T p;
  int n, i;

  n = Uinttable_length(high_basetable);
  keys = Uinttable_keys(high_basetable,/*sortp*/true);

  /* Initialize maxbasehigh */
  base = (Base_T) Uinttable_get(high_basetable,keys[n-1]);
  for (p = base->splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    splice->score_floor = splice->total_down;
    splice->bestprev_floor = (unsigned int) -1U;
  }

  for (i = n-2; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);

    debug2(printf("At base %u, have %d splices\n",keys[i],List_length(base->splices)));
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      
      bestscore = 0;
      bestprevpos = keys[n-1];
      debug2(bestprevsplice = NULL);

      prevpos = base->prevpos;
      nlookback = 0;

      floor = splice->low /*+1*/;
      debug2(printf("  Splice %u",splice->low));
      while (prevpos <= keys[n-1] && nlookback < MAXLOOKBACK) {
	prevbase = (Base_T) Uinttable_get(high_basetable,prevpos);
	if ((prevsplice = apply_floor(prevbase->splices,floor)) != NULL) {
	  debug2(printf(" ... prev:%u, score:%ld+%ld = %ld",
			prevpos,prevsplice->score_floor,splice->total_down,
			prevsplice->score_floor+splice->total_down));

	  if ((score = prevsplice->score_floor + splice->total_down) > bestscore) {
	    debug2(printf("*"));
	    bestscore = score;
	    bestprevpos = prevpos;
#ifdef DEBUG2
	    bestprevsplice = prevsplice;
#endif
	  }
	}

	prevpos = prevbase->prevpos;
	nlookback++;
      }
      debug2(printf("\n"));

      splice->score_floor = bestscore;
      splice->bestprev_floor = bestprevpos;

#ifdef DEBUG2
      if (bestprevsplice == NULL) {
	printf(" no prevsplice, bestprev is %u, score is %ld\n",
	       bestprevpos,splice->score_floor);
      } else {
	printf(" bestprev is pos %u, splice %u, score is %ld\n",
	       bestprevpos,bestprevsplice->low,splice->score_floor);
      }
#endif
    }
  }

  /* Get best overall splice */
  bestscore = 0;
  bestprevpos = (unsigned int) -1U;

  for (i = n-1; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      
      if (splice->score_floor > bestscore) {
	bestscore = splice->score_floor;
	bestprevpos = keys[i];
      }
    }
  }

  FREE(keys);
  return bestprevpos;
}


static void
traceback_ceilings (Uinttable_T low_basetable, Genomicpos_T prevpos) {
  Genomicpos_T ceiling;
  Splice_T end_splice;
  Base_T base, prevbase;
  Genomicpos_T *keys;
  int n, i;
  
  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  ceiling = (unsigned int) -1U;

  i = n-1;
  while (prevpos > keys[0]) {
    debug3c(printf("traceback from endpos %u, back to %u\n",keys[i],prevpos));
    while (/*startpos*/keys[i] > prevpos) {
      base = (Base_T) Uinttable_get(low_basetable,/*startpos*/keys[i]);
      base->maxextent = ceiling;
      debug3c(printf("At low %u, maxextent is %u\n",/*startpos*/keys[i],ceiling));
      i--;
    }

    prevbase = (Base_T) Uinttable_get(low_basetable,prevpos);
    if ((end_splice = apply_ceiling(prevbase->splices,ceiling)) == NULL) {
      prevpos = keys[0];	/* Ends loop */
    } else {
      ceiling = end_splice->high /*-1*/;
      prevpos = end_splice->bestprev_ceiling;
    }
  }

  debug3c(printf("End of loop\n"));
  while (i >= 0) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    base->maxextent = ceiling;
    debug3c(printf("At low %u, maxextent is %u\n",keys[i],ceiling));
    i--;
  }

  FREE(keys);

  return;
}


static void
traceback_floors (Uinttable_T high_basetable, Genomicpos_T prevpos) {
  Genomicpos_T floor;
  Splice_T start_splice;
  Base_T base, prevbase;
  Genomicpos_T *keys;
  int n, i;

  n = Uinttable_length(high_basetable);
  keys = Uinttable_keys(high_basetable,/*sortp*/true);

  floor = 0U;

  i = 0;
  while (prevpos < keys[n-1]) {
    debug3f(printf("traceback from startpos %u, forward to %u\n",keys[i],prevpos));
    while (/*endpos*/keys[i] < prevpos) {
      base = (Base_T) Uinttable_get(high_basetable,/*endpos*/keys[i]);
      base->minextent = floor;
      debug3f(printf("At high %u, minextent is %u\n",/*endpos*/keys[i],floor));
      i++;
    }

    prevbase = (Base_T) Uinttable_get(high_basetable,prevpos);
    if ((start_splice = apply_floor(prevbase->splices,floor)) == NULL) {
      prevpos = keys[n-1];	/* Ends loop */
    } else {
      floor = start_splice->low /*+1*/;
      prevpos = start_splice->bestprev_floor;
    }
  }

  debug3f(printf("End of loop\n"));
  while (i < n) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    base->minextent = floor;
    debug3f(printf("At high %u, minextent is %u\n",/*endpos*/keys[i],floor));
    i++;
  }

  FREE(keys);

  return;
}


static void
bound_splices (Uinttable_T low_basetable, Uinttable_T high_basetable) {
  Genomicpos_T minextent, maxextent;
  Base_T base_low, base_high;
  Splice_T splice;
  List_T p;
  Genomicpos_T *keys;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  for (i = 0; i < n; i++) {
    base_low = (Base_T) Uinttable_get(low_basetable,keys[i]);
    maxextent = base_low->maxextent;
    for (p = base_low->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      base_high = (Base_T) Uinttable_get(high_basetable,splice->high);
      minextent = base_high->minextent;
      if (splice->low < minextent || splice->high > maxextent) {
	/* printf("\t#%ld %s:%u..%u\n",splice->count,chr,pos,splice->pos); */
      } else {
	splice->boundedp = true;
      }
    }
  }

  FREE(keys);

  return;
}


static void
validate_splices (Uinttable_T low_basetable, Uinttable_T high_basetable) {
  Base_T base, base_low, base_high;
  Splice_T splice;
  List_T p;
  Genomicpos_T *keys;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->boundedp == true && splice->count > 0.05 * splice->maxcount) {
	splice->significantp = true;
	base_low = (Base_T) Uinttable_get(low_basetable,splice->low);
	base_high = (Base_T) Uinttable_get(high_basetable,splice->high);
	base_low->usedp = true;
	base_high->usedp = true;
      }
    }
  }

  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->boundedp == true) {
	base_low = (Base_T) Uinttable_get(low_basetable,splice->low);
	base_high = (Base_T) Uinttable_get(high_basetable,splice->high);
	if (base_low->usedp == true && base_high->usedp == true) {
	  splice->validp = true;
	}
      }
    }
  }

  FREE(keys);

  return;
}



/************************************************************************
 *   Significance
 ************************************************************************/

#if 0

static long int
total_n (List_T splices, int n) {
  long int total = 0;
  Splice_T splice;
  List_T p;
  int i;

  for (p = splices, i = 0; p != NULL && i < n; p = List_next(p), i++) {
    splice = (Splice_T) List_head(p);
    if (splice->validp == true) {
      total += splice->count;
    }
  }

  return total;
}


/* Assumes splices are arranged from high count to low count.  Works
   iteratively. */
static int
number_significant (long int *total_significant, List_T splices, double percentile) {
  List_T p;
  Splice_T prevsplice = NULL, splice;
  long int total, cum;
  int oldn, newn, i;

  oldn = 0;
  newn = List_length(splices);
  cum = 0;

  while (newn != oldn) {
    oldn = newn;

    cum = 0;
    newn = 0;
    total = total_n(splices,oldn);
    for (p = splices, i = 0; p != NULL && i < oldn; p = List_next(p), i++) {
      splice = (Splice_T) List_head(p);
      if (splice->validp == true) {
	if ((double) cum/(double) total < percentile) {
	  newn++;
	}
	cum += splice->count;
      }
    }
  }

  *total_significant = cum;
  return newn;
}

static void
apply_significance (Base_T *bases_low, Genomicpos_T minbaselow, Genomicpos_T maxbaselow,
		    Base_T *bases_high, Genomicpos_T minbasehigh, Genomicpos_T maxbasehigh,
		    double percentile) {
  Genomicpos_T minextent, maxextent, low, high;
  Splice_T splice;
  List_T splices, p;
  long int total_significant;
  int n, i;

  /* Apply percentiles */
  for (low = minbaselow; low <= maxbaselow; low++) {
    if (bases_low[low] != NULL && bases_low[low]->splices != NULL) {
      splices = bases_low[low]->splices = Splice_count_sort(bases_low[low]->splices);

      n = number_significant(&total_significant,splices,percentile);
      bases_low[low]->total_significant = total_significant;

      p = splices;
      i = 0;
      while (i < n) {
	splice = (Splice_T) List_head(p);
	if (splice->validp == true) {
	  splice->significant_low = true;
	  i++;
	}
	p = List_next(p);
      }
    }
  }

  for (high = minbasehigh; high <= maxbasehigh; high++) {
    if (bases_high[high] != NULL && bases_high[high]->splices != NULL) {
      splices = bases_high[high]->splices = Splice_count_sort(bases_high[high]->splices);

      n = number_significant(&total_significant,splices,percentile);
      bases_high[high]->total_significant = total_significant;

      p = splices;
      i = 0;
      while (i < n) {
	splice = (Splice_T) List_head(p);
	if (splice->validp == true) {
	  splice->significant_high = true;
	  i++;
	}
	p = List_next(p);
      }
    }
  }

  return;
}

#endif


/************************************************************************/


#if 0
static void
quantify_splices (Uinttable_T low_basetable) {
  Genomicpos_T low, pos;
  Splice_T splice;
  List_T splices, p;
  long int *counts, level;
  Genomicpos_T *keys;
  int n, i;

  counts = (long int *) CALLOC(chrlength+1,sizeof(long int));

  for (low = minbaselow; low <= maxbaselow; low++) {
    if (bases_low[low] != NULL && (splices = bases_low[low]->splices) != NULL) {
      for (p = splices; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (splice->boundedp == true) {
	  counts[low] += splice->count;
	  /* This is splice->high + 1, because we want to keep cum at splice->pos */
	  counts[splice->high + 1] -= splice->count;
	}
      }
    }
  }

  level = 0;
  for (pos = 0; pos < chrlength; pos++) {
    if (counts[pos] != 0) {
      level += counts[pos];
    }
    counts[pos] = level;
  }

  for (low = minbaselow; low <= maxbaselow; low++) {
    if (bases_low[low] != NULL && (splices = bases_low[low]->splices) != NULL) {
      for (p = splices; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (splice->boundedp == true) {
	  Tally_range(&splice->mincount,&splice->maxcount,counts,splice->low,splice->high);
	}
      }
    }
  }

  FREE(counts);

  return;
}
#endif


static void
dump_splices (Uinttable_T low_basetable, Uinttable_T high_basetable, char *chr) {
  Genomicpos_T minextent, maxextent;
  Genomicpos_T *keys;
  Base_T base;
  List_T p;
  Splice_T splice;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    if (base->splices != NULL) {
      maxextent = base->maxextent;
      printf(">%u maxextent:%u\n",keys[i],maxextent);
      for (p = base->splices; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (splice->sign > 0) {
	  printf("\t%ld %s:%u..%u",splice->count,chr,splice->low,splice->high);
	} else {
	  printf("\t%ld %s:%u..%u",splice->count,chr,splice->high,splice->low);
	}

	printf(" total_up:%ld",splice->total_up);
	printf(" score_ceiling:%ld",splice->score_ceiling);
	printf(" prev_ceiling:%u",splice->bestprev_ceiling);

	if (splice->boundedp == false) {
	  printf(" X");
	} else {
	  printf(" *");
	}
	printf(" %ld %ld",splice->mincount,splice->maxcount);
	printf("\n");
      }
    }
  }
  FREE(keys);

  n = Uinttable_length(high_basetable);
  keys = Uinttable_keys(high_basetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    if (base->splices != NULL) {
      minextent = base->minextent;
      printf("<%u minextent:%u\n",keys[i],minextent);
      for (p = base->splices; p != NULL; p = List_next(p)) {
	splice = (Splice_T) List_head(p);
	if (splice->sign > 0) {
	  printf("\t%ld %s:%u..%u",splice->count,chr,splice->low,splice->high);
	} else {
	  printf("\t%ld %s:%u..%u",splice->count,chr,splice->high,splice->low);
	}

	printf(" total_down:%ld",splice->total_down);
	printf(" score_floor:%ld",splice->score_floor);
	printf(" prev_floor:%u",splice->bestprev_floor);

	if (splice->boundedp == false) {
	  printf(" X");
	} else {
	  printf(" *");
	}
	printf(" %ld %ld",splice->mincount,splice->maxcount);
	printf("\n");
      }
    }
  }
  FREE(keys);

  return;
}


static void
print_splicesites (Uinttable_T low_basetable, char *chr) {
  Genomicpos_T *keys;
  Base_T base;
  List_T p;
  Splice_T splice;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->boundedp == print_included_p) {
	if (splice->sign > 0) {
	  printf(">%ld %s:%u..%u donor\n",splice->count,chr,splice->low,splice->low+1);
	  printf(">%ld %s:%u..%u acceptor\n",splice->count,chr,splice->high-1,splice->high);
	} else {
	  printf(">%ld %s:%u..%u donor\n",splice->count,chr,splice->high,splice->high-1);
	  printf(">%ld %s:%u..%u acceptor\n",splice->count,chr,splice->low+1,splice->low);
	}
#if 0
	if (splice->restofheader == NULL) {
	  printf("\n");
	} else {
	  printf("%s",splice->restofheader);
	  /* printf("\n"); -- should contain linefeed */
	}
	if (splice->annotation != NULL) {
	  printf("%s",splice->annotation);
	  /* printf("\n"); -- not necessary since annotation contains linefeeds */
	}
#endif
      }
    }
  }

  return;
}



static void
check_genebounds (Uinttable_T low_basetable, char *chr, int sign) {
  Genomicpos_T *keys;
  Base_T base;
  List_T p;
  Splice_T splice;
  int n, i;
  int divno;
  int *low_matches, nlowmatches, *high_matches, nhighmatches;

  divno = IIT_divint(genebounds_iit,chr);

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->sign == sign) {
	low_matches = IIT_get_signed_with_divno(&nlowmatches,genebounds_iit,divno,splice->low,splice->low,
						/*sortp*/false,sign);
	high_matches = IIT_get_signed_with_divno(&nhighmatches,genebounds_iit,divno,splice->high,splice->high,
						 /*sortp*/false,sign);
	if (nlowmatches >= 1 && nhighmatches >= 1) {
	  if (nlowmatches > 1) {
	    fprintf(stderr,"Unexpected: multiple genebounds at %s:%u\n",chr,splice->low);
	  } else if (nhighmatches > 1) {
	    fprintf(stderr,"Unexpected: multiple genebounds at %s:%u\n",chr,splice->high);
	  } else if (low_matches[0] != high_matches[0]) {
	    if (sign > 0) {
	      fprintf(stderr,"Splice +%s:%u..%u is between genes\n",chr,splice->low,splice->high);
	    } else {
	      fprintf(stderr,"Splice -%s:%u..%u is between genes\n",chr,splice->high,splice->low);
	    }
	    splice->boundedp = false;
	  } else {
#if 0
	    if (sign > 0) {
	      fprintf(stderr,"Restoring +%s:%u..%u as an alternate splice\n",chr,splice->low,splice->high);
	    } else if (sign < 0) {
	      fprintf(stderr,"Restoring -%s:%u..%u as an alternate splice\n",chr,splice->high,splice->low);
	    }
#endif
	  }
	}
	  
	FREE(high_matches);
	FREE(low_matches);
      }
    }
  }

  FREE(keys);

  return;
}



static void
print_splices (Uinttable_T low_basetable, char *chr) {
  Genomicpos_T *keys;
  Base_T base;
  List_T p;
  Splice_T splice;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->boundedp == print_included_p) {

	if (splice->label != NULL) {
	  printf(">%s",splice->label);
	} else {
	  printf(">%ld",splice->count);
	}

	if (splice->sign > 0) {
	  printf(" %s:%u..%u",chr,splice->low,splice->high);
	} else {
	  printf(" %s:%u..%u",chr,splice->high,splice->low);
	}

	if (splice->restofheader == NULL) {
	  printf("\n");
	} else {
	  printf("%s",splice->restofheader);
	  /* printf("\n"); -- should contain linefeed */
	}

	if (splice->annotation != NULL) {
	  printf("%s",splice->annotation);
	  /* printf("\n"); -- not necessary since annotation contains linefeeds */
	}

      }
    }
  }
  FREE(keys);

  return;
}


static List_T
accumulate_splices (List_T splices, Uinttable_T low_basetable) {
  Genomicpos_T *keys;
  Base_T base;
  List_T p;
  Splice_T splice;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->boundedp == true) {
	splices = List_push(splices,(void *) splice);
      }
    }
  }
  FREE(keys);

  return splices;
}



static void
free_splices (Uinttable_T low_basetable, Uinttable_T high_basetable) {
  Base_T *bases;
  List_T p;
  Splice_T splice;
  int n, i;

  n = Uinttable_length(low_basetable);
  bases = (Base_T *) Uinttable_values(low_basetable);
  for (i = 0; i < n; i++) {
    for (p = bases[i]->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      Splice_free(&splice);
    }
    Base_free(&(bases[i]));
  }
  FREE(bases);

  n = Uinttable_length(high_basetable);
  bases = (Base_T *) Uinttable_values(high_basetable);
  for (i = 0; i < n; i++) {

#if 0
    /* splices already freed from low_basetable */
    for (p = bases[i]->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      Splice_free(&splice);
    }
#endif

    Base_free(&(bases[i]));
  }
  FREE(bases);

  return;
}



static void
add_pairings (long int *cum, Genomicpos_T chrlength, List_T splices) {
  List_T p;
  Splice_T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
#if 0
    if (splice->high + 1 > chrlength) {
      fprintf(stderr,"splice high + 1 is %u, but chrlength is %u\n",splice->high + 1,chrlength);
      abort();
    }
#endif
    if (splice->boundedp == true) {
      cum[splice->low] += splice->count;
      /* This is at splice->high + 1, because we want to keep cum at splice->high */
      cum[splice->high + 1] -= splice->count;
    }
  }

  return;
}


static void
add_midpoints (long int *cum, Genomicpos_T chrlength, List_T splices) {
  List_T p;
  Splice_T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (splice->boundedp == true) {
      cum[(splice->low + splice->high)/2] += splice->count;
    }
  }

  return;
}

static void
add_endpoints (long int *low_points, long int *high_points, Genomicpos_T chrlength, List_T splices) {
  List_T p;
  Splice_T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (splice->boundedp == true) {
      low_points[splice->low] += splice->count;
      high_points[splice->high] += splice->count;
    }
  }

  return;
}


static void
print_runlengths (Uinttable_T low_basetable, char *chr) {
  List_T p;
  Splice_T splice;
  Genomicpos_T chrlength = 0, lastpos, pos;
  long int *cum, level;
  Genomicpos_T *keys;
  Base_T base;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->high + 1 > chrlength) {
	chrlength = splice->high + 1;
      }
    }
  }

  cum = (long int *) CALLOC(chrlength+1,sizeof(long int));
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    add_pairings(cum,chrlength,base->splices);
  }

  FREE(keys);

  /* dump_cum(cum,chrlength); */

  /* Print runlengths */
#if 0
  lastpos = 0U;
  level = 0;
  for (pos = 0; pos < chrlength; pos++) {
    if (cum[pos] != 0) {
      if (lastpos != 0) {
	if (level == 0) {
	  /* Marks intergenic regions */
	  printf(">0 %s:%u\n",chr,lastpos);
	} else {
	  printf(">%d %s:%u..%u\n",level,chr,lastpos,pos-1);
	}
      }
      level += cum[pos];
      lastpos = pos;
    }
  }

  if (lastpos != 0) {
    if (level == 0) {
      /* printf(">0 %s:%u\n",chr,lastpos); */
    } else {
      fprintf(stderr,"Ended with a non-zero level\n");
      abort();
    }
  }
#else
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
#endif


  FREE(cum);

  return;
}


static void
print_midpoints (Uinttable_T low_basetable, char *chr) {
  List_T p;
  Splice_T splice;
  Genomicpos_T chrlength = 0, pos;
  long int *cum;
  Genomicpos_T *keys;
  Base_T base;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->high + 1 > chrlength) {
	chrlength = splice->high + 1;
      }
    }
  }

  cum = (long int *) CALLOC(chrlength+1,sizeof(long int));
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    add_midpoints(cum,chrlength,base->splices);
  }

  FREE(keys);

  /* dump_cum(cum,chrlength); */

  /* Print midpoints */
  for (pos = 1; pos <= chrlength; pos++) {
    if (cum[pos] != 0) {
      /* printf("cum at pos %u is %d\n",pos,cum[pos]); */
      printf(">%ld %s:%u\n",cum[pos],chr,pos);
    }
  }

  FREE(cum);

  return;
}


static void
print_endpoints (Uinttable_T low_basetable, char *chr) {
  List_T p;
  Splice_T splice;
  Genomicpos_T chrlength = 0, pos;
  long int *low_points, *high_points;
  Genomicpos_T *keys;
  Base_T base;
  int n, i;

  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);
  
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      if (splice->high + 1 > chrlength) {
	chrlength = splice->high + 1;
      }
    }
  }

  low_points = (long int *) CALLOC(chrlength+1,sizeof(long int));
  high_points = (long int *) CALLOC(chrlength+1,sizeof(long int));
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    add_endpoints(low_points,high_points,chrlength,base->splices);
  }

  FREE(keys);

  /* dump_cum(cum,chrlength); */

  /* Print midpoints */
  for (pos = 1; pos <= chrlength; pos++) {
    if (low_points[pos] != 0) {
      printf(">%ld %s:%u up\n",low_points[pos],chr,pos);
    }
    if (high_points[pos] != 0) {
      printf(">%ld %s:%u down\n",high_points[pos],chr,pos);
    }
  }

  FREE(high_points);
  FREE(low_points);

  return;
}



static void
clean_splices (Uinttable_T low_basetable, Uinttable_T high_basetable) {
  Genomicpos_T prevpos;
  Genomicpos_T *keys;
  Base_T base;
  int n, i;
	       
  /* Prepare splices and base links */
  n = Uinttable_length(low_basetable);
  keys = Uinttable_keys(low_basetable,/*sortp*/true);

  prevpos = 0U;
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    base->splices = Splice_sort_high_ascending(base->splices);
    compute_total_up(base->splices);

    base->prevpos = prevpos;
    prevpos = keys[i];
  }
  FREE(keys);

  n = Uinttable_length(high_basetable);
  keys = Uinttable_keys(high_basetable,/*sortp*/true);
  prevpos = (unsigned int) -1U;
  for (i = n-1; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    base->splices = Splice_sort_low_descending(base->splices);
    compute_total_down(base->splices);

    base->prevpos = prevpos;
    prevpos = keys[i];
  }
  FREE(keys);


  prevpos = compute_ceilings(low_basetable);
  traceback_ceilings(low_basetable,prevpos);
  
  prevpos = compute_floors(high_basetable);
  traceback_floors(high_basetable,prevpos);

  bound_splices(low_basetable,high_basetable);

  /* 
  quantify_splices(bases_low,minbaselow,maxbaselow,chrlength);
  validate_splices(bases_low,minbaselow,maxbaselow,bases_high);
  apply_significance(bases_low,minbaselow,maxbaselow,
		     bases_high,minbasehigh,maxbasehigh,percentile);
  */

  return;
}


static int
strongest (Splice_T fwd_splice, Splice_T rev_splice, 
	   long int *extents_fwd, long int *extents_rev,
	   Splice_T *array, int a, int b, int nsplices, char *chr,
	   Genomicpos_T chrlength) {
  long int fwd_max, rev_max;
  int i;

  fwd_max = Tally_maxcount(extents_fwd,chrlength,fwd_splice->low,fwd_splice->high);
  rev_max = Tally_maxcount(extents_rev,chrlength,rev_splice->low,rev_splice->high);
  if (fwd_max > 10*rev_max) {
    fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld)\n",
	    chr,fwd_splice->low,fwd_splice->high,fwd_max,
	    chr,rev_splice->high,rev_splice->low,rev_max);
    return +1;
  } else if (rev_max > 10*fwd_max) {
    fprintf(stderr,"Choose -%s:%u..%u (%ld) over +%s:%u..%u (%ld)\n",
	    chr,rev_splice->high,rev_splice->low,rev_max,
	    chr,fwd_splice->low,fwd_splice->high,fwd_max);
    return -1;
  } else {
    i = a-1;
    while (i >= 0 && array[i]->high + 500 > fwd_splice->low && array[i]->high > rev_splice->low) {
      fwd_max = Tally_maxcount(extents_fwd,chrlength,array[i]->low,array[i]->high);
      rev_max = Tally_maxcount(extents_rev,chrlength,array[i]->low,array[i]->high);
      if (fwd_max > 10*rev_max) {
	fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld), secondary left\n",
		chr,fwd_splice->low,fwd_splice->high,fwd_max,
		chr,rev_splice->high,rev_splice->low,rev_max);
	return +1;
      } else if (rev_max > 10*fwd_max) {
	fprintf(stderr,"Choose -%s:%u..%u (%ld) over +%s:%u..%u (%ld), secondary left\n",
		chr,rev_splice->high,rev_splice->low,rev_max,
		chr,fwd_splice->low,fwd_splice->high,fwd_max);
	return -1;
      }
      i--;
    }

    i = b+1;
    while (i < nsplices && array[i]->low < fwd_splice->high + 500 && array[i]->low < rev_splice->high + 500) {
      fwd_max = Tally_maxcount(extents_fwd,chrlength,array[i]->low,array[i]->high);
      rev_max = Tally_maxcount(extents_rev,chrlength,array[i]->low,array[i]->high);
      if (fwd_max > 10*rev_max) {
	fprintf(stderr,"Choose +%s:%u..%u (%ld) over -%s:%u..%u (%ld), secondary right\n",
		chr,fwd_splice->low,fwd_splice->high,fwd_max,
		chr,rev_splice->high,rev_splice->low,rev_max);
	return +1;
      } else if (rev_max > 10*fwd_max) {
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


static void
resolve_splices (Table_T fwd_low_table, Table_T rev_low_table,
		 Chrom_T *fwd_keys, int fwd_n, Chrom_T *rev_keys, int rev_n) {
  int i, j, a, b, nsplices;
  int cmp, result;
  Chrom_T chrom;
  Genomicpos_T chrlength, fwd_chrlength, rev_chrlength;
  long int *extents_fwd, *extents_rev;
  List_T splices;
  Splice_T splice1, splice2, *array;
  char *chr;

  i = j = 0;
  while (i < fwd_n && j < rev_n) {
    if ((cmp = Chrom_compare_chrom(&(fwd_keys[i]),&(rev_keys[j]))) < 0) {
      i++;
    } else if (cmp > 0) {
      j++;
    } else {
      chrom = fwd_keys[i];
      chr = Chrom_string(chrom);
      fprintf(stderr,"Processing chr %s...",chr);

      fwd_chrlength = IIT_divlength(extents_fwd_iit,chr);
      rev_chrlength = IIT_divlength(extents_rev_iit,chr);
      chrlength = (fwd_chrlength > rev_chrlength) ? fwd_chrlength : rev_chrlength;

      extents_fwd = (long int *) CALLOC(chrlength+1,sizeof(long int)); /* 1 */
      Tally_store_runlength(extents_fwd,extents_fwd_iit,chr,/*coordstart*/1,/*coordend*/chrlength);

      extents_rev = (long int *) CALLOC(chrlength+1,sizeof(long int)); /* 1 */
      Tally_store_runlength(extents_rev,extents_rev_iit,chr,/*coordstart*/1,/*coordend*/chrlength);

      splices = accumulate_splices((List_T) NULL,Table_get(fwd_low_table,(void *) chrom));
      splices = accumulate_splices(splices,Table_get(rev_low_table,(void *) chrom));

      nsplices = List_length(splices);
      array = (Splice_T *) List_to_array(splices,NULL);
      fprintf(stderr,"%d splices\n",nsplices);

      qsort(array,nsplices,sizeof(Splice_T),Splice_low_cmp);
      for (a = 0; a < nsplices; a++) {
	splice1 = array[a];
	for (b = a+1; b < nsplices; b++) {
	  splice2 = array[b];
	  if (splice2->low > splice1->low + 20) {
	    b = nsplices;	/* End loop */

	  } else if (splice1->sign == splice2->sign) {
	    /* Don't resolve between splices on the same strand */

	  } else if (splice1->sign > 0) {
	    if ((result = strongest(splice1,splice2,extents_fwd,extents_rev,array,a,b,nsplices,chr,
				    chrlength)) > 0) {
	      splice2->boundedp = false;
	    } else if (result < 0) {
	      splice1->boundedp = false;
	    }
	      
	  } else {
	    if ((result = strongest(splice2,splice1,extents_fwd,extents_rev,array,a,b,nsplices,chr,
				    chrlength)) > 0) {
	      splice1->boundedp = false;
	    } else if (result < 0) {
	      splice2->boundedp = false;
	    }
	  }
	}
      }

      qsort(array,nsplices,sizeof(Splice_T),Splice_high_cmp);
      for (a = 0; a < nsplices; a++) {
	splice1 = array[a];
	for (b = a+1; b < nsplices; b++) {
	  splice2 = array[b];
	  if (splice2->high < splice1->high - 20) {
	    b = nsplices;	/* End loop */

	  } else if (splice1->sign == splice2->sign) {
	    /* Don't resolve between splices on the same strand */

	  } else if (splice1->sign > 0) {
	    if ((result = strongest(splice1,splice2,extents_fwd,extents_rev,array,a,b,nsplices,chr,
				    chrlength)) > 0) {
	      splice2->boundedp = false;
	    } else if (result < 0) {
	      splice1->boundedp = false;
	    }

	  } else {
	    if ((result = strongest(splice2,splice1,extents_fwd,extents_rev,array,a,b,nsplices,chr,
				    chrlength)) > 0) {
	      splice1->boundedp = false;
	    } else if (result < 0) {
	      splice2->boundedp = false;
	    }
	  }
	}
      }

      FREE(array);
      List_free(&splices);
      FREE(extents_rev);
      FREE(extents_fwd);
      i++;
      j++;
    }
  }

  return;
}


/* Usage: cat splices | spliceclean  */

int
main (int argc, char *argv[]) {
  Table_T fwd_low_table, fwd_high_table, rev_low_table, rev_high_table;
  Uinttable_T low_basetable, high_basetable;
  char *chr;
  Chrom_T *fwd_keys, *rev_keys, chrom;
  int fwd_n, rev_n, i;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"Df:r:g:u:c:s:NXREMS9v?",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"default-count")) {
	default_count = atoi(optarg);
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run '<program> --help'",long_name);
	exit(9);
      }
      break;

    case 'D': resolve_directions_p = true; break;

    case 'f':
      if ((extents_fwd_iit = IIT_read(optarg,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"Could not open extents_fwd IIT file %s\n",optarg);
	exit(9);
      }
      break;
    case 'r':
      if ((extents_rev_iit = IIT_read(optarg,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"Could not open extents_rev IIT file %s\n",optarg);
	exit(9);
      }
      break;

    case 'g':
      if ((genebounds_iit = IIT_read(optarg,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"Could not open genebounds IIT file %s\n",optarg);
	exit(9);
      }
      break;

    case 'u': required_nunique = atoi(optarg); filterp = true; break;
    case 'c': required_nconcordant = atoi(optarg); filterp = true; break;
    case 's': required_maxminsupport = atoi(optarg); filterp = true; break;
    case 'N': noclean_p = true; break;
    case 'X': print_included_p = false; break;
    case 'R': runlength_p = true; break;
    case 'E': endpoints_p = true; break;
    case 'M': midpoints_p = true; break;
    case 'S': print_sites_p = true; break;
    case '9': dump_p = true; break;
    case 'v': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;


  fwd_low_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  fwd_high_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  rev_low_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  rev_high_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  fprintf(stderr,"Reading input...");
  parse_input(fwd_low_table,fwd_high_table,rev_low_table,rev_high_table);
  fprintf(stderr,"done\n");


  /* Get keys */
  fwd_n = Table_length(fwd_low_table);
  fwd_keys = (Chrom_T *) Table_keys(fwd_low_table,NULL);
  qsort(fwd_keys,fwd_n,sizeof(Chrom_T),Chrom_compare_chrom);

  rev_n = Table_length(rev_low_table);
  rev_keys = (Chrom_T *) Table_keys(rev_low_table,NULL);
  qsort(rev_keys,rev_n,sizeof(Chrom_T),Chrom_compare_chrom);


  if (noclean_p == false) {
    /* Clean forward strand */
    for (i = 0; i < fwd_n; i++) {
      chrom = fwd_keys[i];
      chr = Chrom_string(chrom);
      low_basetable = Table_get(fwd_low_table,(void *) chrom);
      high_basetable = Table_get(fwd_high_table,(void *) chrom);
      clean_splices(low_basetable,high_basetable);
    }

    /* Clean reverse strand */
    for (i = 0; i < rev_n; i++) {
      chrom = rev_keys[i];
      chr = Chrom_string(chrom);
      low_basetable = Table_get(rev_low_table,(void *) chrom);
      high_basetable = Table_get(rev_high_table,(void *) chrom);
      clean_splices(low_basetable,high_basetable);
    }
  }


  /* Resolve */
  if (resolve_directions_p == true) {
    if (extents_fwd_iit == NULL || extents_rev_iit == NULL) {
      fprintf(stderr,"To resolve directions (-D), need to specify extents_fwd (-f) and extents_rev (-r)\n");
    } else {
      resolve_splices(fwd_low_table,rev_low_table,fwd_keys,fwd_n,rev_keys,rev_n);
      IIT_free(&extents_rev_iit);
      IIT_free(&extents_fwd_iit);
    }
  }

  /* Assumes that noclean_p is true */
  if (genebounds_iit != NULL) {
    /* Check genebounds, fwd */
    for (i = 0; i < fwd_n; i++) {
      chrom = fwd_keys[i];
      chr = Chrom_string(chrom);
      low_basetable = Table_get(fwd_low_table,(void *) chrom);
      check_genebounds(low_basetable,chr,/*sign*/+1);
    }

    /* Check genebounds, rev */
    for (i = 0; i < rev_n; i++) {
      chrom = rev_keys[i];
      chr = Chrom_string(chrom);
      low_basetable = Table_get(rev_low_table,(void *) chrom);
      check_genebounds(low_basetable,chr,/*sign*/-1);
    }
    IIT_free(&genebounds_iit);
  }


  /* Print forward splices */
  for (i = 0; i < fwd_n; i++) {
    chrom = fwd_keys[i];
    chr = Chrom_string(chrom);
    low_basetable = Table_get(fwd_low_table,(void *) chrom);
    high_basetable = Table_get(fwd_high_table,(void *) chrom);

    if (dump_p == true) {
      dump_splices(low_basetable,high_basetable,chr);
    } else if (runlength_p == true) {
      print_runlengths(low_basetable,chr);
    } else if (endpoints_p == true) {
      print_endpoints(low_basetable,chr);
    } else if (midpoints_p == true) {
      print_midpoints(low_basetable,chr);
    } else if (print_sites_p == true) {
      print_splicesites(low_basetable,chr);
    } else {
      print_splices(low_basetable,chr);
    }
    free_splices(low_basetable,high_basetable);

    Uinttable_free(&low_basetable);
    Uinttable_free(&high_basetable);
    Chrom_free(&chrom);
  }

  FREE(fwd_keys);

  Table_free(&fwd_low_table);
  Table_free(&fwd_high_table);


  /* Print reverse splices */
  for (i = 0; i < rev_n; i++) {
    chrom = rev_keys[i];
    chr = Chrom_string(chrom);
    low_basetable = Table_get(rev_low_table,(void *) chrom);
    high_basetable = Table_get(rev_high_table,(void *) chrom);
    
    if (dump_p == true) {
      dump_splices(low_basetable,high_basetable,chr);
    } else if (runlength_p == true) {
      print_runlengths(low_basetable,chr);
    } else if (endpoints_p == true) {
      print_endpoints(low_basetable,chr);
    } else if (midpoints_p == true) {
      print_midpoints(low_basetable,chr);
    } else if (print_sites_p == true) {
      print_splicesites(low_basetable,chr);
    } else {
      print_splices(low_basetable,chr);
    }
    free_splices(low_basetable,high_basetable);

    Uinttable_free(&low_basetable);
    Uinttable_free(&high_basetable);
    Chrom_free(&chrom);
  }

  FREE(rev_keys);
  
  Table_free(&rev_low_table);
  Table_free(&rev_high_table);

  return 0;
}



