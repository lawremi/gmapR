static char rcsid[] = "$Id: gsnap_splices.c 46994 2011-09-12 17:44:59Z twu $";
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
#include "complement.h"
#include "chrom.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "table.h"
#include "uinttable.h"
#include "datadir.h"
#include "getopt.h"

#ifdef SAM_INPUT
#include "samflags.h"
#include "samread.h"
#include "genome.h"
#elif defined(BAM_INPUT)
#include "samflags.h"		/* For flags */
#include "bamread.h"
#include "samread.h"
#include "genome.h"
#else
#include "gsnapread.h"
#endif


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


static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *splicesites_file = (char *) NULL;
static IIT_T splicesites_iit = NULL;
static int donor_typeint;	/* for splicesites_iit */
static int acceptor_typeint;	/* for splicesites_iit */

static bool trust_sam_p = true;
static bool need_canonical_p = false;

static bool need_concordant_p = false;
static bool uniquep = false;

static int mincount = 1;
static int minsupport = 1;

static Genomicpos_T shortsplicedist = 200000;


typedef enum {NONCONCORDANT, DONOR_SIDE, ACCEPTOR_SIDE} Extension_T;


typedef struct Splice_T *Splice_T;
struct Splice_T {
  int count;
  int nunique;
  int maxminsupport;

  Genomicpos_T donorpos;
  Genomicpos_T acceptorpos;

  int donor_extensions;		/* For paired-end reads */
  int acceptor_extensions;	/* For paired-end reads */

  bool canonicalp;
  char donor1;
  char donor2;
  char acceptor1;
  char acceptor2;
};

static void
Splice_free (Splice_T *old) {
  FREE(*old);
  return;
}

static int
Splice_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

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

static List_T
Splice_add_at_donor (List_T list, Genomicpos_T donorpos, Genomicpos_T acceptorpos,
		     Extension_T extension, int nhits, int support1, int support2,
		     bool canonicalp, char donor1, char donor2, char acceptor1, char acceptor2) {
  List_T p;
  Splice_T splice;
  int minsupport;

  if (support1 < support2) {
    minsupport = support1;
  } else {
    minsupport = support2;
  }

  debug1(printf("List has %d elements currently\n",List_length(p)));
  for (p = list; p != NULL; p = List_next(p)) {
    splice = List_head(p);
    debug1(printf("  Comparing acceptorpos %u with list acceptorpos %u",
		  acceptorpos,splice->acceptorpos));
    if (splice->acceptorpos == acceptorpos) {
      splice->count += 1;
      if (extension == DONOR_SIDE) {
	splice->donor_extensions += 1;
      } else if (extension == ACCEPTOR_SIDE) {
	splice->acceptor_extensions += 1;
      }
      if (nhits == 1) {
	splice->nunique += 1;
      }
      if (minsupport > splice->maxminsupport) {
	splice->maxminsupport = minsupport;
      }
      debug1(printf("  equal, so count is now %d\n",splice->count));
      return list;
    }
    debug1(printf("\n"));
  }

  /* Not found, so add to list */
  debug1(printf("Acceptorpos not found, so creating new splice\n"));
  splice = (Splice_T) MALLOC(sizeof(*splice));
  splice->count = 1;
  splice->nunique = 0;
  if (extension == DONOR_SIDE) {
    splice->donor_extensions = 1;
    splice->acceptor_extensions = 0;
  } else if (extension == ACCEPTOR_SIDE) {
    splice->donor_extensions = 0;
    splice->acceptor_extensions = 1;
  } else {
    splice->donor_extensions = 0;
    splice->acceptor_extensions = 0;
  }
  if (nhits == 1) {
    splice->nunique += 1;
  }
  splice->maxminsupport = minsupport;
  splice->donorpos = donorpos;
  splice->acceptorpos = acceptorpos;
  splice->canonicalp = canonicalp;
  splice->donor1 = donor1;
  splice->donor2 = donor2;
  splice->acceptor1 = acceptor1;
  splice->acceptor2 = acceptor2;

  return List_push(list,splice);
}


typedef struct Site_T *Site_T;
struct Site_T {
  Genomicpos_T chrpos;
  List_T splices;
};


static void
Site_free (Site_T *old) {
  List_T p;
  Splice_T splice;

  for (p = (*old)->splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    Splice_free(&splice);
  }
  List_free(&(*old)->splices);

  FREE(*old);
  return;
}


static Site_T
Site_new (Genomicpos_T chrpos) {
  Site_T new = (Site_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->splices = (List_T) NULL;

  return new;
}



static void
print_splicesite_labels (char *chr, Genomicpos_T splicesitelow, Genomicpos_T splicesitehigh,
			 IIT_T splicesites_iit, int typeint) {
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  /* fprintf(stderr,"Looking for %u..%u\n",splicesitelow,splicesitehigh); */
  splicesites = IIT_get_exact_multiple(&nsplicesites,splicesites_iit,chr,
				       splicesitelow,splicesitehigh,typeint);

  if (nsplicesites > 0) {
    label = IIT_label(splicesites_iit,splicesites[0],&allocp);
    printf("%s",label);
    if (allocp) FREE(label);

    for (i = 1; i < nsplicesites; i++) {
      label = IIT_label(splicesites_iit,splicesites[i],&allocp);
      printf("|%s",label);
      if (allocp) FREE(label);
    }
    FREE(splicesites);
  }

  printf("\n");

  return;
}


static void
Site_print (Site_T this, char *chr, IIT_T splicesites_iit, 
	    int donor_typeint, int acceptor_typeint) {
  Splice_T splice, *array;
  int n, i;
  
  n = List_length(this->splices);
  if (n > 0) {
    array = (Splice_T *) List_to_array(this->splices,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_cmp);
    for (i = 0; i < n; i++) {
      splice = array[i];
      if (splice->count >= mincount && splice->maxminsupport >= minsupport) {
	printf(">%d %s:%u..%u nunique:%d",
	       splice->count,chr,splice->donorpos,splice->acceptorpos,
	       splice->nunique);
	printf(" nconcordant:%d",splice->donor_extensions + splice->acceptor_extensions);
	printf(" maxminsupport:%d",splice->maxminsupport);
	if (splice->canonicalp == false) {
	  printf(" %c%c-%c%c",splice->donor1,splice->donor2,splice->acceptor1,splice->acceptor2);
	}
	printf("\n");

	if (splicesites_iit != NULL) {
	  if (splice->donorpos < splice->acceptorpos) {
	    print_splicesite_labels(chr,splice->donorpos,splice->donorpos+1U,splicesites_iit,
				    donor_typeint);
	    print_splicesite_labels(chr,splice->acceptorpos-1U,splice->acceptorpos,splicesites_iit,
				    acceptor_typeint);
	  } else if (splice->donorpos > splice->acceptorpos) {
	    print_splicesite_labels(chr,splice->donorpos-1U,splice->donorpos,splicesites_iit,
				    donor_typeint);
	    print_splicesite_labels(chr,splice->acceptorpos,splice->acceptorpos+1U,splicesites_iit,
				    acceptor_typeint);

	  } else {
	    fprintf(stderr,"Error: splice donorpos and acceptorpos are equal at %s:%u\n",
		    chr,splice->donorpos);
	    exit(9);
	  }
	}
      }
    }
    FREE(array);
  }

  return;
}


static void
Site_add_at_donor (Table_T chr_table, char *chr, Genomicpos_T donorpos, Genomicpos_T acceptorpos,
		   Extension_T extension, int nhits, int support1, int support2, bool canonicalp,
		   char donor1, char donor2, char acceptor1, char acceptor2) {
  Site_T site;
  Uinttable_T sitetable;
  Chrom_T chrom;
  
  debug1(printf("Called Site_add_at_donor with chr %s, donorpos %u, acceptorpos %u\n",
		chr,donorpos,acceptorpos));

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((sitetable = (Uinttable_T) Table_get(chr_table,(void *) chrom)) == NULL) {
    debug1(printf("Made new sitetable for chr %s\n",chr));
    sitetable = Uinttable_new(65522); /* estimate 65522 splice sites per chromosome */
    Table_put(chr_table,(void *) chrom,(void *) sitetable);
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
  site->splices = Splice_add_at_donor(site->splices,donorpos,acceptorpos,extension,nhits,
				      support1,support2,canonicalp,donor1,donor2,acceptor1,acceptor2);

  return;
}


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"splicesites", required_argument, 0, 's'}, /* splicesites_iit */
  {"pairmax", required_argument, 0, 'p'},     /* shortsplicedist */

#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif

  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* uniquep */
  {"canonical", required_argument, 0, 'N'}, /* need_canonical_p */

  {"mincount", required_argument, 0, 'n'}, /* mincount */
  {"minsupport", required_argument, 0, 'w'}, /* minsupport */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
#ifdef SAM_INPUT
  fprintf(stdout,"SAM_SPLICES\n");
#elif defined(BAM_INPUT)
  fprintf(stdout,"BAM_SPLICES\n");
#else
  fprintf(stdout,"GSNAP_SPLICES\n");
#endif
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: cat <GSNAP output> | gsnap_splices -d <genome> [OPTIONS...]\n\
\n\
Input options\n\
  -C, --concordant=INT           Concordant hits only (0=no [default], 1=yes)\n\
  -U, --unique=INT               Unique hits only (0=no [default], 1=yes)\n\
  -n, --mincount=INT             Minimum observed count (default 1)\n\
  -w, --minsupport=INT           Threshold for maxminsupport (default 1)\n\
\n\
Output options\n\
  -9, --dump                     Dump graph\n\
  -N, --canonical=INT            Canonical introns (GT-AG, GC-AG, or AT-AC) only (0=no [default], 1=yes)\n\
\n\
");
    return;
}


/************************************************************************
 *   Output
 ************************************************************************/

static void
print_splices (Uinttable_T sitetable, char *chr, IIT_T splicesites_iit,
	       int donor_typeint, int acceptor_typeint) {
  Genomicpos_T *keys;
  int n, i;
  Site_T site;

  n = Uinttable_length(sitetable);
  keys = (Genomicpos_T *) Uinttable_keys(sitetable,/*sortp*/true);
  for (i = 0; i < n; i++) {
    site = Uinttable_get(sitetable,keys[i]);
    if (site == NULL) {
      fprintf(stderr,"key is %u, value is NULL\n",keys[i]);
      abort();
    } else {
      Site_print(site,chr,splicesites_iit,donor_typeint,acceptor_typeint);
      Site_free(&site);
    }
  }
  FREE(keys);

  return;
}


/************************************************************************
 *   Input
 ************************************************************************/

#ifndef BAM_INPUT
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
#endif


#if defined(SAM_INPUT) || defined(BAM_INPUT)

static char complCode[128] = COMPLEMENT_LC;

static char
find_strand (bool *canonicalp, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
	     Genomicpos_T firstpos, Genomicpos_T secondpos, char *chr,
	     Genome_T genome, IIT_T chromosome_iit, char *auxinfo) {
  Chrnum_T chrnum;
  Genomicpos_T chroffset;
  char nt1, nt2, nt3, nt4;
  char truestrand;

  if (auxinfo != NULL && (truestrand = Samread_splice_strand(auxinfo)) != ' ') {
    if (trust_sam_p == true) {
      *canonicalp = true;
      *donor1 = *donor2 = *acceptor1 = *acceptor2 = ' ';
      return truestrand;
    } else {
      chrnum = IIT_find_one(chromosome_iit,chr);
      chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum)) - 1U;

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

      return truestrand;
    }

  } else if (chromosome_iit == NULL) {
    fprintf(stderr,"Strand is not present in auxinfo %s\n",auxinfo);
    fprintf(stderr,"To determine strand, need to provide index file with -d flag\n");
    exit(9);

  } else {
    chrnum = IIT_find_one(chromosome_iit,chr);
    chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum)) - 1U;

    /* Look at genome inside of firstpos and secondpos to determine truestrand */
    nt1 = Genome_get_char(genome,chroffset+firstpos+1);
    nt2 = Genome_get_char(genome,chroffset+firstpos+2);
    nt3 = Genome_get_char(genome,chroffset+secondpos-2);
    nt4 = Genome_get_char(genome,chroffset+secondpos-1);

    debug(printf("Got splice from %u to %u\n",firstpos,secondpos));
    debug(printf("Dinucleotides are %c%c to %c%c\n",nt1,nt2,nt3,nt4));

    if (nt1 == 'G' && (nt2 == 'T' || nt2 == 'C') && nt3 == 'A' && nt4 == 'G') {
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
      return '+';
    } else if (nt1 == 'C' && nt2 == 'T' && (nt3 == 'A' || nt3 == 'G') && nt4 == 'C') {
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
      return '-';
    } else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A' && nt4 == 'C') {
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
      return '+';
    } else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A' && nt4 == 'T') {
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
      return '-';
    } else {
      /* In GSNAP, will want to output sense information in SAM output. */
#if 0
      fprintf(stderr,"Splice %s:%u..%u is not (semi-)canonical: %c%c...%c%c.  Cannot determine sense.\n",
	      chr,firstpos,secondpos,nt1,nt2,nt3,nt4);
#endif
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = false;
      return ' ';
    }
  }
}



static void
add_splice (Genomicpos_T firstpos, Genomicpos_T secondpos,
	    char *acc, char *chr, int nhits, int support1, int support2,
	    Table_T chr_fwd_table, Table_T chr_rev_table, char *auxinfo,
	    Genome_T genome, IIT_T chromosome_iit, bool pairedp, bool concordantp) {
  /* truestrand can be reported by SAM, even though it is not canonical */
  char truestrand;
  bool canonicalp;
  char donor1, donor2, acceptor1, acceptor2;
  Genomicpos_T donorpos, acceptorpos;
  Extension_T extension;

  truestrand = find_strand(&canonicalp,&donor1,&donor2,&acceptor1,&acceptor2,
			   firstpos,secondpos,chr,genome,chromosome_iit,auxinfo);

  if (truestrand == '+') {
    donorpos = firstpos;
    acceptorpos = secondpos;
  } else if (truestrand == '-') {
    donorpos = secondpos;
    acceptorpos = firstpos;
  } else {
    /* donorpos = firstpos; */
    /* acceptorpos = secondpos; */
  }

  debug0(printf("Got splice from %u to %u\n",donorpos,acceptorpos));
  if (truestrand == '+' && acceptorpos > donorpos + shortsplicedist) {
    /* Distance too far */
  } else if (truestrand == '-' && donorpos > acceptorpos + shortsplicedist) {
    /* Distance too far */
  } else if (truestrand == ' ' && secondpos > firstpos + shortsplicedist) {
    /* Distance too far */
  } else {
    debug(printf("%c%s:%u..%u\n",truestrand,chr,donorpos,acceptorpos));
    
    if (concordantp == false) {
      extension = NONCONCORDANT;
    } else {
      extension = DONOR_SIDE; /* Not distinguishing donor and acceptor sides */
    }
    
    if (truestrand == '+') {
      Site_add_at_donor(chr_fwd_table,chr,donorpos,acceptorpos,extension,nhits,
			support1,support2,canonicalp,donor1,donor2,acceptor1,acceptor2);
    } else if (truestrand == '-') {
      Site_add_at_donor(chr_rev_table,chr,donorpos,acceptorpos,extension,nhits,
			support1,support2,canonicalp,donor1,donor2,acceptor1,acceptor2);
    } else if (truestrand == ' ') {
      if (need_canonical_p == false) {
	fprintf(stderr,"Adding splice at %s:%u..%u to both strands\n",chr,firstpos,secondpos);
	/* Add to both sides */
	Site_add_at_donor(chr_fwd_table,chr,/*donorpos*/firstpos,/*acceptorpos*/secondpos,
			  extension,nhits,support1,support2,canonicalp,donor1,donor2,acceptor1,acceptor2);
	Site_add_at_donor(chr_rev_table,chr,/*donorpos*/secondpos,/*acceptorpos*/firstpos,
			  extension,nhits,support1,support2,canonicalp,donor1,donor2,acceptor1,acceptor2);
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

static void
parse_splice (int nhits, char *acc, char *chr, Genomicpos_T chrpos_low,
	      Intlist_T types, Uintlist_T npositions,
	      char *auxinfo, Table_T chr_fwd_table, Table_T chr_rev_table,
	      Genome_T genome, IIT_T chromosome_iit, bool pairedp, bool concordantp) {
  Genomicpos_T firstpos, secondpos, chrpos;
  int type;
  int support1, support2;
  Intlist_T p;
  Uintlist_T q;

#if 0
  /* Doesn't hold for hard clipping */
  if (cigar_readlength != readlength) {
    fprintf(stderr,"Cigar readlength = %d, but read has length %d\n",cigar_readlength,readlength);
    exit(9);
  }
#endif

#if 0
  validlength = Samread_get_query_coordinates(&query5,&query3,types,npositions,readlength,cigar);
  debug(printf("validlength %d, min_readlength %d\n",validlength,min_readlength));
  if (validlength < min_readlength) {
    FREE(chr);
    FREE(acc);
    FREE(cigar);
    FREE(read);
    FREE(quality_string);
    return;
  }
#endif

  /* Get splice coordinates */
  chrpos = chrpos_low;
  support1 = 0;
  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos += Uintlist_head(q);
      support1 += Uintlist_head(q);

    } else if (type == 'N') {
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;
      support2 = get_support2(/*types*/Intlist_next(p),/*npositions*/Uintlist_next(q));

      add_splice(firstpos,secondpos,acc,chr,nhits,support1,support2,
		 chr_fwd_table,chr_rev_table,auxinfo,
		 genome,chromosome_iit,pairedp,concordantp);

      support1 = 0;

    } else if (type == 'I') {
      /* Do nothing */
      support1 += Uintlist_head(q);

    } else if (type == 'D') {
      /* CHECK */
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;
      if (secondpos - firstpos > 50) {
	support2 = get_support2(/*types*/Intlist_next(p),/*npositions*/Uintlist_next(q));

	add_splice(firstpos,secondpos,acc,chr,nhits,support1,support2,
		   chr_fwd_table,chr_rev_table,auxinfo,
		   genome,chromosome_iit,pairedp,concordantp);

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

#endif

#ifdef SAM_INPUT

/* SAM_INPUT */
static void
process_lines (Table_T chr_fwd_table, Table_T chr_rev_table, List_T lines, int nhits,
	       bool pairedp, bool concordantp, Genome_T genome, IIT_T chromosome_iit) {
  List_T ptr;
  char *line;
  Genomicpos_T chrpos;
  unsigned int flag;
  char *chr, *acc, *cigar, *read, *quality_string, *auxinfo;
  int readlength, cigar_readlength;
  int mapq;
  Intlist_T types;
  Uintlist_T npositions;
  
  debug0(printf("Entering process_lines with nhits %d\n",nhits));

  for (ptr = lines; ptr != NULL; ptr = List_next(ptr)) {
    line = (char *) List_head(ptr);
    debug(printf("process_lines, line is %s\n",line));
    
    auxinfo = Samread_parse_line(&acc,&flag,&mapq,&chr,&chrpos,&cigar,&readlength,&read,&quality_string,line);
    types = Samread_parse_cigar(&npositions,&cigar_readlength,cigar);

    parse_splice(nhits,acc,chr,chrpos,types,npositions,auxinfo,
		 chr_fwd_table,chr_rev_table,genome,chromosome_iit,pairedp,concordantp);

    Intlist_free(&types);
    Uintlist_free(&npositions);

    FREE(chr);
    FREE(acc);
    FREE(cigar);
    FREE(read);
    FREE(quality_string);
  }
	
  return;
}

/* Modifies global tally variables */
static void
parse_sam_input (Table_T chr_fwd_table, Table_T chr_rev_table, Genome_T genome,
		 IIT_T chromosome_iit) {
  List_T lines_first = NULL, lines_second = NULL;
  char line[1024000], *copy;
  int nhits_first = 0, nhits_second = 0;
  char *lastacc, *acc;
  unsigned int flag;
  bool pairedp, concordantp;

  lastacc = (char *) CALLOC(1,sizeof(char));
  lastacc[0] = '\0';

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '@') {
      /* Skip */
    } else {
      acc = Samread_get_acc(&flag,line);
      if (strcmp(acc,lastacc)) {
	if (lastacc[0] != '\0') {
	  if ((need_concordant_p == false || concordantp == true) &&
	      (uniquep == false || nhits_first == 1)) {
	    process_lines(chr_fwd_table,chr_rev_table,lines_first,nhits_first,
			  pairedp,concordantp,genome,chromosome_iit);
	  }
	  if ((need_concordant_p == false || concordantp == true) &&
	      (uniquep == false || nhits_second == 1)) {
	    process_lines(chr_fwd_table,chr_rev_table,lines_second,nhits_second,
			  pairedp,concordantp,genome,chromosome_iit);
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
		  flag,flag & PAIRED_READ);
	  abort();
	}
      }
    }
  }

  if (lastacc[0] != '\0') {
    if ((need_concordant_p == false || concordantp == true) &&
      (uniquep == false || nhits_first == 1)) {
      process_lines(chr_fwd_table,chr_rev_table,lines_first,nhits_first,
		    pairedp,concordantp,genome,chromosome_iit);
    }
    if ((need_concordant_p == false || concordantp == true) &&
	(uniquep == false || nhits_second == 1)) {
      process_lines(chr_fwd_table,chr_rev_table,lines_second,nhits_second,
		    pairedp,concordantp,genome,chromosome_iit);
    }
    lines_first = lines_gc(&lines_first);
    lines_second = lines_gc(&lines_second);
  }
  FREE(lastacc);

  return;
}

#elif defined(BAM_INPUT)

static List_T
bamlines_gc (List_T *lines) {
  Bamline_T bamline;
  void *item;
  
  while (*lines != NULL) {
    *lines = List_pop(*lines,&item);
    bamline = (Bamline_T) item;
    Bamline_free(&bamline);
  }
  return NULL;
}


/* BAM_INPUT */
static void
process_lines (Table_T chr_fwd_table, Table_T chr_rev_table, List_T lines, int nhits,
	       bool pairedp, bool concordantp, Genome_T genome, IIT_T chromosome_iit) {
  List_T ptr;
  Bamline_T bamline;
  
  for (ptr = lines; ptr != NULL; ptr = List_next(ptr)) {
    bamline = (Bamline_T) List_head(ptr);
    parse_splice(nhits,Bamline_acc(bamline),Bamline_chr(bamline),Bamline_chrpos_low(bamline),
		 Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),/*auxinfo*/NULL,
		 chr_fwd_table,chr_rev_table,genome,chromosome_iit,pairedp,concordantp);
  }
	
  return;
}


/* Modifies global tally variables */
static void
parse_bam_input (Bamreader_T bamreader, Table_T chr_fwd_table, Table_T chr_rev_table, Genome_T genome,
		 IIT_T chromosome_iit) {
  Bamline_T bamline;
  List_T bamlines_first = NULL, bamlines_second = NULL;
  int nhits_first = 0, nhits_second = 0;
  char *lastacc, *acc;
  unsigned int flag;
  bool pairedp, concordantp;

  lastacc = (char *) CALLOC(1,sizeof(char));
  lastacc[0] = '\0';

  while ((bamline = Bamread_next_bamline(bamreader)) != NULL) {
    acc = Bamline_acc(bamline);
    if (strcmp(acc,lastacc)) {
      if (lastacc[0] != '\0') {
	if ((need_concordant_p == false || concordantp == true) &&
	    (uniquep == false || nhits_first == 1)) {
	  process_lines(chr_fwd_table,chr_rev_table,bamlines_first,nhits_first,
			pairedp,concordantp,genome,chromosome_iit);
	}
	if ((need_concordant_p == false || concordantp == true) &&
	    (uniquep == false || nhits_second == 1)) {
	  process_lines(chr_fwd_table,chr_rev_table,bamlines_second,nhits_second,
			pairedp,concordantp,genome,chromosome_iit);
	}
	bamlines_first = bamlines_gc(&bamlines_first);
	bamlines_second = bamlines_gc(&bamlines_second);
	nhits_first = 0;
	nhits_second = 0;
      }
      FREE(lastacc);
      lastacc = (char *) CALLOC(strlen(acc)+1,sizeof(char));
      strcpy(lastacc,acc);
    }

    flag = Bamline_flag(bamline);
    if (flag & QUERY_UNMAPPED) {
      Bamline_free(&bamline);
    } else {
      if (!(flag & PAIRED_READ)) {
	pairedp = false;
	concordantp = false;
	bamlines_first = List_push(bamlines_first,(void *) bamline);
	nhits_first++;
      } else if (flag & FIRST_READ_P) {
	pairedp = true;
	concordantp = (flag & PAIRED_MAPPING) ? true : false;
	bamlines_first = List_push(bamlines_first,(void *) bamline);
	nhits_first++;
      } else if (flag & SECOND_READ_P) {
	pairedp = true;
	concordantp = (flag & PAIRED_MAPPING) ? true : false;
	bamlines_second = List_push(bamlines_second,(void *) bamline);
	nhits_second++;
      } else {
	pairedp = true;
	fprintf(stderr,"Flag %u is paired (%u), but contains neither first_read nor second_read flag\n",
		flag,flag & PAIRED_READ);
	abort();
      }
    }
  }

  if (lastacc[0] != '\0') {
    if ((need_concordant_p == false || concordantp == true) &&
      (uniquep == false || nhits_first == 1)) {
      process_lines(chr_fwd_table,chr_rev_table,bamlines_first,nhits_first,
		    pairedp,concordantp,genome,chromosome_iit);
    }
    if ((need_concordant_p == false || concordantp == true) &&
	(uniquep == false || nhits_second == 1)) {
      process_lines(chr_fwd_table,chr_rev_table,bamlines_second,nhits_second,
		    pairedp,concordantp,genome,chromosome_iit);
    }
    bamlines_first = bamlines_gc(&bamlines_first);
    bamlines_second = bamlines_gc(&bamlines_second);
  }

  FREE(lastacc);

  return;
}

#else


/* GSNAP_INPUT */
static void
process_lines (Table_T chr_fwd_table, Table_T chr_rev_table, char *header, List_T lines, int nhits) {
  List_T ptr;
  char *line, *prevline, end;
  int support1, support2, nmismatches1, nmismatches2, query5, query3;
  char strand1, strand2, s, *chr1, *chr2, firstchar1, secondchar1, firstchar2, secondchar2;
  Genomicpos_T donorpos, acceptorpos, donorpos1, acceptorpos1, donorpos2, acceptorpos2, firstpos, secondpos;
  bool sensep1, sensep2;
  Extension_T extension;
  
  debug(printf("Entering process_lines with nhits %d\n",nhits));

  if (header == NULL) {
    fprintf(stderr,"Header is NULL\n");
    abort();
  } else {
    end = header[0];
    if (end != '>' && end != '<') {
      fprintf(stderr,"Header does not begin with '<' or '>'\n");
      abort();
    }
  }

  prevline = (char *) NULL;
  for (ptr = lines; ptr != NULL; ptr = List_next(ptr)) {
    line = (char *) List_head(ptr);
    debug(printf("process_lines, line is %s\n",line));
    if (line[0] == ',') {
      if (prevline == NULL) {
	fprintf(stderr,"Got a comma line without a prevline\n");
	abort();
      }
      chr1 = chr2 = (char *) NULL;
      Gsnapread_parse_line(prevline,&query5,&query3,&firstchar1,&secondchar1,
			   &support1,&nmismatches1,&s,&chr1,&firstpos,&secondpos,
			   &donorpos1,&acceptorpos1,&strand1,&sensep1);
      Gsnapread_parse_line(line,&query5,&query3,&firstchar2,&secondchar2,
			   &support2,&nmismatches2,&s,&chr2,&firstpos,&secondpos,
			   &donorpos2,&acceptorpos2,&strand2,&sensep2);
      debug(printf("support1 %d, nmismatches1 %d\n",support1,nmismatches1));
      debug(printf("support2 %d, nmismatches2 %d\n",support2,nmismatches2));
      debug(printf("%c%c %u %u  %c%c %u %u\n",
		   firstchar1,secondchar1,donorpos1,acceptorpos1,
		   firstchar2,secondchar2,donorpos2,acceptorpos2));

      donorpos = acceptorpos = 0U;
      if (secondchar1 == '5' && firstchar2 == '3') {
	donorpos = donorpos1;
	acceptorpos = acceptorpos2;
      } else if (secondchar1 == '3' && firstchar2 == '5') {
	donorpos = donorpos2;
	acceptorpos = acceptorpos1;
      }

      if (donorpos > 0U || acceptorpos > 0U) {
	if (strcmp(chr1,chr2)) {
	  /* Translocation */
	  debug(printf("Translocation\n"));
	} else if (strand1 != strand2) {
	  /* Inversion */
	  debug(printf("Inversion\n"));
	} else if (strand1 == '+' && acceptorpos < donorpos) {
	  /* Scramble */
	  debug(printf("Scramble\n"));
	} else if (strand1 == '+' && acceptorpos > donorpos + shortsplicedist) {
	  /* Distance too far -- broken */
	  debug(printf("Too far\n"));
	} else if (strand1 == '-' && donorpos < acceptorpos) {
	  /* Scramble */
	  debug(printf("Scramble\n"));
	} else if (strand1 == '-' && donorpos > acceptorpos + shortsplicedist) {
	  /* Distance too far -- broken */
	  debug(printf("Too far\n"));
	} else {
	  debug(printf("Splice at %c%s:%u..%u\n",strand1,chr1,donorpos,acceptorpos));
	  if (sensep1 != sensep2) {
	    fprintf(stderr,"sense is inconsistent between\n");
	    fprintf(stderr,"%s",prevline);
	    fprintf(stderr,"%s",line);
	    abort();
	  } else {
	    if (Gsnapread_concordantp(header) == false) {
	      extension = NONCONCORDANT;
	    } else if (end == '>' && sensep1 == true) {
	      extension = ACCEPTOR_SIDE;
	    } else if (end == '>' && sensep1 == false) {
	      extension = DONOR_SIDE;
	    } else if (end == '<' && sensep1 == true) {
	      extension = DONOR_SIDE;
	    } else if (end == '<' && sensep1 == false) {
	      extension = ACCEPTOR_SIDE;
	    }
	  }
	  
	  /* Note: nmismatches is not supported using sam_splices */
	  if (strand1 == '+') {
	    Site_add_at_donor(chr_fwd_table,chr1,donorpos,acceptorpos,extension,nhits,
			      support1 /*- nmismatches1*/,support2 /*- nmismatches2*/,
			      /*concordantp*/true,/*donor1*/'X',/*donor2*/'X',/*acceptor1*/'X',/*acceptor2*/'X');
	  } else if (strand1 == '-') {
	    Site_add_at_donor(chr_rev_table,chr1,donorpos,acceptorpos,extension,nhits,
			      support1 /*- nmismatches1*/,support2 /*- nmismatches2*/,
			      /*concordantp*/true,/*donor1*/'X',/*donor2*/'X',/*acceptor1*/'X',/*acceptor2*/'X');
	  }
	}
      }

      FREE(chr1);
      FREE(chr2);
    }

    prevline = line;
  }
	
  return;
}

/* Modifies global tally variables */
static void
parse_gsnap_input (Table_T chr_fwd_table, Table_T chr_rev_table) {
  List_T lines;
  char line[1024000], *copy, *header = NULL;
  int nhits = 0;
  bool concordantp;

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '>' || line[0] == '<') {
      header = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header,line);
      concordantp = Gsnapread_concordantp(header);
      lines = (List_T) NULL;
      nhits = 0;

    } else if (line[0] == '\n') {
      if (header != NULL &&
	  (need_concordant_p == false || concordantp == true) &&
	  (uniquep == false || nhits == 1)) {
	lines = List_reverse(lines);
	process_lines(chr_fwd_table,chr_rev_table,header,lines,nhits);
      }
      lines = lines_gc(&lines);
      FREE(header);

    } else {
      if (line[0] != ',') {
	nhits++;
      }
      copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(copy,line);
      debug(printf("Pushing %s",copy));
      lines = List_push(lines,(void *) copy);
    }
  }

  if (header != NULL &&
      (need_concordant_p == false || concordantp == true) &&
      (uniquep == false || nhits == 1)) {
    lines = List_reverse(lines);
    process_lines(chr_fwd_table,chr_rev_table,header,lines,nhits);
  }
  lines = lines_gc(&lines);
  FREE(header);

  return;
}

#endif



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot, *mapdir = NULL, *iitfile;
  Table_T chr_fwd_table, chr_rev_table;
  Uinttable_T sitetable;
  Chrom_T *keys, chrom;
  int n, i;

#ifdef SAM_INPUT
  IIT_T chromosome_iit;
  Genome_T genome;
#elif defined(BAM_INPUT)
  IIT_T chromosome_iit;
  Genome_T genome;
  Bamreader_T bamreader;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"D:d:s:p:C:U:N:w:n:V?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 's': splicesites_file = optarg; break;

    case 'p': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'U':
      switch (atoi(optarg)) {
      case 0: uniquep = false; break;
      case 1: uniquep = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'N':
      switch (atoi(optarg)) {
      case 0: need_canonical_p = false; fprintf(stderr,"Allowing non-canonical splices\n"); break;
      case 1: need_canonical_p = true; fprintf(stderr,"Allowing only canonical splices\n"); break;
      default: fprintf(stderr,"Canonical mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'n': mincount = atoi(optarg); break;
    case 'w': minsupport = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      

  if (splicesites_file != NULL) {
    if (dbroot == NULL) {
      fprintf(stderr,"If you specify the -s flag, you need to specify the -d flag\n");
      print_program_usage();
      exit(9);
    } else {
      genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    }

    mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(splicesites_file)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,splicesites_file);
    fprintf(stderr,"Reading splicesite file %s...",splicesites_file);
    if ((splicesites_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Splicesite file %s.iit not found in %s.  Available files:\n",splicesites_file,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",splicesites_file);
      fprintf(stderr,"using the -D flag to gsnap.\n");
      exit(9);
    }
    if ((donor_typeint = IIT_typeint(splicesites_iit,"donor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'donor'\n",splicesites_file);
      exit(9);
    }
    if ((acceptor_typeint = IIT_typeint(splicesites_iit,"acceptor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'acceptor'\n",splicesites_file);
      exit(9);
    }

    fprintf(stderr,"done\n");
    FREE(iitfile);
    FREE(mapdir);
  }

#if defined(SAM_INPUT) || defined(BAM_INPUT)
  if (genomesubdir == NULL) {
    if (dbroot == NULL) {
      fprintf(stderr,"If you specify the -s flag, you need to specify the -d flag\n");
      print_program_usage();
      exit(9);
    } else {
      genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    }
  }

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);

  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);
#endif
  
  chr_fwd_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);
  chr_rev_table = Table_new(100,Chrom_compare_table,Chrom_hash_table);

#ifdef SAM_INPUT
  parse_sam_input(chr_fwd_table,chr_rev_table,genome,chromosome_iit);
#elif defined(BAM_INPUT)
  bamreader = Bamread_new(argv[0]);
  parse_bam_input(bamreader,chr_fwd_table,chr_rev_table,genome,chromosome_iit);
  Bamread_free(&bamreader);
#else
  parse_gsnap_input(chr_fwd_table,chr_rev_table);
#endif


  /* Forward table */
  n = Table_length(chr_fwd_table);
  keys = (Chrom_T *) Table_keys(chr_fwd_table,NULL);
  qsort(keys,n,sizeof(Chrom_T),Chrom_compare_chrom);

  fprintf(stderr,"Saw splices in %d fwd chromosomes:",n);
  for (i = 0; i < n; i++) {
    chrom = keys[i];
    fprintf(stderr," %s",Chrom_string(chrom));
    sitetable = Table_get(chr_fwd_table,(void *) chrom);
    print_splices(sitetable,Chrom_string(chrom),splicesites_iit,donor_typeint,acceptor_typeint);
    Uinttable_free(&sitetable);
  }
  fprintf(stderr,"\n");
  for (i = 0; i < n; i++) {
    Chrom_free(&(keys[i]));
  }
  FREE(keys);
  Table_free(&chr_fwd_table);


  /* Reverse table */
  n = Table_length(chr_rev_table);
  keys = (Chrom_T *) Table_keys(chr_rev_table,NULL);
  qsort(keys,n,sizeof(Chrom_T),Chrom_compare_chrom);

  fprintf(stderr,"Saw splices in %d rev chromosomes:",n);
  for (i = 0; i < n; i++) {
    chrom = keys[i];
    fprintf(stderr," %s",Chrom_string(chrom));
    sitetable = Table_get(chr_rev_table,(void *) chrom);
    print_splices(sitetable,Chrom_string(chrom),splicesites_iit,donor_typeint,acceptor_typeint);
    Uinttable_free(&sitetable);
  }
  fprintf(stderr,"\n");
  for (i = 0; i < n; i++) {
    Chrom_free(&(keys[i]));
  }
  FREE(keys);
  Table_free(&chr_rev_table);


#if defined(SAM_INPUT) || defined(BAM_INPUT)
  Genome_free(&genome);
  IIT_free(&chromosome_iit);
#endif

  if (splicesites_iit != NULL) {
    IIT_free(&splicesites_iit);
  }

  if (genomesubdir != NULL) {
    FREE(fileroot);
    FREE(dbversion);
    FREE(genomesubdir);
  }

  if (dbroot != NULL) {
    FREE(dbroot);
  }

  return 0;
}

