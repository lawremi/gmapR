static char rcsid[] = "$Id: singlemap.c 143437 2014-08-05 21:45:04Z twu $";
/* Note: Handles only paired-end data */
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

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "tableuint.h"
#include "tableint.h"
#include "uinttable.h"

#ifdef BAM_INPUT
#include "bamread.h"
#include "multimap.h"
#endif

#include "samflags.h"
#include "samread.h"
#include "genome.h"
#include "complement.h"

#include "bamtally.h"
#include "splice.h"
#include "gstruct.h"
#include "cappaths.h"		/* For Cappaths_setup */

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"


#define TEN_MILLION 10000000


/************************************************************************
 *   Global variables 
 ************************************************************************/

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static int filter_maximum_nhits = 10;

/* For tally */
static int alloclength = 200000;
static int blocksize = 1000;

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;


/* Output */
static bool show_excluded_p = false;


static struct option long_options[] = {
  /* Input options */
#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif

  {"read-group", required_argument, 0, 0},   /* desired_read_group */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */

  {"excluded", no_argument, 0, 'E'}, /* show_excluded_p */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_MULTIMAP\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
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
Usage: bam_multimap [-d <genome>] [OPTIONS...] <bam file>\n\
\n\
Requirements for initial tally\n\
  --read-group=STRING            Require alignments to have this read group in the RG field of\n\
                                   the BAM line\n\
  -q, --mapq=INT                 Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
\n\
Output options\n\
  -E, --excluded                 Show excluded lines, prefixed with '#'\n\
\n\
");

  return;
}



static long int
get_total_tally (IIT_T tally_iit, char *chr, Genomicpos_T coordstart, Genomicpos_T coordend) {
  long int subtotal = 0;
  int *matches;
  int nmatches, i;
  Interval_T interval;
  Genomicpos_T intervalstart, intervalend, chrpos;

  int *counts;
  int k;

  matches = IIT_get(&nmatches,tally_iit,chr,coordstart,coordend,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    counts = (int *) IIT_data(tally_iit,matches[i]);

    interval = IIT_interval(tally_iit,matches[i]);
    intervalstart = Interval_low(interval);
    intervalend = Interval_high(interval);

    if (coordstart < intervalstart) {
      k = 0;
      chrpos = intervalstart;
    } else {
      k = coordstart - intervalstart;
      chrpos = coordstart;
    }
    
    while (chrpos <= intervalend && chrpos <= coordend) {
      subtotal += counts[k];
      k++;
      chrpos++;
    }
  }

  FREE(matches);

  return subtotal;
}


long int
cigar_tally (IIT_T tally_iit, char *chr, Genomicpos_T chrpos_low,
	     Intlist_T cigar_types, Uintlist_T cigar_npositions) {
  long int total = 0U;
  Intlist_T p;
  Uintlist_T q;
  Genomicpos_T length;
  int type;
  Genomicpos_T chrpos;

  
  chrpos = chrpos_low;
  for (p = cigar_types, q = cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      length = Uintlist_head(q);
      total += get_total_tally(tally_iit,chr,/*coordstart*/chrpos,/*coordend*/chrpos + length - 1U);
      chrpos += length;

    } else if (type == 'N') {
      chrpos += Uintlist_head(q);

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      chrpos += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
  }

  return total;
}


/* Taken from Gstruct_bam_input */
static void
resolve_multimappers (Tableuint_T *resolve_low_table, Tableuint_T *resolve_high_table, Bamreader_T bamreader,
		      IIT_T tally_iit, IIT_T chromosome_iit) {
  Tableint_T bestcount_table, bestcount2_table;

  Bamline_T bamline, bamline_low;
  int nhits;
  Table_T bamstore_chrtable;
  long int total1, total2;
  char *chr, *acc, *key1, *key2;

  int index;
  Genomicpos_T chroffset;
  Interval_T interval;

  Uinttable_T bamstore_table;
  Chrom_T *chroms, chrom;
  int n, i;


  *resolve_low_table = Tableuint_new(TEN_MILLION,Table_string_compare,Table_string_hash);
  *resolve_high_table = Tableuint_new(TEN_MILLION,Table_string_compare,Table_string_hash);

  bestcount_table = Tableint_new(TEN_MILLION,Table_string_compare,Table_string_hash);
  bestcount2_table = Tableint_new(TEN_MILLION,Table_string_compare,Table_string_hash);
  bamstore_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  while ((bamline = Bamread_next_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,ignore_duplicates_p,
					 need_concordant_p)) != NULL) {
    if ((chr = Bamline_chr(bamline)) == NULL) {
      /* Skip */
      Bamline_free(&bamline);

    } else if ((nhits = Bamline_nhits(bamline)) == 1) {
      /* Unique, so skip */
      Bamline_free(&bamline);

    } else if (nhits > filter_maximum_nhits) {
      /* To be filtered, so skip */
      Bamline_free(&bamline);

    } else if (Bamline_concordantp(bamline) == false) {
      /* Handle now */
      acc = Bamline_acc(bamline);
      if (Bamline_firstend_p(bamline) == true) {
	total1 = cigar_tally(tally_iit,chr,Bamline_chrpos_low(bamline),
			     Bamline_cigar_types(bamline),
			     Bamline_cigar_npositions(bamline));
	if (total1 > Tableint_get(bestcount_table,(void *) acc)) {
	  index = IIT_find_linear(chromosome_iit,chr);
	  interval = IIT_interval(chromosome_iit,index);
	  chroffset = Interval_low(interval);
	  if (Tableuint_get(*resolve_low_table,(void *) acc) > 0) {
	    key1 = acc;
	  } else {
	    key1 = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key1,acc);
	  }
	  Tableuint_put(*resolve_low_table,(void *) key1,chroffset+Bamline_chrpos_low(bamline));
	  Tableint_put(bestcount_table,(void *) key1,total1);
	}
      } else {
	total2 = cigar_tally(tally_iit,chr,Bamline_chrpos_low(bamline),
			     Bamline_cigar_types(bamline),
			     Bamline_cigar_npositions(bamline));
	if (total2 > Tableint_get(bestcount2_table,(void *) acc)) {
	  index = IIT_find_linear(chromosome_iit,chr);
	  interval = IIT_interval(chromosome_iit,index);
	  chroffset = Interval_low(interval);
	  if (Tableuint_get(*resolve_high_table,(void *) acc) > 0) {
	    key2 = acc;
	  } else {
	    key2 = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key2,acc);
	  }
	  Tableuint_put(*resolve_high_table,(void *) key2,chroffset+Bamline_chrpos_low(bamline));
	  Tableint_put(bestcount2_table,(void *) key2,total2);
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
	total1 = cigar_tally(tally_iit,chr,Bamline_chrpos_low(bamline_low),
			     Bamline_cigar_types(bamline_low),
			     Bamline_cigar_npositions(bamline_low));
	total2 = cigar_tally(tally_iit,chr,Bamline_chrpos_low(bamline),
			     Bamline_cigar_types(bamline),
			     Bamline_cigar_npositions(bamline));
	if (total1 + total2 > Tableint_get(bestcount_table,(void *) acc)) {
#if 0
	  printf("%s %s:%u %s:%u %ld %ld\n",
		 acc,chr,Bamline_chrpos_low(bamline_low),chr,Bamline_chrpos_low(bamline),
		 total1,total2);
#endif

	  index = IIT_find_linear(chromosome_iit,chr);
	  interval = IIT_interval(chromosome_iit,index);
	  chroffset = Interval_low(interval);
	  if (Tableuint_get(*resolve_low_table,(void *) acc) > 0) {
	    key1 = key2 = acc;
	  } else {
	    key1 = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key1,acc);
	    key2 = (char *) CALLOC(strlen(acc) + 1,sizeof(char));
	    strcpy(key2,acc);
	  }
	  Tableuint_put(*resolve_low_table,(void *) key1,chroffset+Bamline_chrpos_low(bamline_low));
	  Tableuint_put(*resolve_high_table,(void *) key2,chroffset+Bamline_chrpos_low(bamline));
	  Tableint_put(bestcount_table,(void *) key1,total1 + total2);
	}
      
	Bamline_free(&bamline_low);
      }

      Bamline_free(&bamline);
    }
  }


  Tableint_free(&bestcount2_table);
  Tableint_free(&bestcount_table);

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
print_resolved (Bamreader_T bamreader, Tableuint_T resolve_low_table, Tableuint_T resolve_high_table,
		IIT_T chromosome_iit) {
  Genomicpos_T genomicpos_low, genomicpos_high, genomicpos;
  char *chr, *acc;
  unsigned int flag;
  Bamline_T bamline;

  int index;
  Interval_T interval;
  Genomicpos_T chroffset;

  while ((bamline = Bamread_next_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,ignore_duplicates_p,
					 need_concordant_p)) != NULL) {
    if (Bamline_nhits(bamline) > filter_maximum_nhits) {
      /* Skip */
    } else if ((chr = Bamline_chr(bamline)) == NULL) {
      /* Skip */
    } else {
      acc = Bamline_acc(bamline);
      flag = Bamline_flag(bamline);
      genomicpos_low = Tableuint_get(resolve_low_table,(void *) acc);
      genomicpos_high = Tableuint_get(resolve_high_table,(void *) acc);

      if (genomicpos_low == 0 && genomicpos_high == 0) {
	/* Was already unique.  Could also check nhits. */
	Bamline_print(stdout,bamline,flag,/*quality_score_adj*/0);
      } else {
	index = IIT_find_linear(chromosome_iit,chr);
	interval = IIT_interval(chromosome_iit,index);
	chroffset = Interval_low(interval);
	genomicpos = chroffset + Bamline_chrpos_low(bamline);
	if (genomicpos == genomicpos_low) {
	  Bamline_print(stdout,bamline,flag & SET_PRIMARY,/*quality_score_adj*/0);
	} else if (genomicpos == genomicpos_high) {
	  Bamline_print(stdout,bamline,flag & SET_PRIMARY,/*quality_score_adj*/0);
	} else if (show_excluded_p == true) {
	  printf("#");
	  Bamline_print(stdout,bamline,flag,/*quality_score_adj*/0);
	}
      }
    }

    Bamline_free(&bamline);
  }

  return;
}



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL;

  char *iitfile;
  IIT_T chromosome_iit = NULL;
  Genome_T genome;

  Bamreader_T bamreader;
  char *bamfile;

  IIT_T tally_iit;
  Tableuint_T resolve_low_table, resolve_high_table;
  char **keys;
  int n, i;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:q:n:C:U:P:EV?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;

      } else if (!strcmp(long_name,"read-group")) {
	desired_read_group = optarg;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'q': minimum_mapq = atoi(optarg); break;
    case 'n': maximum_nhits = atoi(optarg); break;
    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'U':
      switch (atoi(optarg)) {
      case 0: need_unique_p = false; break;
      case 1: need_unique_p = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'P':
      switch (atoi(optarg)) {
      case 0: need_primary_p = false; break;
      case 1: need_primary_p = true; break;
      default: fprintf(stderr,"Primary mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'E': show_excluded_p = true; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;


  if (dbroot != NULL) {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
			/*access*/USE_MMAP_ONLY);
  }

  bamfile = argv[0];
  fprintf(stderr,"Bamfile is %s\n",bamfile);

  fprintf(stderr,"Creating tally...\n");
  bamreader = Bamread_new(bamfile);
  tally_iit = Bamtally_iit(bamreader,/*chr*/NULL,/*bam_lacks_chr*/NULL,/*chrstart*/0,/*chrend*/0,
			   /*genome*/NULL,chromosome_iit,/*map_iit*/NULL,alloclength,
			   desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
			   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
			   /*min_depth*/1,/*variant_strands*/0,
                           /*ignore_query_Ns_p*/true,/*print_indels_p*/false,
			   blocksize,/*verbosep*/true,/*readlevel_p*/false,/*max_softclip*/0,
			   /*print_xs_scores_p*/false,/*print_noncovered_p*/false);
  Bamread_free(&bamreader);
  fprintf(stderr,"done\n");

  bamreader = Bamread_new(bamfile);
  resolve_multimappers(&resolve_low_table,&resolve_high_table,bamreader,
		       tally_iit,chromosome_iit);
  Bamread_free(&bamreader);

  IIT_free(&tally_iit);

  bamreader = Bamread_new(bamfile);
  Bamread_write_header(bamreader);
  print_resolved(bamreader,resolve_low_table,resolve_high_table,chromosome_iit);
  Bamread_free(&bamreader);

  if ((n = Tableuint_length(resolve_high_table)) > 0) {
    keys = (char **) Tableuint_keys(resolve_high_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
  }
  Tableuint_free(&resolve_high_table);

  if ((n = Tableuint_length(resolve_low_table)) > 0) {
    keys = (char **) Tableuint_keys(resolve_low_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
  }
  Tableuint_free(&resolve_low_table);

  Genome_free(&genome);
  IIT_free(&chromosome_iit);

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
