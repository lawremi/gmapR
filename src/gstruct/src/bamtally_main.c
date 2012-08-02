static char rcsid[] = "$Id: bamtally_main.c 68316 2012-07-06 20:55:23Z twu $";
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
#include <math.h>


#include "except.h"
#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "complement.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "genome.h"

#include "bamtally.h"

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"


static bool via_iit_p = false;


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

/* Filters */
static int min_depth = 1;
static int variant_strands = 0;
static bool genomic_diff_p = false;
static bool signed_counts_p = false;

/* Output options */
static bool verbosep = false;
static Tally_outputtype_T output_type = OUTPUT_BLOCKS;

static bool blockp = true;

static bool ignore_query_Ns_p = true;
static bool print_indels_p = false;
static bool print_cycles_p = false;
static bool print_quality_scores_p = false;
static bool print_mapq_scores_p = false;
static bool want_genotypes_p = false;

static char *chromosome = NULL;

static bool print_totals_p = false;


static bool print_wig_p = false;
static bool print_allele_counts_p = true;

static int blocksize = 1000;

static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = false;

/* static int min_mlength = 0; */

/* For Illumina, subtract 64.  For Sanger, subtract 33 */
/* For BAM, quality is already adjusted down to 0 */
static int quality_score_adj = 0;

#if 0
static int quality_score_constant = -1;
#endif

static int alloclength = 200000;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Filters */
  {"min-depth", required_argument, 0, 'n'}, /* min_depth */
  {"variants", required_argument, 0, 'X'}, /* variant_strands */
  {"diffs-only", no_argument, 0, 0},	   /* genomic_diff_p */

  /* Output options */
  {"verbose", required_argument, 0, 0},	       /* verbosep */
  {"block-format", required_argument, 0, 'B'}, /* blockp */
  {"format", required_argument, 0, 'A'},       /* output_type */
  {"signed-counts", no_argument, 0, 'S'}, /* signed_counts_p */
  {"totals", no_argument, 0, 'T'}, /* print_totals_p */
  {"indels", no_argument, 0, 0}, /* print_indels_p */
  {"cycles", no_argument, 0, 0}, /* print_cycles_p */
  {"quality-scores", no_argument, 0, 'Q'}, /* print_quality_scores_p */
  {"mapq-scores", no_argument, 0, 'M'}, /* print_mapq_scores_p */
  {"genotypes", no_argument, 0, 'G'}, /* want_genotypes_p */

  {"pairmax", required_argument, 0, 'p'}, /* alloclength */
  {"blocksize", required_argument, 0, 'b'}, /* blocksize */

  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {"via-iit", no_argument, 0, 0}, /* for debugging */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_TALLY\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}



static void
print_program_usage ();



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL;
  long int *tally_matches, *tally_mismatches;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL;
  int quality_counts_match[256], quality_counts_mismatch[256], i;

  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile;
  int index;
  IIT_T chromosome_iit;
  bool allocp;

  Genome_T genome = NULL;
  Bamreader_T bamreader;
  IIT_T bamtally_iit = NULL;

  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;


  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:q:n:C:U:P:A:B:p:X:STGQMb:",
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
      } else if (!strcmp(long_name,"verbose")) {
	verbosep = true;
      } else if (!strcmp(long_name,"via-iit")) {
	via_iit_p = true;
      } else if (!strcmp(long_name,"diffs-only")) {
	genomic_diff_p = true;
      } else if (!strcmp(long_name,"indels")) {
	print_indels_p = true;
      } else if (!strcmp(long_name,"cycles")) {
	print_cycles_p = true;
      } else if (!strcmp(long_name,"depth")) {
	min_depth = atoi(optarg);

#if 0
      } else if (!strcmp(long_name,"use-quality-const")) {
	quality_score_constant = atoi(optarg);
	if (quality_score_constant > MAX_QUALITY_SCORE) {
	  fprintf(stderr,"quality substitution score %d is > maximum %d\n",
		  quality_score_constant,MAX_QUALITY_SCORE);
	  exit(9);
	}
#endif
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

    case 'A':
      if (!strcmp(optarg,"blocks")) {
	output_type = OUTPUT_BLOCKS;
      } else if (!strcmp(optarg,"runlengths")) {
	output_type = OUTPUT_RUNLENGTHS;
      } else if (!strcmp(optarg,"tally")) {
	output_type = OUTPUT_TALLY;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
      }
      break;

    case 'B':
      if (!strcmp(optarg,"0")) {
	blockp = false;
      } else if (!strcmp(optarg,"1")) {
	blockp = true;
      } else {
	fprintf(stderr,"Argument to -B flag must be 0 or 1\n");
	exit(9);
      }
      break;

    case 'X': variant_strands = atoi(optarg); break;

    case 'p': alloclength = strtoul(optarg,NULL,10); break;

    case 'S': signed_counts_p = true; break;
    case 'T': print_totals_p = true; break;
    case 'G': want_genotypes_p = true; break;
    case 'Q': print_quality_scores_p = true; break;
    case 'M': print_mapq_scores_p = true; break;

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

    case 'b': blocksize = atoi(optarg); break;

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

  fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);


  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }


  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);

  bamreader = Bamread_new(argv[0]);
  if (argc > 1) {
    if ((Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			      &chroffset,&chrlength,argv[1],genomesubdir,fileroot)) == false) {
      fprintf(stderr,"Chromosome coordinates %s could not be found in the genome\n",argv[1]);
      exit(9);

    } else if (via_iit_p == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
      chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREE(iitfile);

      bamtally_iit = Bamtally_iit(bamreader,/*chr*/chromosome,chrstart,chrend,
				  genome,chromosome_iit,alloclength,
				  minimum_mapq,good_unique_mapq,maximum_nhits,
				  need_concordant_p,need_unique_p,need_primary_p,
				  min_depth,variant_strands,ignore_query_Ns_p,
				  print_indels_p,blocksize,verbosep);
      IIT_free(&chromosome_iit);

    } else {
      Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,chromosome,chroffset,chrstart,chrend,alloclength,
		   /*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   output_type,blockp,blocksize,
		   quality_score_adj,min_depth,variant_strands,
		   genomic_diff_p,signed_counts_p,ignore_query_Ns_p,
		   print_indels_p,print_totals_p,print_cycles_p,print_quality_scores_p,
		   print_mapq_scores_p,want_genotypes_p,verbosep);
      Bamread_unlimit_region(bamreader);
    }

  } else {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    if (via_iit_p == true) {
      bamtally_iit = Bamtally_iit(bamreader,/*chr*/NULL,/*chrstart*/0,/*chrend*/0,
				  genome,chromosome_iit,alloclength,
				  minimum_mapq,good_unique_mapq,maximum_nhits,
				  need_concordant_p,need_unique_p,need_primary_p,
				  min_depth,variant_strands,ignore_query_Ns_p,
				  print_indels_p,blocksize,verbosep);
    } else {
      for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
	chromosome = IIT_label(chromosome_iit,index,&allocp);
	chrstart = 1;
	chrend = Interval_length(IIT_interval(chromosome_iit,index));
	chroffset = Interval_low(IIT_interval(chromosome_iit,index));

	Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
	Bamtally_run(&tally_matches,&tally_mismatches,
		     &intervallist,&labellist,&datalist,
		     quality_counts_match,quality_counts_mismatch,
		     bamreader,genome,chromosome,chroffset,chrstart,chrend,alloclength,
		     /*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		     minimum_mapq,good_unique_mapq,maximum_nhits,
		     need_concordant_p,need_unique_p,need_primary_p,
		     /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		     output_type,blockp,blocksize,
		     quality_score_adj,min_depth,variant_strands,
		     genomic_diff_p,signed_counts_p,ignore_query_Ns_p,
		     print_indels_p,print_totals_p,print_cycles_p,print_quality_scores_p,
		     print_mapq_scores_p,want_genotypes_p,verbosep);
	Bamread_unlimit_region(bamreader);

	if (allocp == true) {
	  FREE(chromosome);
	}
      }
    }

    IIT_free(&chromosome_iit);
  }

  if (bamtally_iit != NULL) {
    fprintf(stderr,"IIT created successfully.  Printing not yet implemented\n");
    IIT_free(&bamtally_iit);
  }

  Bamread_free(&bamreader);

  if (genome != NULL) {
    Genome_free(&genome);
  }

  FREE(fileroot);
  FREE(genomesubdir);


#if 0
  if (quality_score_constant > 0) {
    make_quality_scores_constant(quality_score_constant);
  }
#endif


  for (i = 0; i < 256; i++) {
    if (quality_counts_match[i] > 0 || quality_counts_mismatch[i] > 0) {
      fprintf(stderr,"Quality %d: %d %d\n",i,quality_counts_match[i],quality_counts_mismatch[i]);
    }
  }

#if 0
  /* Should have already been freed during printing */
  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byshift[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byshift[pos]));
  }
  FREE(mismatches_byshift);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byquality[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byquality[pos]));
  }
  FREE(mismatches_byquality);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = matches_byshift[pos]; p != NULL; p = List_next(p)) {
      match = (Match_T) List_head(p);
      Match_free(&match);
    }
    List_free(&(matches_byshift[pos]));
  }
  FREE(matches_byshift);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = matches_byquality[pos]; p != NULL; p = List_next(p)) {
      match = (Match_T) List_head(p);
      Match_free(&match);
    }
    List_free(&(matches_byquality[pos]));
  }
  FREE(matches_byquality);
#endif

  
  fprintf(stderr,"Finished bam_tally on");
  for (i = 0; i < argc; i++) {
    fprintf(stderr," %s",argv[i]);
  }
  fprintf(stderr,"\n");

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_tally [OPTIONS...] bamfile chromosome:range, or\n\
       bam_tally [OPTIONS...] bamfile\n\
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
Input options (must include -d)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
\n\
Compute options\n\
  -q, --mapq=INT                 Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --pairmax=INT                  Expected insert length (discards alignments longer than\n\
                                   this value) [default=200000]\n\
\n\
Filtering of output (options may be combined)\n\
  --depth=INT                    Print only positions with this depth or more\n\
  -X, --variants=INT             Print only positions showing a variant allele\n\
                                   Argument of 0 means all positions (default)\n\
                                   Argument of 1 means an observation of 1 strand is required\n\
                                   Argument of 2 means an observation of 2 strands is required\n\
  --diffs-only                   Print only differences from genome (heterozygous or\n\
                                   homozygous different from reference)\n\
\n\
Output options\n\
  -A, --format=STRING            Output format: blocks, runlengths, bins\n\
  -B, --block-format=INT         Print in block format (0=no, 1=yes (default))\n\
  -b, --blocksize=INT            Block size for printing in block format [default 1000]\n\
  -S, --signed-counts            Print signed allele counts (as plus_count|minus_count)\n\
  -T, --totals                   Print total count (i.e., coverage or depth)\n\
  -G, --genotypes                Print genotype information (probability and likelihood)\n\
  --indels                       Print output for insertions and deletions\n\
  --cycles                       Include details about cycles\n\
  -Q, --quality-scores           Include details about quality scores\n\
  -M, --mapq-scores              Include details about mapping quality (MAPQ) scores\n\
  --verbose                      Print information about problematic reads to stderr\n\
\n\
");
    return;
}
