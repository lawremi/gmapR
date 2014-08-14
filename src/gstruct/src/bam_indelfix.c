static char rcsid[] = "$Id: bam_indelfix.c 141219 2014-07-10 20:40:33Z twu $";
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
#include "dynprog.h"

#include "indelfix.h"

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"


#define BUFFERLEN 1024


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

/* Filters */
static int min_depth = 1;
static int variant_strands = 0;
static bool genomic_diff_p = false;
static bool signed_counts_p = false;

/* Coordinates */
static bool whole_genome_p = false;

/* Compute options */
static int max_softclip = 0;
static bool allow_multiple_p = false;

/* Dynamic programming */
static int max_rlength = 1000;
static int max_glength = 1000;


/* Conversion */
static char *bam_lacks_chr = NULL;
static int bam_lacks_chr_length = 0;

/* Output options */
static bool readlevel_p = false;
static bool verbosep = false;
static bool blockp = true;

static char *chromosome = NULL;

static int neighborhood_size = 20;	/* Establishes neighborhood for competing indels */

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;

/* static int min_mlength = 0; */

/* For Illumina, subtract 64.  For Sanger, subtract 33 */
/* For BAM, quality is already adjusted down to 0 */
static int quality_score_adj = 0;

#if 0
static int quality_score_constant = -1;
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Filters */
  {"min-depth", required_argument, 0, 'n'}, /* min_depth */
  {"variants", required_argument, 0, 'X'}, /* variant_strands */
  {"diffs-only", no_argument, 0, 0},	   /* genomic_diff_p */

  /* Genomic region */
  {"whole-genome", no_argument, 0, 0}, /* whole_genome_p */

  /* Conversion */
  {"bam-lacks-chr", required_argument, 0, 0}, /* bam_lacks_chr */

  /* Compute options */
  {"include-soft-clips", required_argument, 0, 0}, /* max_softclip */
  {"allow-multiple", no_argument, 0, 0}, /* allow_multiple_p */

  /* Output options */
  {"readlevel", no_argument, 0, 0},	       /* readlevel_p */
  {"verbose", required_argument, 0, 0},	       /* verbosep */
  {"block-format", required_argument, 0, 'B'}, /* blockp */
  {"signed-counts", no_argument, 0, 'S'}, /* signed_counts_p */
  {"totals", no_argument, 0, 'T'}, /* print_totals_p */
  {"indels", no_argument, 0, 0}, /* print_indels_p */
  {"cycles", no_argument, 0, 0}, /* print_cycles_p */
  {"quality-scores", no_argument, 0, 'Q'}, /* print_quality_scores_p */
  {"mapq-scores", no_argument, 0, 'M'}, /* print_mapq_scores_p */
  {"genotypes", no_argument, 0, 'G'}, /* want_genotypes_p */

  {"neighborhood-size", required_argument, 0, 'b'}, /* neighborhood_size */

  {"read-group", required_argument, 0, 0}, /* desired_read_group */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {"via-iit", no_argument, 0, 0}, /* for debugging */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_INDELFIX\n");
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
  char Buffer[BUFFERLEN], *p;
  int i;

  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile;
  int index;
  IIT_T chromosome_iit;
  char *chrptr;
  bool allocp;

  Genome_T genome = NULL;
  Bamreader_T bamreader;

  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;


  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:q:n:C:U:P:B:X:Sb:",
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
      } else if (!strcmp(long_name,"readlevel")) {
	readlevel_p = true;
      } else if (!strcmp(long_name,"verbose")) {
	verbosep = true;
      } else if (!strcmp(long_name,"whole-genome")) {
	whole_genome_p = true;
      } else if (!strcmp(long_name,"bam-lacks-chr")) {
	bam_lacks_chr = optarg;
	bam_lacks_chr_length = strlen(bam_lacks_chr);
      } else if (!strcmp(long_name,"include-soft-clips")) {
	max_softclip = atoi(optarg);
      } else if (!strcmp(long_name,"allow-multiple")) {
	allow_multiple_p = true;
      } else if (!strcmp(long_name,"diffs-only")) {
	genomic_diff_p = true;
      } else if (!strcmp(long_name,"depth")) {
	min_depth = atoi(optarg);
      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;
      } else if (!strcmp(long_name,"read-group")) {
	desired_read_group = optarg;

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

    case 'S': signed_counts_p = true; break;

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

    case 'b': neighborhood_size = atoi(optarg); break;

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


  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);
  Dynprog_init(/*mode*/STANDARD);


  bamreader = Bamread_new(argv[0]);
  if (whole_genome_p == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chromosome = IIT_label(chromosome_iit,index,&allocp);
      chrstart = 1;
      chrend = Interval_length(IIT_interval(chromosome_iit,index));
      chroffset = Interval_low(IIT_interval(chromosome_iit,index));

      if (bam_lacks_chr == NULL) {
	chrptr = chromosome;
      } else if (!strncmp(chromosome,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chromosome[bam_lacks_chr_length]);
      } else {
	chrptr = chromosome;
      }

      if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == true) {
	Indelfix_run(bamreader,max_rlength,max_glength,
		     genome,chromosome,chroffset,chrstart,chrend,
		     desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
		     need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		     /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		     neighborhood_size,allow_multiple_p,/*bamfile*/argv[0]);
	Bamread_unlimit_region(bamreader);
      }

      if (allocp == true) {
	FREE(chromosome);
      }
    }

    IIT_free(&chromosome_iit);

  } else if (argc > 1) {
    /* Coordinates given on command line */
    if ((Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			      &chroffset,&chrlength,argv[1],genomesubdir,fileroot)) == false) {
      fprintf(stderr,"Chromosome coordinates %s could not be found in the genome\n",argv[1]);
      exit(9);

    } else {

      if (bam_lacks_chr == NULL) {
	chrptr = chromosome;
      } else if (!strncmp(chromosome,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chromosome[bam_lacks_chr_length]);
      } else {
	chrptr = chromosome;
      }

      if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == true) {
	Indelfix_run(bamreader,max_rlength,max_glength,
		     genome,chromosome,chroffset,chrstart,chrend,
		     desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
		     need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		     /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		     neighborhood_size,allow_multiple_p,/*bamfile*/argv[0]);
	Bamread_unlimit_region(bamreader);
      }
    }

  } else {
    fprintf(stderr,"Expecting coordinates from stdin.  If you want the whole genome instead, use the --whole-genome flag\n");
    /* Expecting coordinates from stdin */
    while (fgets(Buffer,BUFFERLEN,stdin) != NULL) {
      printf("# Query: %s",Buffer);
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
      }

      /* Truncate query at first space character */
      p = Buffer;
      while (*p != '\0' && !isspace(*p)) {
	p++;
      }
      if (isspace(*p)) {
	*p = '\0';
      }

      if ((Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
				&chroffset,&chrlength,Buffer,genomesubdir,fileroot)) == false) {
	fprintf(stderr,"Chromosome coordinates %s could not be found in the genome\n",Buffer);

      } else {
	if (bam_lacks_chr == NULL) {
	  chrptr = chromosome;
	} else if (!strncmp(chromosome,bam_lacks_chr,bam_lacks_chr_length)) {
	  chrptr = &(chromosome[bam_lacks_chr_length]);
	} else {
	  chrptr = chromosome;
	}

	if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == true) {
	  Indelfix_run(bamreader,max_rlength,max_glength,
		       genome,chromosome,chroffset,chrstart,chrend,
		       desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
		       need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		       /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		       neighborhood_size,allow_multiple_p,/*bamfile*/argv[0]);
	  Bamread_unlimit_region(bamreader);
	}
      }
      printf("# End\n");
    }
  }

  Bamread_free(&bamreader);

  Dynprog_term();

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


  fprintf(stderr,"Finished bam_indelfix on");
  for (i = 0; i < argc; i++) {
    fprintf(stderr," %s",argv[i]);
  }
  fprintf(stderr,"\n");

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_indelfix [OPTIONS...] bamfile chromosome:range, or\n\
       bam_indelfix [OPTIONS...] bamfile\n\
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
If chromosome:range not given, then either entire genome is processed\n\
if --whole-genome flag is given, or the program expects to get chromosomal\n\
coordinates via stdin.\n\
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
  --allow-duplicates             Allow alignments even if marked as duplicate (0x400) [default behavior]\n\
  --ignore-duplicates            Ignore alignments marked as duplicate (0x400)\n\
  --read-group=STRING            Process only alignments that have the given read group as an RG field\n\
                                   in their BAM lines\n\
  --include-soft-clips=INT       Include soft clips of up to this length in the tally results\n\
                                   (May want to set to maximum read length) [default = 0]\n\
                                   Soft clips as part of terminal alignments (marked with PG:Z:T\n\
                                   by GSNAP) are excluded\n\
\n\
Coordinates\n\
  --whole-genome                 Compute tally over entire genome\n\
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
  -B, --block-format=INT         Print in block format (0=no, 1=yes (default))\n\
  -b, --neighborhood-size=INT    Neighborhood size for resolving competing indels [default 20]\n\
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
