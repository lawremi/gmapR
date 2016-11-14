static char rcsid[] = "$Id: bamtally_main.c 198589 2016-10-01 04:22:06Z twu $";
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


#define BUFFERLEN 1024

static bool via_iit_p = false;


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL; /* For coding regions */

/* Filters */
static bool print_noncovered_p = false;
static int min_depth = 1;
static int variant_strands = 0;
static double variant_pct = 0.05;
static bool genomic_diff_p = false;
static bool signed_counts_p = false;

/* Coordinates */
static bool whole_genome_p = false;
static char *user_region = NULL;

/* Compute options */
static int max_softclip = 0;

/* Conversion */
static char *bam_lacks_chr = NULL;
static int bam_lacks_chr_length = 0;

/* Output options */
static bool readlevel_p = false;
static bool verbosep = false;
static Tally_outputtype_T output_type = OUTPUT_BLOCKS;

static bool blockp = true;

static bool ignore_query_Ns_p = true;
static bool print_indels_p = false;
static bool print_cycles_p = false;
static bool print_nm_scores_p = false;
static bool print_xs_scores_p = false;
static bool want_genotypes_p = false;

static char *chromosome = NULL;

static bool print_totals_p = false;

static int blocksize = 1;  /* was 1000, but segment-based ref counts is not compatible with block printing */

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int minimum_quality_score = 0;
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

static int alloclength = 200000;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  /* Filters */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"min-quality", required_argument, 0, 0}, /* minimum_quality_score */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */

  {"noncovered", no_argument, 0, 0},	   /* print_noncovered_p */
  {"depth", required_argument, 0, 0}, /* min_depth */
  {"variants", required_argument, 0, 'X'}, /* variant_strands */
  {"variant-pct", required_argument, 0, 0}, /* variant_pct */
  {"diffs-only", no_argument, 0, 0},	   /* genomic_diff_p */

  /* Genomic region */
  {"whole-genome", no_argument, 0, 0}, /* whole_genome_p */
  {"region", required_argument, 0, 'r'}, /* user_region */

  /* Conversion */
  {"bam-lacks-chr", required_argument, 0, 0}, /* bam_lacks_chr */

  /* Compute options */
  {"include-soft-clips", required_argument, 0, 0}, /* max_softclip */

  /* Output options */
  {"readlevel", no_argument, 0, 0},	       /* readlevel_p */
  {"verbose", required_argument, 0, 0},	       /* verbosep */
  {"block-format", required_argument, 0, 'B'}, /* blockp */
  {"format", required_argument, 0, 'A'},       /* output_type */
  {"signed-counts", no_argument, 0, 'S'}, /* signed_counts_p */
  {"totals", no_argument, 0, 'T'}, /* print_totals_p */
  {"indels", no_argument, 0, 0}, /* print_indels_p */
  {"cycles", no_argument, 0, 0}, /* print_cycles_p */
  {"nm-scores", no_argument, 0, 0}, /* print_nm_scores_p */
  {"xs-scores", no_argument, 0, 'I'}, /* print_xs_scores_p */
  {"genotypes", no_argument, 0, 'G'}, /* want_genotypes_p */

  {"pairmax", required_argument, 0, 'p'}, /* alloclength */
  {"blocksize", required_argument, 0, 'b'}, /* blocksize */

  {"read-group", required_argument, 0, 0}, /* desired_read_group */

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
  char Buffer[BUFFERLEN], *p;
  long int *tally_matches, *tally_mismatches, grand_total;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL;
  int quality_counts_match[256], quality_counts_mismatch[256], i;

  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile, *mapdir = NULL;
  int index;
  IIT_T chromosome_iit, map_iit = NULL;
  char *chrptr;
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
  int len;

  while ((opt = getopt_long(argc,argv,"D:d:r:M:m:q:n:C:U:P:A:B:p:X:STGIb:",
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
      } else if (!strcmp(long_name,"via-iit")) {
	via_iit_p = true;

      } else if (!strcmp(long_name,"min-quality")) {
	minimum_quality_score = atoi(optarg);
      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;
      } else if (!strcmp(long_name,"noncovered")) {
	print_noncovered_p = true;
      } else if (!strcmp(long_name,"depth")) {
	min_depth = atoi(optarg);
      } else if (!strcmp(long_name,"diffs-only")) {
	genomic_diff_p = true;

      } else if (!strcmp(long_name,"variant-pct")) {
	variant_pct = atof(optarg);

      } else if (!strcmp(long_name,"indels")) {
	print_indels_p = true;
      } else if (!strcmp(long_name,"cycles")) {
	print_cycles_p = true;
      } else if (!strcmp(long_name,"nm-scores")) {
	print_nm_scores_p = true;

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
    case 'r': user_region = optarg; break;

    case 'M': user_mapdir = optarg; break;
    case 'm': 
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
      break;

    case 'A':
      if (!strcmp(optarg,"blocks")) {
	output_type = OUTPUT_BLOCKS;
      } else if (!strcmp(optarg,"runlengths")) {
	output_type = OUTPUT_RUNLENGTHS;
      } else if (!strcmp(optarg,"tally")) {
	output_type = OUTPUT_TALLY;
      } else if (!strcmp(optarg,"total")) {
	output_type = OUTPUT_TOTAL;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
	exit(9);
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
    case 'I': print_xs_scores_p = true; break;

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

  if (map_iitfile != NULL) {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    if ((map_iit = IIT_read(map_iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
			    /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      iitfile = (char*) CALLOC(strlen(mapdir)+strlen("/")+strlen(map_iitfile)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,map_iitfile);
      if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
			      /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
      }
      FREE(iitfile);
    }
    FREE(map_iitfile);
  }


  bamreader = Bamread_new(argv[0]);
  if (whole_genome_p == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    if (via_iit_p == true) {
      bamtally_iit = Bamtally_iit(bamreader,/*chr*/NULL,/*bam_lacks_chr*/NULL,
				  /*chrstart*/0,/*chrend*/0,
				  genome,chromosome_iit,map_iit,alloclength,
				  desired_read_group,minimum_mapq,good_unique_mapq,
				  minimum_quality_score,maximum_nhits,
				  need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
				  min_depth,variant_strands,variant_pct,ignore_query_Ns_p,
				  print_indels_p,blocksize,verbosep,readlevel_p,max_softclip,
				  print_cycles_p,print_nm_scores_p,print_xs_scores_p,print_noncovered_p);
    } else {
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

	if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == false) {
	  grand_total = 0;
	} else {
	  grand_total = Bamtally_run(&tally_matches,&tally_mismatches,
				     &intervallist,&labellist,&datalist,
				     bamreader,genome,chromosome,chroffset,chrstart,chrend,map_iit,
				     alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
				     desired_read_group,minimum_mapq,good_unique_mapq,
				     minimum_quality_score,maximum_nhits,
				     need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
				     /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
				     output_type,blockp,blocksize,
				     quality_score_adj,min_depth,variant_strands,variant_pct,
				     genomic_diff_p,signed_counts_p,ignore_query_Ns_p,
				     print_indels_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				     want_genotypes_p,verbosep,readlevel_p,
				     max_softclip,print_noncovered_p,/*bamfile*/argv[0]);
	  Bamread_unlimit_region(bamreader);
	}
	if (output_type == OUTPUT_TOTAL) {
	  printf("%ld\n",grand_total);
	}

	if (allocp == true) {
	  FREE(chromosome);
	}
      }
    }

    IIT_free(&chromosome_iit);

  } else if (user_region != NULL) {
    /* Coordinates given on command line */
    if ((Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			      &chroffset,&chrlength,user_region,genomesubdir,fileroot)) == false) {
      fprintf(stderr,"Chromosome coordinates %s could not be found in the genome\n",user_region);
      exit(9);

    } else if (via_iit_p == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
      chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREE(iitfile);

      bamtally_iit = Bamtally_iit(bamreader,/*chr*/chromosome,/*bam_lacks_chr*/NULL,
				  chrstart,chrend,genome,chromosome_iit,map_iit,alloclength,
				  desired_read_group,minimum_mapq,good_unique_mapq,
				  minimum_quality_score,maximum_nhits,
				  need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
				  min_depth,variant_strands,variant_pct,ignore_query_Ns_p,
				  print_indels_p,blocksize,verbosep,readlevel_p,
				  max_softclip,print_cycles_p,print_nm_scores_p,print_xs_scores_p,print_noncovered_p);
      IIT_free(&chromosome_iit);

    } else {

      if (bam_lacks_chr == NULL) {
	chrptr = chromosome;
      } else if (!strncmp(chromosome,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chromosome[bam_lacks_chr_length]);
      } else {
	chrptr = chromosome;
      }

      if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == false) {
	grand_total = 0;
      } else {
	grand_total = Bamtally_run(&tally_matches,&tally_mismatches,
				   &intervallist,&labellist,&datalist,
				   bamreader,genome,chromosome,chroffset,chrstart,chrend,map_iit,
				   alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
				   desired_read_group,minimum_mapq,good_unique_mapq,
				   minimum_quality_score,maximum_nhits,
				   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
				   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
				   output_type,blockp,blocksize,
				   quality_score_adj,min_depth,variant_strands,variant_pct,
				   genomic_diff_p,signed_counts_p,ignore_query_Ns_p,
				   print_indels_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				   want_genotypes_p,verbosep,readlevel_p,
				   max_softclip,print_noncovered_p,/*bamfile*/argv[0]);
	Bamread_unlimit_region(bamreader);
      }
      if (output_type == OUTPUT_TOTAL) {
	printf("%ld\n",grand_total);
      }
    }

  } else {
    fprintf(stderr,"Expecting coordinates from stdin.  If you want the whole genome instead, use the --whole-genome flag\n");
    /* Expecting coordinates from stdin */
    if (via_iit_p == true) {
      fprintf(stderr,"Combination of --via-iit and coordinates from stdin not supported\n");
      exit(9);
    }

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

	if (Bamread_limit_region(bamreader,chrptr,chrstart,chrend) == false) {
	  grand_total = 0;
	} else {
	  grand_total = Bamtally_run(&tally_matches,&tally_mismatches,
				     &intervallist,&labellist,&datalist,
				     bamreader,genome,chromosome,chroffset,chrstart,chrend,map_iit,
				     alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
				     desired_read_group,minimum_mapq,good_unique_mapq,
				     minimum_quality_score,maximum_nhits,
				     need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
				     /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
				     output_type,blockp,blocksize,
				     quality_score_adj,min_depth,variant_strands,variant_pct,
				     genomic_diff_p,signed_counts_p,ignore_query_Ns_p,
				     print_indels_p,print_totals_p,print_cycles_p,print_nm_scores_p,print_xs_scores_p,
				     want_genotypes_p,verbosep,readlevel_p,
				     max_softclip,print_noncovered_p,/*bamfile*/argv[0]);
	  Bamread_unlimit_region(bamreader);
	}
	if (output_type == OUTPUT_TOTAL) {
	  printf("%ld\n",grand_total);
	}
      }
      printf("# End\n");
    }
  }

  if (bamtally_iit != NULL) {
    fprintf(stderr,"IIT created successfully.  Printing not yet implemented\n");
    IIT_free(&bamtally_iit);
  }

  Bamread_free(&bamreader);

  if (genome != NULL) {
    Genome_free(&genome);
  }

  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }
  if (mapdir != NULL) {
    FREE(mapdir);
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
  --min-quality=INT              Require positions to have this quality score and higher\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --allow-duplicates             Allow alignments even if marked as duplicate (0x400) [default behavior]\n\
  --ignore-duplicates            Ignore alignments marked as duplicate (0x400)\n\
  --pairmax=INT                  Expected insert length (discards alignments longer than\n\
                                   this value) [default=200000]\n\
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
  --noncovered                   Print positions even if not covered by reads (somewhat equivalent to --depth=0,\n\
                                   but includes all positions in the given range\n\
                                   This option invalidates all other filtering options\n\
  --depth=INT                    Print only positions with this depth or more\n\
  -X, --variants=INT             Print only positions showing a variant allele\n\
                                   Argument of 0 means all positions (default)\n\
                                   Argument of 1 means an observation of 1 strand is required\n\
                                   Argument of 2 means an observation of 2 strands is required\n\
  --diffs-only                   Print only differences from genome (heterozygous or\n\
                                   homozygous different from reference)\n\
\n\
Output options\n\
  -A, --format=STRING            Output format: blocks, runlengths, bins, total\n\
  -B, --block-format=INT         Print in block format (0=no, 1=yes (default))\n\
  -b, --blocksize=INT            Block size for printing in block format [default 1000]\n\
  -S, --signed-counts            Print signed allele counts (as plus_count|minus_count)\n\
  -T, --totals                   Print total count (i.e., coverage or depth)\n\
  -G, --genotypes                Print genotype information (probability and likelihood)\n\
  --indels                       Print output for insertions and deletions\n\
  --cycles                       Include details about cycles\n\
  --nm-scores                    Include details about mismatches/indels (NM) scores\n\
  -I, --xs-scores                Include details about splice strand (XS) scores\n\
  --verbose                      Print information about problematic reads to stderr\n\
\n\
");
    return;
}
