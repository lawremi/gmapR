static char rcsid[] = "$Id: bam_splices.c 138416 2014-06-06 21:10:59Z twu $";
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
#include <math.h>		/* For qsort, rint */

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "chrom.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "table.h"
#include "uinttable.h"
#include "maxent_hr.h"		/* For Maxent_hr_setup */

#ifdef BAM_INPUT
#include "bamread.h"
#endif

#include "samflags.h"
#include "samread.h"
#include "genome.h"
#include "complement.h"
#include "tally.h"

#include "bamtally.h"
#include "gstruct.h"

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static IIT_T map_iit = NULL;

/* Filters */
static int min_depth = 1;

/* For splicegene */
static int max_exonlength = 2000;

/* Output options */
static bool separate_extents_p = false;


static char *chromosome = NULL;
static bool whole_genome_p = false;
static char *bam_lacks_chr = NULL;
static int bam_lacks_chr_length = 0;

static bool show_invalid_p = false;
static bool need_canonical_p = false;

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;
static bool need_concordant_p = false;

static Genomicpos_T max_pairlength = 1000000;
static Genomicpos_T shortsplicedist = 400000;

static int mincount = 2;
static int minsupport = 8;


/* static int min_mlength = 0; */

#if 0
static int quality_score_constant = -1;
#endif

static int alloclength = 400000;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Known genes */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfiles */

  /* Filters */
  {"depth", required_argument, 0, 0},	   /* min_depth */

  /* Genomic region */
  {"whole-genome", no_argument, 0, 0}, /* whole_genome_p */
  {"bam-lacks-chr", required_argument, 0, 0}, /* bam_lacks_chr */

  /* Output options */
  {"sep-tallies", no_argument, 0, 0}, /* separate_tallies_p */
  {"sep-extents", no_argument, 0, 0}, /* separate_extents_p */

  {"read-group", required_argument, 0, 0},   /* desired_read_group */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"pairmax", required_argument, 0, 'p'}, /* alloclength */
  {"show-invalid", no_argument, 0, 0},	  /* show_invalid_p */
  

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_SPLICES\n");
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
  Splice_T splice;

  List_T p;
  int quality_counts_match[256], quality_counts_mismatch[256], i;

  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile, *mapdir;
  IIT_T chromosome_iit;
  int index;
  bool allocp;

  List_T bamfiles;
  Genome_T genome = NULL;
  int readlength, insertlength;
  Gstruct_T gstruct = NULL;
  List_T splices = NULL;

  long int *fwd_extents = NULL, *rev_extents = NULL, *null_extents = NULL,
    *primary_extents = NULL, *crosshyb_extents = NULL;

  char *chr, *chrptr;
  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;
  int len;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;


  while ((opt = getopt_long(argc,argv,"D:d:q:n:C:U:P:M:m:p:",
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

      } else if (!strcmp(long_name,"whole-genome")) {
	whole_genome_p = true;
      } else if (!strcmp(long_name,"bam-lacks-chr")) {
	bam_lacks_chr = optarg;
	bam_lacks_chr_length = strlen(bam_lacks_chr);

      } else if (!strcmp(long_name,"sep-extents")) {
	separate_extents_p = true;

      } else if (!strcmp(long_name,"depth")) {
	min_depth = atoi(optarg);
      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;

      } else if (!strcmp(long_name,"read-group")) {
	desired_read_group = optarg;

      } else if (!strcmp(long_name,"show-invalid")) {
	show_invalid_p = true;
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

    case 'M': user_mapdir = optarg; break;
    case 'm':
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
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

    case 'p': alloclength = strtoul(optarg,NULL,10); break;

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

  if (map_iitfile == NULL) {
    /* fprintf(stderr,"Must specify genes iit with the -m flag\n"); */
    map_iit = (IIT_T) NULL;
  } else {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
    if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",map_iitfile);
      fprintf(stderr,"using the -M flag to gmap.\n");
      exit(9);
    }
    FREE(mapdir);
    FREE(map_iitfile);
    FREE(iitfile);
  }


  fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);

  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }

  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);
  Maxent_hr_setup(Genome_blocks(genome));

  /* bamreader = Bamread_new(argv[0]); */
  bamfiles = List_push(NULL,argv[0]);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);


  if (whole_genome_p == true) {
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

      gstruct = Gstruct_bam_input(&readlength,&insertlength,bamfiles,chrptr,chrstart,chrend,
				  /*ngoodhits_low_table*/NULL,/*ngoodhits_high_table*/NULL,
				  /*genes_iit*/NULL,shortsplicedist,max_pairlength,genome,chromosome_iit,
				  desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_unique_p,need_primary_p,
				  ignore_duplicates_p,/*trust_sam_p*/true,need_canonical_p,
				  bam_lacks_chr,bam_lacks_chr_length);
      
      if ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				  &primary_extents,&crosshyb_extents,&chrlength,
				  gstruct,max_exonlength,mincount,
				  minsupport,need_canonical_p,bam_lacks_chr,bam_lacks_chr_length)) != NULL) {
	/* Report splices */
	for (p = splices; p != NULL; p = List_next(p)) {
	  splice = (Splice_T) List_head(p);
	  Splice_print(splice,chr/*not chrptr*/,map_iit,show_invalid_p);
	}

	FREE(fwd_extents);
	FREE(rev_extents);
	FREE(null_extents);
	FREE(primary_extents);
	FREE(crosshyb_extents);
	List_free(&splices);
      }

      Gstruct_free(&gstruct);
    }

  } else {
    if (argc <= 1) {
      fprintf(stderr,"Usage: bam_splices -d <genome> [-m <genes iit, comma delimited>] <bam file> <coords>\n");
      exit(9);
    }

    Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			 &chroffset,&chrlength,argv[1],genomesubdir,fileroot);
    fprintf(stderr,"GMAP index says chrlength is %u\n",chrlength);

    if (bam_lacks_chr == NULL) {
      chrptr = chromosome;
    } else if (!strncmp(chromosome,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chromosome[bam_lacks_chr_length]);
    } else {
      chrptr = chromosome;
    }

    gstruct = Gstruct_bam_input(&readlength,&insertlength,bamfiles,chrptr,chrstart,chrend,
				/*ngoodhits_low_table*/NULL,/*ngoodhits_high_table*/NULL,
				/*genes_iit*/NULL,shortsplicedist,max_pairlength,genome,chromosome_iit,
				desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_unique_p,need_primary_p,
				ignore_duplicates_p,/*trust_sam_p*/true,need_canonical_p,
				bam_lacks_chr,bam_lacks_chr_length);
    
    if ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				&primary_extents,&crosshyb_extents,&chrlength,
				gstruct,max_exonlength,mincount,
				minsupport,need_canonical_p,bam_lacks_chr,bam_lacks_chr_length)) == NULL) {
      fprintf(stderr,"No data in given region\n");
      exit(0);
    } else {
      fprintf(stderr,"BAM file says chrlength of chr %s is %u\n",chrptr,chrlength);
    }

    /* Report splices */
    for (p = splices; p != NULL; p = List_next(p)) {
      splice = (Splice_T) List_head(p);
      Splice_print(splice,chr/*not chrptr*/,map_iit,show_invalid_p);
    }

    FREE(fwd_extents);
    FREE(rev_extents);
    FREE(null_extents);
    FREE(primary_extents);
    FREE(crosshyb_extents);
    Gstruct_free(&gstruct);
    List_free(&splices);
  }
  
  IIT_free(&chromosome_iit);

  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }

  if (genome != NULL) {
    Genome_free(&genome);
  }

  List_free(&bamfiles);

  FREE(fileroot);
  FREE(genomesubdir);

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_splices [OPTIONS...] bamfile chromosome:range, or\n\
       bam_splices [OPTIONS...] bamfile\n			      \
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
Input options (must include -d)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -M, --mapdir=STRING            Map file directory.  Program will for map file as given and also here.\n\
  -m, --map=STRING               Map file(s), comma-delimited\n\
  --whole-genome                 Process whole genome\n\
  --bam-lacks-chr=STRING\n\
\n\
Compute options\n\
  --read-group=STRING            Require alignments to have this read group in the RG field of\n\
                                   the BAM line\n\
  -q, --min-mapq=INT             Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default: 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --allow-duplicates             Allow alignments even if marked as duplicate (0x400) [default behavior]\n\
  --ignore-duplicates            Ignore alignments marked as duplicate (0x400)\n\
  --pairmax=INT                  Expected insert length (reserves memory for this amount, so\n\
                                   alignments longer than this value are discarded) [default=400000]\n\
  --show-invalid                 Include invalid splices (in opposite direction of apparent gene)\n\
\n\
Filtering of output (options may be combined)\n\
  --depth=INT                    Print only positions with this depth or more\n\
\n\
");
    return;
}
