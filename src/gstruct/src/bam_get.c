static char rcsid[] = "$Id: bam_get.c 87713 2013-03-01 18:32:34Z twu $";
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


#ifdef BAM_INPUT
#include "bamread.h"
#endif

#include "samflags.h"
#include "samread.h"

#include "parserange.h"
#include "getopt.h"


/* Filters */
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;
static int quality_score_adj = 64;


static struct option long_options[] = {
  /* Output options */
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
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_GET\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}



static void
print_program_usage ();



int
main (int argc, char *argv[]) {
  char *chromosome, *mate_chr;
  Genomicpos_T chrstart, chrend, mate_chrpos_low, mate_chrpos_high;

  Bamreader_T bamreader, bamreader_mate;
  Bamline_T bamline, bamline_mate;

  bool revcomp;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;


  while ((opt = getopt_long(argc,argv,"q:n:C:U:P:",
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

      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
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

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (argc <= 1) {
    fprintf(stderr,"Usage: bam_get <bam file> <coords>\n");
    exit(9);
  }

  Parserange_simple(&chromosome,&revcomp,&chrstart,&chrend,argv[1]);


  bamreader = Bamread_new(argv[0]);
  bamreader_mate = Bamread_new(argv[0]);

  Bamread_write_header(bamreader);
  Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
  while ((bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,ignore_duplicates_p,
					 /*need_concordant_p*/true)) != NULL) {
    Bamline_print(stdout,bamline,Bamline_flag(bamline),quality_score_adj);
    if ((mate_chr = Bamline_mate_chr(bamline)) == NULL) {
      /* No mate.  Need to know how to get an unmapped mate. */
    } else if (strcmp(mate_chr,chromosome) || Bamline_mate_chrpos_low(bamline) > chrend) {
      /* Mate will not be retrieved */
      mate_chrpos_low = Bamline_mate_chrpos_low(bamline);
      bamline_mate = Bamread_get_acc(bamreader_mate,mate_chr,mate_chrpos_low,Bamline_acc(bamline));
      if (bamline_mate != NULL) {
	Bamline_print(stdout,bamline_mate,Bamline_flag(bamline_mate),quality_score_adj);
	Bamline_free(&bamline_mate);
      }
    } else {
      mate_chrpos_low = Bamline_mate_chrpos_low(bamline);
      bamline_mate = Bamread_get_acc(bamreader_mate,mate_chr,mate_chrpos_low,Bamline_acc(bamline));
      if (bamline_mate != NULL) {
	mate_chrpos_high = Bamline_chrpos_high(bamline_mate);
	if (mate_chrpos_high < chrstart) {
	  /* Mate will not be retrieved */
	  Bamline_print(stdout,bamline_mate,Bamline_flag(bamline_mate),quality_score_adj);
	}
	Bamline_free(&bamline_mate);
      }
    }

    Bamline_free(&bamline);
  }
  Bamread_unlimit_region(bamreader);

  Bamread_free(&bamreader_mate);
  Bamread_free(&bamreader);

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_get [OPTIONS...] bamfile chromosome:range\n\
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
Filtering options\n\
  -q, --min-mapq=INT             Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --allow-duplicates             Allow alignments even if marked as duplicate (0x400) [default behavior]\n\
  --ignore-duplicates            Ignore alignments marked as duplicate (0x400)\n\
\n\
");
    return;
}
