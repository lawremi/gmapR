static char rcsid[] = "$Id: bam_get.c 51938 2011-11-08 01:23:41Z twu $";
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


static struct option long_options[] = {
  /* Output options */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */

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
  char *chromosome;
  Genomicpos_T chrstart, chrend;

  Bamreader_T bamreader;
  Bamline_T bamline;

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
  Bamread_write_header(bamreader);
  Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
  while ((bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,/*need_concordant_p*/true)) != NULL) {
    Bamline_print(stdout,bamline,Bamline_flag(bamline));
    Bamline_free(&bamline);
  }
  Bamread_unlimit_region(bamreader);
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
\n\
");
    return;
}
