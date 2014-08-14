static char rcsid[] = "$Id: add-splice-prob.c 143396 2014-08-05 16:02:37Z twu $";
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
#include <ctype.h>		/* For tolower, toupper, islower */
#include <math.h>		/* For qsort */

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "genome.h"

#include "interval.h"
#include "iit-read.h"
#include "maxent_hr.h"

#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"ADD-SPLICE-PROB -- Adds splice probabilities to .exons file\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: cat <exons> | add-splice-prob [-D <genomedir>] -d <genome>\n\
\n\
Options:\n\
");

    return;
}


int
main (int argc, char *argv[]) {
  char *iitfile;

  char *genomesubdir = NULL, *fileroot = NULL;
  IIT_T chromosome_iit = NULL;
  int divno;
  Interval_T interval;
  Genomicpos_T chroffset;
  Genome_T genome;

  char source;
  char Buffer[10240], donor_chr[1000], acceptor_chr[1000], type[1000];
  int nchars;
  char donor_genestrand, acceptor_genestrand;
  Genomicpos_T donor_chrpos, acceptor_chrpos;
  int nspliced, nunpaired;
  double donor_prob, acceptor_prob;


  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:",
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
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'bam_fasta --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

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

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);

  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);
  Maxent_hr_setup(Genome_blocks(genome));

  /* Reads output of filter-exons */
  while (fgets(Buffer,10240,stdin) != NULL) {
    if (Buffer[0] == '#') {
      /* Skip line */

    } else if (sscanf(Buffer,"%c %d %d %s %c %s %u %c %s %u %n",
		      &source,&nspliced,&nunpaired,type,
		      &donor_genestrand,donor_chr,&donor_chrpos,
		      &acceptor_genestrand,acceptor_chr,&acceptor_chrpos,
		      &nchars) < 10) {
      fprintf(stderr,"Skipping %s",Buffer);

    } else {
      debug(printf("Read %s",Buffer));

      if ((divno = IIT_find_one(chromosome_iit,donor_chr)) < 0) {
	fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file...skipping\n",donor_chr);
	donor_prob = 0.0;

      } else {
	interval = IIT_interval(chromosome_iit,divno);
	chroffset = Interval_low(interval);

	if (donor_genestrand == '+') {
	  donor_prob = Maxent_hr_donor_prob(chroffset+donor_chrpos);
	} else if (donor_genestrand == '-') {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+donor_chrpos-1U);
	} else {
	  donor_prob = 0.0;
	}
      }

      if ((divno = IIT_find_one(chromosome_iit,acceptor_chr)) < 0) {
	fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file...skipping\n",acceptor_chr);
	acceptor_prob = 0.0;

      } else {
	interval = IIT_interval(chromosome_iit,divno);
	chroffset = Interval_low(interval);

	if (acceptor_genestrand == '+') {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+acceptor_chrpos-1U);
	} else if (acceptor_genestrand == '-') {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+acceptor_chrpos);
	} else {
	  acceptor_prob = 0.0;
	}
      }

      printf("%c\t%d\t%d\t%s\t%c\t%s\t%u\t%.2f\t%c\t%s\t%u\t%.2f\t%s",
	     source,nspliced,nunpaired,type,
	     donor_genestrand,donor_chr,donor_chrpos,donor_prob,
	     acceptor_genestrand,acceptor_chr,acceptor_chrpos,acceptor_prob,
	     &(Buffer[nchars]));
    }
  }




  return 0;
}

