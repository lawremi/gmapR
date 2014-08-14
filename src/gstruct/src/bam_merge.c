static char rcsid[] = "$Id: bam_merge.c 138417 2014-06-06 21:12:15Z twu $";
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
#include "list.h"
#include "getopt.h"

#include "samflags.h"
#include "bamread.h"
#include "samread.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static bool allp = true;


static struct option long_options[] = {
  /* Input options */
  {"changed-only", no_argument, 0, 0},	/* allp */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_MERGE -- Integrates revised SAM lines into a BAM file\n");
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
Usage: cat <revised sam> | bam_merge [options] <bam file...>\n\
\n\
Options:\n\
");

    return;
}


static bool strip_warning_p = false;
static bool stripp = true;
static bool clip_overlap_p = false;


/* Deletes /\D1/ or /\D2 or 3/ endings. */
static char *
strip_illumina_acc_ending (bool *copyp, char *acc) {
  char *newacc = NULL;
  char *p;
  char slash;
  int strlength;

  *copyp = false;
  if (stripp == true) {
    p = acc;
    while (*p != '\0') {
      p++;
    }

    /* Handle old style Illumina data that ends in ":0" or ":1" */
    if (p[-2] == ':') {
      *copyp = true;
      p -= 2;
    }

    /* Delete "/1" or "/2 or /3" endings */
    slash = p[-2];

    if (slash != '/' && slash != '#') {
      /* fprintf(stderr,"Do not see /1 or /2 endings in header fields %s and %s.  Will no longer look for them.\n",acc1,acc2); */
      stripp = false;

    } else if (p[-1] != '1' && p[-1] != '2' && p[-1] != '3') {
      /* fprintf(stderr,"Do not see /1 or /2 endings in header fields %s and %s.  Will no longer look for them.\n",acc1,acc2); */
      stripp = false;

    } else {
      if (strip_warning_p == false) {
	fprintf(stderr,"BAM file is improper, containing /1 or /2 endings for accessions (e.g., %s).  Will strip these endings.\n",acc);
	strip_warning_p = true;
      }
      *copyp = true;
      p -= 2;
    }
  }

  if (*copyp == false) {
    return acc;
  } else {
    strlength = (p - acc)/sizeof(char);
    newacc = (char *) CALLOC(strlength+1,sizeof(char));
    strncpy(newacc,acc,strlength);
    return newacc;
  }
}


/************************************************************************
 *   Input
 ************************************************************************/

static char *
string_copy (char *string) {
  char *copy;

  copy = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(copy,string);
  return copy;
}


static void
parse_bam_input (Bamreader_T bamreader, Table_T revision_table_0, Table_T revision_table_1, Table_T revision_table_2) {
  Bamline_T bamline;
  char *acc;
  char *revised_line, *mate_line;
  bool copyp;
  char *chr, *mate_chr, *cigar, *mate_cigar;
  Genomicpos_T chrpos, mate_chrpos;
  int insert_length;
  int readlength, mate_readlength;

  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					 /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/true,
					 /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    acc = Bamline_acc(bamline);
    acc = strip_illumina_acc_ending(&copyp,acc);
    debug(printf("Acc is %s",acc));
    if (Bamline_paired_read_p(bamline) == false) {
      debug(printf(" single-end\n"));
      if ((revised_line = (char *) Table_get(revision_table_0,acc)) == NULL) {
	if (allp == true) {
	  /* Print line as is */
	  Bamline_print(/*fp*/stdout,bamline,/*newflag*/Bamline_flag(bamline),/*quality_score_adj*/0);
	}
      } else {
	/* Replace with revised line */
	printf("%s",revised_line);
      }

    } else {

      if (Bamline_firstend_p(bamline) == true) {
	debug(printf(" paired-end 1\n"));
	revised_line = (char *) Table_get(revision_table_1,acc);
	mate_line = (char *) Table_get(revision_table_2,acc);
      } else {
	debug(printf(" paired-end 2\n"));
	revised_line = (char *) Table_get(revision_table_2,acc);
	mate_line = (char *) Table_get(revision_table_1,acc);
      }

      if (revised_line == NULL && mate_line == NULL) {
	if (allp == true) {
	  /* Print line as is */
	  Bamline_print(/*fp*/stdout,bamline,/*newflag*/Bamline_flag(bamline),/*quality_score_adj*/0);
	}

      } else if (revised_line != NULL && mate_line == NULL) {
	/* Replace with revised line.  No change to mate information. */
	if (clip_overlap_p == false) {
	  printf("%s",revised_line);
	} else {
	  /* Unfortunately, we can't get the mate_readlength or CIGAR, so we can't compute the insert length or overlap */
	  printf("%s",revised_line);
	}

      } else if (revised_line == NULL && mate_line != NULL) {
	/* Print line as is, except with change to mate information */
	mate_chr = Samread_chrinfo(&mate_chrpos,&mate_cigar,mate_line);

	chr = Bamline_chr(bamline);
	chrpos = Bamline_chrpos_low(bamline);
	cigar = Bamline_cigar_string(bamline);

	if (strcmp(chr,mate_chr)) {
	  insert_length = 0;
	} else {
	  insert_length = Samread_compute_insert_length(&readlength,&mate_readlength,cigar,chrpos,mate_cigar,mate_chrpos);
	}

	if (clip_overlap_p == false) {
	  Bamline_print_new_mate(/*fp*/stdout,bamline,mate_chr,mate_chrpos,insert_length);
	} else if (insert_length == 0) {
	  Bamline_print_new_mate(/*fp*/stdout,bamline,mate_chr,mate_chrpos,insert_length);
	} else if (insert_length >= readlength + mate_readlength) {
	  Bamline_print_new_mate(/*fp*/stdout,bamline,mate_chr,mate_chrpos,insert_length);
	} else {
	  fprintf(stderr,"Not supported yet\n");
	}

	FREE(cigar);
	FREE(mate_cigar);
	FREE(mate_chr);

      } else {
	/* Replace with revised line, and with change to mate information */
	mate_chr = Bamline_chr(bamline);
	mate_chrpos = Bamline_chrpos_low(bamline);
	mate_cigar = Bamline_cigar_string(bamline);

	chr = Samread_chrinfo(&chrpos,&cigar,revised_line);
	if (strcmp(chr,mate_chr)) {
	  insert_length = 0;
	} else {
	  insert_length = Samread_compute_insert_length(&mate_readlength,&readlength,mate_cigar,mate_chrpos,cigar,chrpos);
	}

	if (clip_overlap_p == false || insert_length == 0) {
	  Samread_print_altered_mate(/*fp*/stdout,chr,chrpos,mate_chr,mate_chrpos,insert_length,revised_line);
	} else if (insert_length == 0) {
	  Samread_print_altered_mate(/*fp*/stdout,chr,chrpos,mate_chr,mate_chrpos,insert_length,revised_line);
	} else if (insert_length >= readlength + mate_readlength) {
	  Samread_print_altered_mate(/*fp*/stdout,chr,chrpos,mate_chr,mate_chrpos,insert_length,revised_line);
	} else {
	  fprintf(stderr,"Not supported yet\n");
	}

	FREE(cigar);
	FREE(chr);
      }
    }

    if (copyp) {
      FREE(acc);
    }

    Bamline_free(&bamline);
  }

  return;
}

static void
parse_sam_input (Table_T revision_table_0, Table_T revision_table_1, Table_T revision_table_2) {
  char line[1024000];
  char *acc, *acc_orig;
  unsigned int flag;
  bool copyp;

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '@') {
      /* Skip */
    } else {
      acc_orig = Samread_get_acc(&flag,line);
      acc = strip_illumina_acc_ending(&copyp,acc_orig);
      if (copyp == true) {
	FREE(acc_orig);
      } else {
	acc = acc_orig;
      }
      debug(printf("Acc is %s",acc));

      if (flag & NOT_PRIMARY) {
	/* Skip non-primary alignments */
	FREE(acc);

      } else if (!(flag & PAIRED_READ)) {
	debug(printf(" single-end\n"));
	Table_put(revision_table_0,acc,string_copy(line));

      } else if (flag & FIRST_READ_P) {
	debug(printf(" paired-end 1\n"));
	Table_put(revision_table_1,acc,string_copy(line));

      } else if (flag & SECOND_READ_P) {
	debug(printf(" paired-end 2\n"));
	Table_put(revision_table_2,acc,string_copy(line));

      } else {
	fprintf(stderr,"Flag %u is paired (%u), but contains neither first_read nor second_read flag\n",
		flag,flag & PAIRED_READ);
	abort();
      }
    }
  }

  return;
}


static void
empty_table (Table_T revision_table) {
  char *acc, *acc_copy;
  char **accessions;
  char *revised_line;
  int n, i;

  n = Table_length(revision_table);
  accessions = (char **) Table_keys(revision_table,(void *) NULL);
  for (i = 0; i < n; i++) {
    acc = accessions[i];
    revised_line = (char *) Table_get(revision_table,acc);

    acc_copy = Table_remove(revision_table,acc);
    FREE(acc_copy);
    FREE(revised_line);
  }
  FREE(accessions);

  return;
}


int
main (int argc, char *argv[]) {
  Bamreader_T bamreader;
  Table_T revision_table_0, revision_table_1, revision_table_2;
  int i;

  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"",
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

      } else if (!strcmp(long_name,"changed-only")) {
	allp = false;

      } else if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'bam_merge --help'",long_name);
	exit(9);
      }
      break;

    case '?': fprintf(stderr,"For usage, run 'bam_merge --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  revision_table_0 = Table_new(1000000,Table_string_compare,Table_string_hash);
  revision_table_1 = Table_new(1000000,Table_string_compare,Table_string_hash);
  revision_table_2 = Table_new(1000000,Table_string_compare,Table_string_hash);
  parse_sam_input(revision_table_0,revision_table_1,revision_table_2);


  for (i = 0; i < argc; i++) {
    fprintf(stderr,"Processing file %s\n",argv[i]);
    bamreader = Bamread_new(argv[i]);
    parse_bam_input(bamreader,revision_table_0,revision_table_1,revision_table_2);
    Bamread_free(&bamreader);
  }

  empty_table(revision_table_2);
  empty_table(revision_table_1);
  empty_table(revision_table_0);

  Table_free(&revision_table_2);
  Table_free(&revision_table_1);
  Table_free(&revision_table_0);

  return 0;
}

