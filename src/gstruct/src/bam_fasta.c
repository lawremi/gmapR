static char rcsid[] = "$Id: bam_fasta.c 67870 2012-07-02 21:40:10Z twu $";
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
#include "complement.h"


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


bool allp = true;
bool single_end_p = false;
int terminal_minlength = -1;


static struct option long_options[] = {
  /* Input options */
#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif

  {"all", no_argument, 0, 'A'}, /* allp */
  {"single-end", no_argument, 0, '1'}, /* single_end_p */
  {"terminal-min-length", no_argument, 0, 'S'}, /* terminal_minlength */
  {"nomappers", no_argument, 0, 'N'}, /* allp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_FASTA -- Converts BAM results into FASTA format\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_fasta [options] <bam file...>\n\
\n\
Options:\n\
  -1, --single-end                Prints reads as single-end (even if BAM files have paired-end reads)\n\
  -S, --terminal-min-length=INT   Print only reads with terminal ends (cigar S length) equal to\n\
                                     or greater than this value (if not specified, ignores cigar S length)\n\
");
#if 0
    fprintf(stdout,"\
  -A, --all         Prints all query sequences (default)\n\
  -N, --nomappers   Prints only sequences that did not align\n\
\n\
");
#endif

    return;
}


/************************************************************************
 *   Input
 ************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static char *
string_copy (char *string) {
  char *copy;

  copy = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(copy,string);
  return copy;
}


static char *
queryseq_create (bool *completep, char *queryseq, Bamline_T bamline) {
  unsigned int flag;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  char *shortread, *p;
  char strand;
  int type;
  int cigar_querylength, mlength, querypos, i;

  flag = Bamline_flag(bamline);
  shortread = Bamline_read(bamline);
  if (flag & QUERY_UNMAPPED) {
    *completep = true;
    return string_copy(shortread);
  }

  types = Bamline_cigar_types(bamline);
  npositions = Bamline_cigar_npositions(bamline);
  cigar_querylength = Bamline_cigar_querylength(bamline);

  if (queryseq == NULL) {
    queryseq = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    for (i = 0; i < cigar_querylength; i++) {
      queryseq[i] = '*';
    }
  }

  if (flag & QUERY_MINUSP) {
    strand = '-';
    querypos = cigar_querylength - 1;
  } else {
    strand = '+';
    querypos = 0;
  }

  p = shortread;
  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);

    if (type == 'H') {
      querypos += ((strand == '+') ? +mlength : -mlength);

    } else if (type == 'N') {
      /* Ignore */
      
    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'D') {
      /* Ignore */

    } else if (strand == '+') {
      if (terminal_minlength >= 0 && type == 'S') {
	while (mlength-- > 0) {
	  queryseq[querypos] = tolower(*p);
	  querypos++;
	  p++;
	}
      } else {
      /* I and M */
	while (mlength-- > 0) {
	  queryseq[querypos] = toupper(*p);
	  querypos++;
	  p++;
	}
      }

    } else {
      if (terminal_minlength >= 0 && type == 'S') {
	while (mlength-- > 0) {
	  queryseq[querypos] = tolower(complCode[(int) *p]);
	  querypos--;
	  p++;
	}
      } else {
	/* I and M */
	while (mlength-- > 0) {
	  queryseq[querypos] = toupper(complCode[(int) *p]);
	  querypos--;
	  p++;
	}
      }
    }
  }
  
  *completep = true;
  for (i = 0; i < cigar_querylength; i++) {
    if (queryseq[i] == '*') {
      *completep = false;
    }
  }

  return queryseq;
}



static bool
queryseq_complete (char *queryseq) {
  unsigned int i;

  for (i = 0; i < strlen(queryseq); i++) {
    if (queryseq[i] == '*') {
      return false;
    }
  }

  return true;
}


/* Modifies queryseq */
static int
terminal_length (char *queryseq) {
  int terminal_length = 0, i;

  for (i = 0; i < strlen(queryseq); i++) {
    if (islower(queryseq[i])) {
      terminal_length += 1;
      queryseq[i] = toupper(queryseq[i]);
    }
  }
  return terminal_length;
}


static void
print_fasta (char *acc, char *queryseq1, char *queryseq2) {
  int terminal_length_1, terminal_length_2;

  if (single_end_p == true) {
    /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
    if (terminal_minlength < 0 || terminal_length(queryseq1) >= terminal_minlength) {
      printf(">%s/1\n",acc);
      printf("%s\n",queryseq1);
    }
    if (terminal_minlength < 0 || terminal_length(queryseq2) >= terminal_minlength) {
      printf(">%s/2\n",acc);
      printf("%s\n",queryseq2);
    }
  } else {
    if (terminal_minlength < 0) {
      printf(">%s\n",acc);
      printf("%s\n",queryseq1);
      printf("%s\n",queryseq2);
    } else {
      terminal_length_1 = terminal_length(queryseq1);
      terminal_length_2 = terminal_length(queryseq2);
      if (terminal_length_1 >= terminal_minlength || terminal_length_2 >= terminal_minlength) {
	printf(">%s\n",acc);
	printf("%s\n",queryseq1);
	printf("%s\n",queryseq2);
      }
    }
  }

  return;
}


static void
parse_bam_input (Bamreader_T bamreader, Table_T queryseq_table1, Table_T queryseq_table2) {
  Bamline_T bamline;
  char *acc, *copy;
  char *queryseq1, *queryseq2;
  bool found1p, found2p, complete1p, complete2p;

  while ((bamline = Bamread_next_bamline(bamreader,/*minimum_mapq*/0,/*good_unique_mapq*/0,/*maximum_nhits*/1000000,
					 /*need_unique_p*/false,/*need_primary_p*/true,/*need_concordant_p*/false)) != NULL) {
    acc = Bamline_acc(bamline);
    if (Bamline_firstend_p(bamline) == true) {
      queryseq1 = (char *) Table_get(queryseq_table1,acc);
      found1p = (queryseq1 ? true : false);
      queryseq1 = queryseq_create(&complete1p,queryseq1,bamline);
      if (complete1p == false) {
	Table_put(queryseq_table1,found1p ? acc : string_copy(acc),queryseq1);
      } else if ((queryseq2 = Table_get(queryseq_table2,acc)) == NULL) {
	Table_put(queryseq_table1,found1p ? acc : string_copy(acc),queryseq1);
      } else if (queryseq_complete(queryseq2) == false) {
	Table_put(queryseq_table1,found1p ? acc : string_copy(acc),queryseq1);
      } else {
	print_fasta(acc,queryseq1,queryseq2);
	if (found1p == true) {
	  copy = Table_remove(queryseq_table1,acc);
	  FREE(copy);
	}
	copy = Table_remove(queryseq_table2,acc);
	FREE(copy);
	FREE(queryseq2);
	FREE(queryseq1);
      }

    } else {
      queryseq2 = (char *) Table_get(queryseq_table2,acc);
      found2p = (queryseq2 ? true : false);
      queryseq2 = queryseq_create(&complete2p,queryseq2,bamline);
      if (complete2p == false) {
	Table_put(queryseq_table2,found2p ? acc : string_copy(acc),queryseq2);
      } else if ((queryseq1 = Table_get(queryseq_table1,acc)) == NULL) {
	Table_put(queryseq_table2,found2p ? acc : string_copy(acc),queryseq2);
      } else if (queryseq_complete(queryseq1) == false) {
	Table_put(queryseq_table2,found2p ? acc : string_copy(acc),queryseq2);
      } else {
	print_fasta(acc,queryseq1,queryseq2);
	if (found2p == true) {
	  copy = Table_remove(queryseq_table2,acc);
	  FREE(copy);
	}
	copy = Table_remove(queryseq_table1,acc);
	FREE(copy);
	FREE(queryseq1);
	FREE(queryseq2);
      }
    }

    Bamline_free(&bamline);
  }

  return;
}



int
main (int argc, char *argv[]) {
  Bamreader_T bamreader;
  Table_T queryseq_table1, queryseq_table2;
  int i;

  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"A1S:N",
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

    case 'A': allp = true; break;
    case '1': single_end_p = true; break;
    case 'S': terminal_minlength = atoi(optarg); break;
    case 'N': allp = false; break;
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  for (i = 0; i < argc; i++) {
    fprintf(stderr,"Processing file %s\n",argv[i]);
    queryseq_table1 = Table_new(1000000,Table_string_compare,Table_string_hash);
    queryseq_table2 = Table_new(1000000,Table_string_compare,Table_string_hash);

    bamreader = Bamread_new(argv[i]);
    parse_bam_input(bamreader,queryseq_table1,queryseq_table2);
    Bamread_free(&bamreader);

    Table_free(&queryseq_table2);
    Table_free(&queryseq_table1);
  }

  return 0;
}

