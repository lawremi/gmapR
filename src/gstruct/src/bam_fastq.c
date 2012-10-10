static char rcsid[] = "$Id: bam_fastq.c 73122 2012-09-01 03:51:47Z twu $";
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


static char *output_root = "out";
static bool fastap = false;
static bool allp = true;
static bool single_end_output_p = false; /* print paired-end BAM as single-end file */
static int terminal_minlength = -1;

/* For Illumina, add 64.  For Sanger, add 33 */
static int quality_shift = 33;


static struct option long_options[] = {
  /* Input options */
#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif

  {"output", required_argument, 0, 'o'},	/* output_root */
  {"fasta", no_argument, 0, 0},			/* fastap */
  {"all", no_argument, 0, 'A'}, /* allp */
  {"single-end-output", no_argument, 0, '1'}, /* single_end_output_p */
  {"terminal-min-length", no_argument, 0, 'S'}, /* terminal_minlength */
  {"nomappers", no_argument, 0, 'N'}, /* allp */
  {"quality-protocol", required_argument, 0, 0}, /* quality_shift */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_FASTQ -- Converts BAM results into FASTQ format\n");
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
Usage: bam_fastq [options] [-o <output>] <bam file...>\n\
\n\
Options:\n\
  --fasta                         Writes to stdout in FASTA format\n\
  -o, --output=STRING             Writes FASTQ to <output>_1.fq and <output>_2.fq (default=out)\n\
  -1, --single-end-output         Prints reads as single-end (even if BAM files have paired-end reads)\n\
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


static char complCode[128] = COMPLEMENT_LC;


/************************************************************************
 *   Read_T
 ************************************************************************/


#define T Read_T
typedef struct T *T;
struct T {
  char *queryseq;
  char *quality_string;
};


static void
Read_free (T *old) {
  FREE((*old)->queryseq);
  if ((*old)->quality_string != NULL) {
    FREE((*old)->quality_string);
  }
  FREE(*old);
  return;
}


static T
Read_create (bool *completep, T read, Bamline_T bamline) {
  unsigned int flag;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  char *shortread, *quality_string, *p, *q;
  char strand;
  int type;
  int cigar_querylength, mlength, querypos, i;

  flag = Bamline_flag(bamline);
  shortread = Bamline_read(bamline);
  quality_string = Bamline_quality_string(bamline);
  if (flag & QUERY_UNMAPPED) {
    read = (T) MALLOC(sizeof(struct T));
    read->queryseq = (char *) CALLOC(strlen(shortread)+1,sizeof(char));
    strcpy(read->queryseq,shortread);
    if (quality_string == NULL) {
      read->quality_string = (char *) NULL;
    } else {
      read->quality_string = (char *) CALLOC(strlen(quality_string)+1,sizeof(char));
      strcpy(read->quality_string,quality_string);
    }
    *completep = true;
    return read;
  }

  types = Bamline_cigar_types(bamline);
  npositions = Bamline_cigar_npositions(bamline);
  cigar_querylength = Bamline_cigar_querylength(bamline);

  if (read == NULL) {
    read = (T) MALLOC(sizeof(struct T));
    read->queryseq = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    for (i = 0; i < cigar_querylength; i++) {
      read->queryseq[i] = '*';
    }
    if (Bamline_quality_string(bamline) == NULL) {
      read->quality_string = (char *) NULL;
    } else {
      read->quality_string = (char *) CALLOC(cigar_querylength+1,sizeof(char));
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
  q = quality_string;
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
	  read->queryseq[querypos] = tolower(*p++);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos++;
	}
      } else {
      /* I and M */
	while (mlength-- > 0) {
	  read->queryseq[querypos] = toupper(*p++);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos++;
	}
      }

    } else {
      if (terminal_minlength >= 0 && type == 'S') {
	while (mlength-- > 0) {
	  read->queryseq[querypos] = tolower(complCode[(int) *p++]);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos--;
	}
      } else {
	/* I and M */
	while (mlength-- > 0) {
	  read->queryseq[querypos] = toupper(complCode[(int) *p++]);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos--;
	}
      }
    }
  }
  
  *completep = true;
  for (i = 0; i < cigar_querylength; i++) {
    if (read->queryseq[i] == '*') {
      *completep = false;
    }
  }

  return read;
}

static bool
Read_complete (T read) {
  unsigned int i;

  for (i = 0; i < strlen(read->queryseq); i++) {
    if (read->queryseq[i] == '*') {
      return false;
    }
  }

  return true;
}

#if 0
static int
Read_string_compare (const void *x, const void *y) {
  T a = (T) x;
  T b = (T) y;

  return strcmp(a->queryseq,b->queryseq);
}

/* This is the X31 hash function */
static unsigned int
Read_string_hash (const void *x) {
  T a = (T) x;
  unsigned int h = 0U;
  const char *p;
  
  for (p = a->queryseq; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}
#endif



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


/* Modifies queryseq */
static int
terminal_length (char *queryseq) {
  int terminal_length = 0, i;

  for (i = 0; i < (int) strlen(queryseq); i++) {
    if (islower(queryseq[i])) {
      terminal_length += 1;
      queryseq[i] = toupper(queryseq[i]);
    }
  }
  return terminal_length;
}


static void
print_fasta_single (FILE *fp, char *acc, T read) {
  if (terminal_minlength < 0 || terminal_length(read->queryseq) >= terminal_minlength) {
    fprintf(fp,">%s\n",acc);
    fprintf(fp,"%s\n",read->queryseq);
  }
  return;
}


static void
print_fasta_paired (FILE *fp, char *acc, T read1, T read2) {
  int terminal_length_1, terminal_length_2;

  if (single_end_output_p == true) {
    /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
    if (terminal_minlength < 0 || terminal_length(read1->queryseq) >= terminal_minlength) {
      fprintf(fp,">%s/1\n",acc);
      fprintf(fp,"%s\n",read1->queryseq);
    }
    if (terminal_minlength < 0 || terminal_length(read2->queryseq) >= terminal_minlength) {
      fprintf(fp,">%s/2\n",acc);
      fprintf(fp,"%s\n",read2->queryseq);
    }
  } else {
    if (terminal_minlength < 0) {
      fprintf(fp,">%s\n",acc);
      fprintf(fp,"%s\n",read1->queryseq);
      fprintf(fp,"%s\n",read2->queryseq);
    } else {
      terminal_length_1 = terminal_length(read1->queryseq);
      terminal_length_2 = terminal_length(read2->queryseq);
      if (terminal_length_1 >= terminal_minlength || terminal_length_2 >= terminal_minlength) {
	fprintf(fp,">%s\n",acc);
	fprintf(fp,"%s\n",read1->queryseq);
	fprintf(fp,"%s\n",read2->queryseq);
      }
    }
  }

  return;
}


static void
print_quality_string (FILE *fp, char *quality_string) {
  int i;

  for (i = 0; i < (int) strlen(quality_string); i++) {
    fprintf(fp,"%c",quality_string[i] + quality_shift);
  }
  fprintf(fp,"\n");

  return;
}



static void
print_fastq_single (FILE *fp, char *acc, T read) {

  if (terminal_minlength < 0 || terminal_length(read->queryseq) >= terminal_minlength) {
    fprintf(fp,"@%s\n",acc);
    fprintf(fp,"%s\n",read->queryseq);
    if (read->quality_string != NULL) {
      fprintf(fp,"+\n");
      print_quality_string(fp,read->quality_string);
    }
  }

  return;
}


static void
print_fastq_paired (FILE *fp1, FILE *fp2, char *acc, T read1, T read2) {
  int terminal_length_1, terminal_length_2;

  if (single_end_output_p == true) {
    /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
    if (terminal_minlength < 0 || terminal_length(read1->queryseq) >= terminal_minlength) {
      fprintf(fp1,"@%s/1\n",acc);
      fprintf(fp1,"%s\n",read1->queryseq);
      if (read1->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read1->quality_string);
      }
    }
    if (terminal_minlength < 0 || terminal_length(read2->queryseq) >= terminal_minlength) {
      fprintf(fp1,"@%s/2\n",acc);
      fprintf(fp1,"%s\n",read2->queryseq);
      if (read2->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read2->quality_string);
      }
    }

  } else {
    if (terminal_minlength < 0) {
      fprintf(fp1,"@%s/1\n",acc);
      fprintf(fp1,"%s\n",read1->queryseq);
      if (read1->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read1->quality_string);
      }
      fprintf(fp2,"@%s/2\n",acc);
      fprintf(fp2,"%s\n",read2->queryseq);
      if (read2->quality_string != NULL) {
	fprintf(fp2,"+\n");
	print_quality_string(fp2,read2->quality_string);
      }

    } else {
      terminal_length_1 = terminal_length(read1->queryseq);
      terminal_length_2 = terminal_length(read2->queryseq);
      if (terminal_length_1 >= terminal_minlength || terminal_length_2 >= terminal_minlength) {
	fprintf(fp1,"@%s/1\n",acc);
	fprintf(fp1,"%s\n",read1->queryseq);
	if (read1->quality_string != NULL) {
	  fprintf(fp1,"+\n");
	  print_quality_string(fp1,read1->quality_string);
	}

	fprintf(fp2,"@%s/2\n",acc);
	fprintf(fp2,"%s\n",read2->queryseq);
	if (read2->quality_string != NULL) {
	  fprintf(fp2,"+\n");
	  print_quality_string(fp2,read2->quality_string);
	}
      }

    }
  }

  return;
}


static void
parse_bam_input (FILE *fp1, FILE *fp2, Bamreader_T bamreader, Table_T read_table1, Table_T read_table2) {
  Bamline_T bamline;
  char *acc, *acc_copy;
  Read_T read1, read2;
  bool found1p, found2p, complete1p, complete2p;

  while ((bamline = Bamread_next_bamline(bamreader,/*minimum_mapq*/0,/*good_unique_mapq*/0,/*maximum_nhits*/1000000,
					 /*need_unique_p*/false,/*need_primary_p*/true,/*need_concordant_p*/false)) != NULL) {
    acc = Bamline_acc(bamline);
    if (Bamline_paired_read_p(bamline) == false) {
      read1 = (T) Table_get(read_table1,acc);
      found1p = (read1 ? true : false);
      read1 = Read_create(&complete1p,read1,bamline);
      if (complete1p == false) {
	Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
      } else {
	if (fastap == true) {
	  print_fasta_single(fp1,acc,read1);
	} else {
	  print_fastq_single(fp1,acc,read1);
	}
	if (found1p == true) {
	  acc_copy = Table_remove(read_table1,acc);
	  FREE(acc_copy);
	}
	Read_free(&read1);
      }

    } else if (Bamline_firstend_p(bamline) == true) {
      read1 = (T) Table_get(read_table1,acc);
      found1p = (read1 ? true : false);
      read1 = Read_create(&complete1p,read1,bamline);
      if (complete1p == false) {
	Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
      } else if ((read2 = (T) Table_get(read_table2,acc)) == NULL) {
	Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
      } else if (Read_complete(read2) == false) {
	Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
      } else {
	if (fastap == true) {
	  print_fasta_paired(fp1,acc,read1,read2);
	} else {
	  print_fastq_paired(fp1,fp2,acc,read1,read2);
	}
	if (found1p == true) {
	  acc_copy = Table_remove(read_table1,acc);
	  FREE(acc_copy);
	}
	acc_copy = Table_remove(read_table2,acc);
	FREE(acc_copy);
	Read_free(&read2);
	Read_free(&read1);
      }

    } else {
      read2 = (T) Table_get(read_table2,acc);
      found2p = (read2 ? true : false);
      read2 = Read_create(&complete2p,read2,bamline);
      if (complete2p == false) {
	Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
      } else if ((read1 = Table_get(read_table1,acc)) == NULL) {
	Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
      } else if (Read_complete(read1) == false) {
	Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
      } else {
	if (fastap == true) {
	  print_fasta_paired(fp1,acc,read1,read2);
	} else {
	  print_fastq_paired(fp1,fp2,acc,read1,read2);
	}
	if (found2p == true) {
	  acc_copy = Table_remove(read_table2,acc);
	  FREE(acc_copy);
	}
	acc_copy = Table_remove(read_table1,acc);
	FREE(acc_copy);
	Read_free(&read1);
	Read_free(&read2);
      }
    }

    Bamline_free(&bamline);
  }

  return;
}



int
main (int argc, char *argv[]) {
  FILE *fp1, *fp2 = NULL;
  char *filename;
  Bamreader_T bamreader;
  Table_T read_table1, read_table2;
  int i;

  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"o:A1S:N",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"fasta")) {
	fastap = true;

      } else if (!strcmp(long_name,"quality-protocol")) {
	if (!strcmp(optarg,"illumina")) {
	  quality_shift = 64;
	} else if (!strcmp(optarg,"sanger")) {
	  quality_shift = 33;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"version")) {
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

    case 'o': output_root = optarg; break;
    case 'A': allp = true; break;
    case '1': single_end_output_p = true; break;
    case 'S': terminal_minlength = atoi(optarg); break;
    case 'N': allp = false; break;
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (fastap == true) {
    fp1 = stdout;

  } else {
    filename = (char *) CALLOC(strlen(output_root) + strlen("_1.fq") + 1,sizeof(char));
    sprintf(filename,"%s_1.fq",output_root);
    fp1 = fopen(filename,"w");
    FREE(filename);

    if (single_end_output_p == false) {
      filename = (char *) CALLOC(strlen(output_root) + strlen("_2.fq") + 1,sizeof(char));
      sprintf(filename,"%s_2.fq",output_root);
      fp2 = fopen(filename,"w");
      FREE(filename);
    }
  }

  for (i = 0; i < argc; i++) {
    fprintf(stderr,"Processing file %s\n",argv[i]);
    read_table1 = Table_new(1000000,Table_string_compare,Table_string_hash);
    read_table2 = Table_new(1000000,Table_string_compare,Table_string_hash);

    bamreader = Bamread_new(argv[i]);
    parse_bam_input(fp1,fp2,bamreader,read_table1,read_table2);
    Bamread_free(&bamreader);

    Table_free(&read_table2);
    Table_free(&read_table1);
  }

  if (fp2 != NULL) {
    fclose(fp2);
  }
  if (fastap == false) {
    fclose(fp1);
  }

  return 0;
}

