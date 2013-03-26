static char rcsid[] = "$Id: bam_fastq.c 87713 2013-03-01 18:32:34Z twu $";
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
static char *error_root = NULL;
static bool fastap = false;
static bool allp = true;
static bool single_end_output_p = false; /* print paired-end BAM as single-end file */
static int terminal_minlength = -1;
static bool want_splicedp = false;

/* For Illumina, add 64.  For Sanger, add 33 */
static int quality_shift = 33;


static struct option long_options[] = {
  /* Input options */
#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif

  {"output", required_argument, 0, 'o'},	/* output_root */
  {"error", required_argument, 0, 'e'},	/* error_root */
  {"fasta", no_argument, 0, 0},			/* fastap */
  {"all", no_argument, 0, 'A'}, /* allp */
  {"single-end-output", no_argument, 0, '1'}, /* single_end_output_p */
  {"spliced", no_argument, 0, 0}, /* splicedp */
  {"terminal-minlength", required_argument, 0, 'S'}, /* terminal_minlength */
  {"nomappers", no_argument, 0, 'N'}, /* allp */
  {"quality-protocol", required_argument, 0, 0}, /* quality_shift */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_FASTQ -- Converts BAM results into FASTQ format\n");
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
Usage: bam_fastq [options] [-o <output>] <bam file...>\n\
\n\
Options:\n\
  --fasta                         Writes to stdout in FASTA format\n\
  -o, --output=STRING             Writes FASTQ to <output>_1.fq and <output>_2.fq (default=out)\n\
  -1, --single-end-output         Prints reads as single-end (even if BAM files have paired-end reads)\n\
  -S, --terminal-minlength=INT    Print only reads with terminal ends (cigar S length) equal to\n\
                                     or greater than this value (if not specified, ignores cigar S length)\n\
  --spliced                       Also print reads that contain a splice or are not mapped\n\
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


static bool strip_warning_p = false;
static bool stripp = true;

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


static char complCode[128] = COMPLEMENT_LC;


/************************************************************************
 *   Read_T
 ************************************************************************/


#define T Read_T
typedef struct T *T;
struct T {
  char *queryseq;
  char *quality_string;
  bool splicedp;
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
Read_create_bam (bool *completep, T read, Bamline_T bamline) {
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
    read->splicedp = true;
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
    read->splicedp = false;
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
      read->splicedp = true;
      
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


static T
Read_create_sam (bool *completep, T read, char *line) {
  unsigned int flag;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  char *shortread, *quality_string, *p, *q;
  char strand;
  int type;
  int cigar_querylength, mlength, querypos, i;

  Genomicpos_T chrpos_low;
  char *chr, *acc, *cigar, *auxinfo;
  int mapq;
  int readlength;

  auxinfo = Samread_parse_line(&acc,&flag,&mapq,&chr,&chrpos_low,&cigar,&readlength,&shortread,&quality_string,line);
  FREE(acc);
  FREE(chr);

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
    read->splicedp = true;
    *completep = true;

    FREE(cigar);
    FREE(quality_string);
    FREE(shortread);

    return read;
  }

  types = Samread_parse_cigar(&npositions,&cigar_querylength,cigar);
  FREE(cigar);

  if (read == NULL) {
    read = (T) MALLOC(sizeof(struct T));
    read->queryseq = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    for (i = 0; i < cigar_querylength; i++) {
      read->queryseq[i] = '*';
    }
    if (quality_string == NULL) {
      read->quality_string = (char *) NULL;
    } else {
      read->quality_string = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    }
    read->splicedp = false;
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
      read->splicedp = true;
      
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
  
  Uintlist_free(&npositions);
  Intlist_free(&types);
  FREE(quality_string);
  FREE(shortread);

  
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

/* Needed if terminal_minlength is set (which introduces the
   lower-case characters), and we don't call terminal_length first. */
static void
print_upper_case (FILE *fp, char *queryseq) {
  int i;

  for (i = 0; i < (int) strlen(queryseq); i++) {
    fputc(toupper(queryseq[i]),fp);
  }
  fputc('\n',fp);

  return;
}


static void
print_fasta_single (FILE *fp, char *acc, T read) {
  if (want_splicedp == true && read->splicedp == true) {
    fprintf(fp,">%s\n",acc);
    print_upper_case(fp,read->queryseq);
  } else if (terminal_minlength < 0 || terminal_length(read->queryseq) >= terminal_minlength) {
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
    if (want_splicedp == true && read1->splicedp == true) {
      fprintf(fp,">%s/1\n",acc);
      print_upper_case(fp,read1->queryseq);
    } else if (terminal_minlength < 0 || terminal_length(read1->queryseq) >= terminal_minlength) {
      fprintf(fp,">%s/1\n",acc);
      fprintf(fp,"%s\n",read1->queryseq);
    }
    if (want_splicedp == true && read2->splicedp == true) {
      fprintf(fp,">%s/2\n",acc);
      print_upper_case(fp,read2->queryseq);
    } else if (terminal_minlength < 0 || terminal_length(read2->queryseq) >= terminal_minlength) {
      fprintf(fp,">%s/2\n",acc);
      fprintf(fp,"%s\n",read2->queryseq);
    }
  } else {
    if (want_splicedp == true && (read1->splicedp == true || read2->splicedp == true)) {
      fprintf(fp,">%s\n",acc);
      print_upper_case(fp,read1->queryseq);
      print_upper_case(fp,read2->queryseq);
    } else if (terminal_minlength < 0) {
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


#define MONITOR_INTERVAL 1000000
static int nprinted = 0;
static int nprinted_single = 0;
static int nprinted_paired = 0;

static void
print_fastq_single (FILE *fp, char *acc, T read) {
  bool printp = false;

  if (want_splicedp == true && read->splicedp == true) {
    printp = true;
  } else if (terminal_minlength < 0 || terminal_length(read->queryseq) >= terminal_minlength) {
    printp = true;
  }

  if (printp == true) {
    nprinted_single++;
    if (++nprinted % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Output %d reads (%d single-end, %d paired-end)\n",
	      nprinted,nprinted_single,nprinted_paired);
    }

    fprintf(fp,"@%s\n",acc);
    print_upper_case(fp,read->queryseq);
    if (read->quality_string != NULL) {
      fprintf(fp,"+\n");
      print_quality_string(fp,read->quality_string);
    }
  }

  return;
}

static void
print_fastq_single_nomonitor (FILE *fp, char *acc, T read, char *suffix) {
  bool printp = false;

  if (want_splicedp == true && read->splicedp == true) {
    printp = true;
  } else if (terminal_minlength < 0 || terminal_length(read->queryseq) >= terminal_minlength) {
    printp = true;
  }

  if (printp == true) {
    fprintf(fp,"@%s%s\n",acc,suffix);
    print_upper_case(fp,read->queryseq);
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
  bool printp;


  if (single_end_output_p == true) {
    /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
    printp = false;
    if (want_splicedp == true && read1->splicedp == true) {
      printp = true;
    } else if (terminal_minlength < 0 || terminal_length(read1->queryseq) >= terminal_minlength) {
      printp = true;
    }
    if (printp == true) {
      nprinted_single++;
      if (++nprinted % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Output %d reads (%d single-end, %d paired-end)\n",
		nprinted,nprinted_single,nprinted_paired);
      }

      fprintf(fp1,"@%s/1\n",acc);
      print_upper_case(fp1,read1->queryseq);
      if (read1->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read1->quality_string);
      }
    }

    printp = false;
    if (want_splicedp == true && read2->splicedp == true) {
      printp = true;
    } else if (terminal_minlength < 0 || terminal_length(read2->queryseq) >= terminal_minlength) {
      printp = true;
    }
    if (printp == true) {
      nprinted_single++;
      if (++nprinted % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Output %d reads (%d single-end, %d paired-end)\n",
		nprinted,nprinted_single,nprinted_paired);
      }

      fprintf(fp1,"@%s/2\n",acc);
      print_upper_case(fp1,read2->queryseq);
      if (read2->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read2->quality_string);
      }
    }

  } else {
    printp = false;
    if (want_splicedp == true && (read1->splicedp == true || read2->splicedp == true)) {
      printp = true;
    } else if (terminal_minlength < 0) {
      printp = true;
    } else {
      terminal_length_1 = terminal_length(read1->queryseq);
      terminal_length_2 = terminal_length(read2->queryseq);
      if (terminal_length_1 >= terminal_minlength || terminal_length_2 >= terminal_minlength) {
	printp = true;
      }
    }

    if (printp == true) {
      nprinted_paired++;
      if (++nprinted % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Output %d reads (%d single-end, %d paired-end)\n",
		nprinted,nprinted_single,nprinted_paired);
      }

      fprintf(fp1,"@%s/1\n",acc);
      print_upper_case(fp1,read1->queryseq);
      if (read1->quality_string != NULL) {
	fprintf(fp1,"+\n");
	print_quality_string(fp1,read1->quality_string);
      }
      fprintf(fp2,"@%s/2\n",acc);
      print_upper_case(fp2,read2->queryseq);
      if (read2->quality_string != NULL) {
	fprintf(fp2,"+\n");
	print_quality_string(fp2,read2->quality_string);
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
  bool copyp;

  while ((bamline = Bamread_next_bamline(bamreader,/*minimum_mapq*/0,/*good_unique_mapq*/0,/*maximum_nhits*/1000000,
					 /*need_unique_p*/false,/*need_primary_p*/true,/*ignore_duplicates_p*/false,
					 /*need_concordant_p*/false)) != NULL) {
    acc = Bamline_acc(bamline);
    acc = strip_illumina_acc_ending(&copyp,acc);
    debug(printf("Acc is %s",acc));
    if (Bamline_paired_read_p(bamline) == false) {
      debug(printf(" single-end\n"));
      read1 = (T) Table_get(read_table1,acc);
      found1p = (read1 ? true : false);
      read1 = Read_create_bam(&complete1p,read1,bamline);
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
      debug(printf(" paired-end 1\n"));
      read1 = (T) Table_get(read_table1,acc);
      found1p = (read1 ? true : false);
      read1 = Read_create_bam(&complete1p,read1,bamline);
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
      debug(printf(" paired-end 2\n"));
      read2 = (T) Table_get(read_table2,acc);
      found2p = (read2 ? true : false);
      read2 = Read_create_bam(&complete2p,read2,bamline);
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

    if (copyp) {
      FREE(acc);
    }

    Bamline_free(&bamline);
  }

  return;
}

static void
parse_sam_input (Table_T read_table1, Table_T read_table2) {
  Read_T read1, read2;
  char line[1024000];
  char *acc, *acc_orig, *acc_copy;
  bool found1p, found2p, complete1p, complete2p;
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

      if (!(flag & PAIRED_READ)) {
	debug(printf(" single-end\n"));
	read1 = (T) Table_get(read_table1,acc);
	found1p = (read1 ? true : false);
	read1 = Read_create_sam(&complete1p,read1,line);
	if (complete1p == false) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else {
	  if (fastap == true) {
	    print_fasta_single(stdout,acc,read1);
	  } else {
	    print_fastq_single(stdout,acc,read1);
	  }
	  if (found1p == true) {
	    acc_copy = Table_remove(read_table1,acc);
	    FREE(acc_copy);
	  }
	  Read_free(&read1);
	}

      } else if (flag & FIRST_READ_P) {
	debug(printf(" paired-end 1\n"));
	read1 = (T) Table_get(read_table1,acc);
	found1p = (read1 ? true : false);
	read1 = Read_create_sam(&complete1p,read1,line);
	if (complete1p == false) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else if ((read2 = (T) Table_get(read_table2,acc)) == NULL) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else if (Read_complete(read2) == false) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else {
	  if (fastap == true) {
	    print_fasta_paired(stdout,acc,read1,read2);
	  } else {
	    print_fastq_paired(stdout,stdout,acc,read1,read2);
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

      } else if (flag & SECOND_READ_P) {
	debug(printf(" paired-end 2\n"));
	read2 = (T) Table_get(read_table2,acc);
	found2p = (read2 ? true : false);
	read2 = Read_create_sam(&complete2p,read2,line);
	if (complete2p == false) {
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
	} else if ((read1 = Table_get(read_table1,acc)) == NULL) {
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
	} else if (Read_complete(read1) == false) {
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
	} else {
	  if (fastap == true) {
	    print_fasta_paired(stdout,acc,read1,read2);
	  } else {
	    print_fastq_paired(stdout,stdout,acc,read1,read2);
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

      } else {
	fprintf(stderr,"Flag %u is paired (%u), but contains neither first_read nor second_read flag\n",
		flag,flag & PAIRED_READ);
	abort();
      }

      FREE(acc_orig);
    }
  }

  return;
}


static void
empty_table (FILE *fp, Table_T read_table, char *suffix) {
  char *acc, *acc_copy;
  char **accessions;
  Read_T read;
  int n, i;

  n = Table_length(read_table);
  accessions = (char **) Table_keys(read_table,(void *) NULL);
  for (i = 0; i < n; i++) {
    acc = accessions[i];
    read = (T) Table_get(read_table,acc);
    if (fastap == true) {
      print_fasta_single(fp,acc,read);
    } else if (nprinted_paired > 0) {
      print_fastq_single_nomonitor(fp,acc,read,suffix);
    } else {
      print_fastq_single_nomonitor(fp,acc,read,"");
    }
    acc_copy = Table_remove(read_table,acc);
    FREE(acc_copy);
    Read_free(&read);
  }
  FREE(accessions);

  return;
}



int
main (int argc, char *argv[]) {
  FILE *fp1, *fp2 = NULL, *err1 = NULL, *err2 = NULL;
  char *filename;
  Bamreader_T bamreader;
  Table_T read_table1, read_table2;
  int i, n;

  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"o:e:A1S:N",
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

      } else if (!strcmp(long_name,"fasta")) {
	fastap = true;

      } else if (!strcmp(long_name,"terminal-minlength")) {
	terminal_minlength = atoi(optarg);

      } else if (!strcmp(long_name,"spliced")) {
	want_splicedp = true;

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
    case 'e': error_root = optarg; break;
    case 'A': allp = true; break;
    case '1': single_end_output_p = true; break;
    case 'S': terminal_minlength = atoi(optarg); break;
    case 'N': allp = false; break;

    case '?': fprintf(stderr,"For usage, run 'bam_fastq --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (error_root != NULL) {
    filename = (char *) CALLOC(strlen(error_root) + strlen("_1.fq") + 1,sizeof(char));
    sprintf(filename,"%s_1.fq",error_root);
    err1 = fopen(filename,"w");
    FREE(filename);

    filename = (char *) CALLOC(strlen(error_root) + strlen("_2.fq") + 1,sizeof(char));
    sprintf(filename,"%s_2.fq",error_root);
    err2 = fopen(filename,"w");
    FREE(filename);
  }

  if (argc == 0) {
    fastap = true;
    fp1 = stdout;
    read_table1 = Table_new(1000000,Table_string_compare,Table_string_hash);
    read_table2 = Table_new(1000000,Table_string_compare,Table_string_hash);
    parse_sam_input(read_table1,read_table2);

  } else {
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
    }
  }

  if ((n = Table_length(read_table1)) > 0) {
    fprintf(stderr,"Warning: for file %s, %d first-end reads were not completed or paired",argv[i],n);
    if (fastap == true) {
      empty_table(stdout,read_table1,"/1");
    } else if (err1 == NULL) {
      fprintf(stderr," (to write to file, use the -e flag)\n");
    } else {
      fprintf(stderr," (written to file %s_1.fq)\n",error_root);
      empty_table(err1,read_table1,"/1");
    }
  }
  if ((n = Table_length(read_table2)) > 0) {
    fprintf(stderr,"Warning: for file %s, %d second-end reads were not completed or paired",argv[i],n);
    if (fastap == true) {
      empty_table(stdout,read_table2,"/2");
    } else if (err2 == NULL) {
      fprintf(stderr," (to write to file, use the -e flag)\n");
    } else {
      fprintf(stderr," (written to file %s_2.fq)\n",error_root);
      empty_table(err2,read_table2,"/2");
    }
  }

  Table_free(&read_table2);
  Table_free(&read_table1);


  if (err2 != NULL) {
    fclose(err2);
  }
  if (err1 != NULL) {
    fclose(err1);
  }

  if (fp2 != NULL) {
    fclose(fp2);
  }
  if (fastap == false) {
    fclose(fp1);
  }

  if (fastap == false) {
    fprintf(stderr,"Output %d reads total (%d single-end, %d paired-end)\n",
	    nprinted,nprinted_single,nprinted_paired);
  }

  return 0;
}

