static char rcsid[] = "$Id: extract-concordant.c 143394 2014-08-05 16:01:50Z twu $";
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

#include "iit-read.h"
#include "samflags.h"
#include "bamread.h"
#include "table.h"
#include "complement.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

static List_T chr_omit = NULL;
static char *chr_add_prefix = NULL;
static char *unpaired_filename = NULL;
static char *splices_filename = NULL;
static bool fastap = false;
static char *read_output_filename = NULL;
static bool add_bounds_p = false;

/* For Illumina, add 64.  For Sanger, add 33 */
static int quality_shift = 33;

static IIT_T transcript_regions_iit = NULL;


static struct option long_options[] = {
  /* Input options */
  {"transcript", required_argument, 0, 't'}, /* transcript_regions_iitfile */
  {"omit", required_argument, 0, 'x'}, /* chr_omit */
  {"add-prefix", required_argument, 0, 'P'}, /* chr_add_prefix */
  {"unpaired", required_argument, 0, 'u'},   /* unpaired_filename */
  {"splices", required_argument, 0, 's'},    /* splices_filename */
  {"fasta", no_argument, 0, 0},	     /* fastap */
  {"read-output", required_argument, 0, 'f'},	     /* read_output_filename */
  {"add-bounds", no_argument, 0, 0},	     /* add_bounds_p */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"EXTRACT-CONCORDANT -- Finds intertranscript splices in concordant GSNAP output\n");
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
Usage: extract-concordant [options] [-o <output>] <bam file...>\n\
\n\
Options:\n\
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

static void
print_substring (FILE *fp, char *queryseq, int querystart, int queryend) {
  int i;

  for (i = querystart; i < queryend; i++) {
    fprintf(fp,"%c",queryseq[i]);
  }
  return;
}

static char complCode[128] = COMPLEMENT_LC;

static void
print_substring_revcomp (FILE *fp, char *queryseq, int querystart, int queryend) {
  int i;

  for (i = queryend - 1; i >= querystart; i--) {
    fprintf(fp,"%c",complCode[(int) queryseq[i]]);
  }
  return;
}


static void
print_quality (FILE *fp, char *quality_string, int querystart, int queryend) {
  int i;

  for (i = querystart; i < queryend; i++) {
    fprintf(fp,"%c",quality_string[i] + quality_shift);
  }
  return;
}

static void
print_quality_reverse (FILE *fp, char *quality_string, int querystart, int queryend) {
  int i;

  for (i = queryend - 1; i >= querystart; i--) {
    fprintf(fp,"%c",quality_string[i] + quality_shift);
  }
  return;
}




/************************************************************************
 *   Splice_T (simplified from bam_fastq)
 ************************************************************************/

typedef struct Splice_T *Splice_T;
struct Splice_T {
  int splice_querypos;
  Genomicpos_T splicestart;
  Genomicpos_T spliceend;
};

static void
Splice_free (Splice_T *old) {
  FREE(*old);
  return;
}

static Splice_T
Splice_new (int splice_querypos, Genomicpos_T splicestart, Genomicpos_T spliceend) {
  Splice_T new = (Splice_T) MALLOC(sizeof(*new));

  new->splice_querypos = splice_querypos;
  new->splicestart = splicestart;
  new->spliceend = spliceend;
  return new;
}


/************************************************************************
 *   Read_T (simplified from bam_fastq)
 ************************************************************************/


#define T Read_T
typedef struct T *T;
struct T {
  int flags;
  char *queryseq;
  char *quality_string;
  char *cigar_string;

  /* A read could provide evidence for multiple types of splices */
  List_T splices;

  char genestrand;
  int querylength;

  char *chr;
  char readstrand;
  Genomicpos_T chrpos_low;
  Genomicpos_T chrpos_high;
  Genomicpos_T querystart_chrpos;
  Genomicpos_T queryend_chrpos;
};


static void
Read_free (T *old) {
  List_T p;
  Splice_T splice;

  for (p = (*old)->splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    Splice_free(&splice);
  }
  List_free(&(*old)->splices);

  FREE((*old)->cigar_string);
  FREE((*old)->chr);
  FREE((*old)->queryseq);
  if ((*old)->quality_string != NULL) {
    FREE((*old)->quality_string);
  }
  FREE(*old);
  return;
}


static void
Read_add_chrinfo (Read_T read, int flags, char *cigar_string, char *chr, Genomicpos_T chrpos_low, Genomicpos_T chrpos_high) {

  read->flags = flags;
  read->cigar_string = cigar_string;
  if (chr == NULL) {
    read->readstrand = '.';
    read->chr = (char *) CALLOC(2,sizeof(char));
    strcpy(read->chr,"*");
    read->chrpos_low = read->chrpos_high = 0;
    read->querystart_chrpos = read->queryend_chrpos = 0;
  } else {
    read->chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(read->chr,chr);
    read->chrpos_low = chrpos_low;
    read->chrpos_high = chrpos_high;
    if (flags & QUERY_MINUSP) {
      read->readstrand = '-';
      read->querystart_chrpos = chrpos_high;
      read->queryend_chrpos = chrpos_low;
    } else {
      read->readstrand = '+';
      read->querystart_chrpos = chrpos_low;
      read->queryend_chrpos = chrpos_high;
    }
  }

  return;
}


static T
Read_create_bam (bool *completep, T read, Bamline_T bamline, int divnoT) {
  unsigned int flag;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  Genomicpos_T chrpos;
  char *shortread, *quality_string, *cigar_string, *p, *q;
  char strand;
  int type;
  int cigar_querylength, mlength, querypos, i;
  /* bool insidep; */

  flag = Bamline_flag(bamline);

  shortread = Bamline_read(bamline);
  quality_string = Bamline_quality_string(bamline);
  if (flag & QUERY_UNMAPPED) {
    abort();
    read = (T) MALLOC(sizeof(struct T));
    read->genestrand = Bamline_splice_strand(bamline);
    read->queryseq = (char *) CALLOC(strlen(shortread)+1,sizeof(char));
    strcpy(read->queryseq,shortread);
    if (quality_string == NULL) {
      read->quality_string = (char *) NULL;
    } else {
      read->quality_string = (char *) CALLOC(strlen(quality_string)+1,sizeof(char));
      strcpy(read->quality_string,quality_string);
    }
    read->splices = (List_T) NULL;
    cigar_string = Bamline_cigar_string(bamline);
    Read_add_chrinfo(read,flag,cigar_string,Bamline_chr(bamline),
		     Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    *completep = true;
    return read;
  }

  types = Bamline_cigar_types(bamline);
  npositions = Bamline_cigar_npositions(bamline);
  cigar_querylength = Bamline_cigar_querylength(bamline);

  if (read == NULL) {
    /* Initialize */
    read = (T) MALLOC(sizeof(struct T));
    read->genestrand = Bamline_splice_strand(bamline);
    read->queryseq = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    for (i = 0; i < cigar_querylength; i++) {
      read->queryseq[i] = '*';
    }
    read->querylength = cigar_querylength;

    if (Bamline_quality_string(bamline) == NULL) {
      read->quality_string = (char *) NULL;
    } else {
      read->quality_string = (char *) CALLOC(cigar_querylength+1,sizeof(char));
    }
    read->splices = (List_T) NULL;
    cigar_string = Bamline_cigar_string(bamline);
    Read_add_chrinfo(read,flag,cigar_string,Bamline_chr(bamline),
		     Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
  }

  if (flag & QUERY_MINUSP) {
    strand = '-';
    querypos = cigar_querylength;
  } else {
    strand = '+';
    querypos = 0;
  }



  p = shortread;
  q = quality_string;
  chrpos = Bamline_chrpos_low(bamline);

  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);

    if (type == 'H') {
      if (strand == '+') {
	querypos += mlength;
      } else {
	querypos -= mlength;
      }
      /* Does not affect chrpos */

    } else if (type == 'N') {
      if (read->genestrand == '+') {
	/* Subtract 1U from chrpos to go from intronstart to exonend */
	if (IIT_contains_region_with_divno_signed(transcript_regions_iit,divnoT,
						  chrpos-1U,chrpos+mlength,/*sign*/+1) == true) {
	  /* Splice is within a transcript interval.  Ignore. */

	} else {
	  /* altexon, readthrough, or distant */
	  read->splices = List_push(read->splices,
				    (void *) Splice_new(querypos,/*splicestart*/chrpos - 1U,/*spliceend*/chrpos + mlength));
	}

      } else if (read->genestrand == '-') {
	if (IIT_contains_region_with_divno_signed(transcript_regions_iit,divnoT,
						  chrpos-1U,chrpos+mlength,/*sign*/-1) == true) {
	  /* Splice is within a transcript interval.  Ignore. */

	} else {
	  /* altexon, readthrough, or distant */
	  read->splices = List_push(read->splices,
				    (void *) Splice_new(querypos,/*splicestart*/chrpos + mlength,/*spliceend*/chrpos - 1U));
	}
      }
      chrpos += mlength;

    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'D') {
      chrpos += mlength;

    } else if (strand == '+') {
      if (type == 'M') {
	chrpos += mlength;
	while (mlength-- > 0) {
	  read->queryseq[querypos] = toupper(*p++);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos++;
	}
      } else {
	/* I and S.  Do not affect chrpos */
	while (mlength-- > 0) {
	  read->queryseq[querypos] = tolower(*p++);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	  querypos++;
	}
      }

    } else {
      if (type == 'M') {
	chrpos += mlength;
	while (mlength-- > 0) {
	  querypos--;
	  read->queryseq[querypos] = toupper(complCode[(int) *p++]);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	}
      } else {
	/* I and S.  Do not affect chrpos */
	while (mlength-- > 0) {
	  querypos--;
	  read->queryseq[querypos] = tolower(complCode[(int) *p++]);
	  if (q) {
	    read->quality_string[querypos] = *q++;
	  }
	}
	/* Does not affect chrpos */
      }
    }
  }

  
  *completep = true;
  for (i = 0; i < cigar_querylength; i++) {
    if (read->queryseq[i] == '*') {
      *completep = false;
    }
  }

#if 0
  if (insidep == true && Bamline_paired_read_p(bamline) == true) {
    if (Bamline_lowend_p(bamline) == true) {
      read->lowp = true;
    } else {
      read->highp = true;
    }
  }
#endif

  return read;
}


static void
print_splices (FILE *fp, char *acc, Read_T read, List_T splices, char *acc_ending) {
  List_T p;
  Splice_T splice;
  char readstrand;

  if (read->flags & SECOND_READ_P) {
    /* Invert second read */
    if (read->readstrand == '+') {
      readstrand = '-';
    } else {
      readstrand = '+';
    }
  } else {
    readstrand = read->readstrand;
  }
    
  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);

    if (read->genestrand == '+') {
      fprintf(fp,"%s\t%c\t%u\t%u\t+\t",read->chr,readstrand,read->chrpos_low,splice->splicestart);
      fprintf(fp,"%s\t",read->cigar_string);

      fprintf(fp,"%s\t%c\t%u\t%u\t+\t",read->chr,readstrand,read->chrpos_high,splice->spliceend);
      fprintf(fp,"%s\t",read->cigar_string);

    } else if (read->genestrand == '-') {
      fprintf(fp,"%s\t%c\t%u\t%u\t-\t",read->chr,readstrand,read->chrpos_high,splice->splicestart);
      fprintf(fp,"%s\t",read->cigar_string);

      fprintf(fp,"%s\t%c\t%u\t%u\t-\t",read->chr,readstrand,read->chrpos_low,splice->spliceend);
      fprintf(fp,"%s\t",read->cigar_string);
    }

    fprintf(fp,"%s%s",acc,acc_ending);
    fprintf(fp,"\t");

    /* Print from distal end towards splice */
    if (read->genestrand == '+') {
      if (read->readstrand == '+') {
	print_substring(fp,read->queryseq,/*querystart*/0,/*queryend*/splice->splice_querypos);
	fprintf(fp,"\t");
	print_substring_revcomp(fp,read->queryseq,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);

	fprintf(fp,"\t");
	if (read->quality_string == NULL) {
	  fprintf(fp,"*");
	  fprintf(fp,"\t");
	  fprintf(fp,"*");
	} else {
	  print_quality(fp,read->quality_string,/*querystart*/0,/*queryend*/splice->splice_querypos);
	  fprintf(fp,"\t");
	  print_quality_reverse(fp,read->quality_string,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	}

      } else {
	print_substring_revcomp(fp,read->queryseq,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	fprintf(fp,"\t");
	print_substring(fp,read->queryseq,/*querystart*/0,/*queryend*/splice->splice_querypos);

	fprintf(fp,"\t");
	if (read->quality_string == NULL) {
	  fprintf(fp,"*");
	  fprintf(fp,"\t");
	  fprintf(fp,"*");
	} else {
	  print_quality_reverse(fp,read->quality_string,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	  fprintf(fp,"\t");
	  print_quality(fp,read->quality_string,/*querystart*/0,/*queryend*/splice->splice_querypos);
	}

      }

    } else {
      if (read->readstrand == '+') {
	print_substring_revcomp(fp,read->queryseq,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	fprintf(fp,"\t");
	print_substring(fp,read->queryseq,/*querystart*/0,/*queryend*/splice->splice_querypos);

	fprintf(fp,"\t");
	if (read->quality_string == NULL) {
	  fprintf(fp,"*");
	  fprintf(fp,"\t");
	  fprintf(fp,"*");
	} else {
	  print_quality_reverse(fp,read->quality_string,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	  fprintf(fp,"\t");
	  print_quality(fp,read->quality_string,/*querystart*/0,/*queryend*/splice->splice_querypos);
	}

      } else {
	print_substring(fp,read->queryseq,/*querystart*/0,/*queryend*/splice->splice_querypos);
	fprintf(fp,"\t");
	print_substring_revcomp(fp,read->queryseq,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);

	fprintf(fp,"\t");
	if (read->quality_string == NULL) {
	  fprintf(fp,"*");
	  fprintf(fp,"\t");
	  fprintf(fp,"*");
	} else {
	  print_quality(fp,read->quality_string,/*querystart*/0,/*queryend*/splice->splice_querypos);
	  fprintf(fp,"\t");
	  print_quality_reverse(fp,read->quality_string,/*querystart*/splice->splice_querypos,/*queryend*/read->querylength);
	}
      }
    }

    fprintf(fp,"\n");
  }

  return;
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
print_quality_string (FILE *fp, char *quality_string) {
  int i;

  for (i = 0; i < (int) strlen(quality_string); i++) {
    fprintf(fp,"%c",quality_string[i] + quality_shift);
  }
  fprintf(fp,"\n");

  return;
}


static void
print_fasta_single (FILE *fp, char *acc, T read) {
  fprintf(fp,">%s",acc);
  if (add_bounds_p == true) {
    fprintf(fp,"(%c%s:%u,%c%s:%u)",
	    read->readstrand,read->chr,read->querystart_chrpos,
	    read->readstrand,read->chr,read->queryend_chrpos);
  }
  fprintf(fp,"\n");
  print_upper_case(fp,read->queryseq);

  return;
}

static void
print_fastq_single (FILE *fp, char *acc, T read) {
  fprintf(fp,"@%s",acc);
  if (add_bounds_p == true) {
    fprintf(fp,"(%c%s:%u,%c%s:%u)",
	    read->readstrand,read->chr,read->querystart_chrpos,
	    read->readstrand,read->chr,read->queryend_chrpos);
  }
  fprintf(fp,"\n");
  print_upper_case(fp,read->queryseq);
  if (read->quality_string != NULL) {
    fprintf(fp,"+\n");
    print_quality_string(fp,read->quality_string);
  }

  return;
}


/* Taken from single_end_output branch from bam_fastq.c */
/* Bounds are from given read querystart to mate querystart */
static void
print_fasta_paired (FILE *fp, char *acc, T read, T mate, char *acc_ending) {
  /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
  fprintf(fp,">%s",acc);
  if (add_bounds_p == true) {
    fprintf(fp,"(%c%s:%u,%c%s:%u..%c%s:%u)",
	    read->readstrand,read->chr,read->querystart_chrpos,
	    read->readstrand,read->chr,read->queryend_chrpos,
	    mate->readstrand,mate->chr,mate->querystart_chrpos);
  }
  fprintf(fp,"%s\n",acc_ending);
  print_upper_case(fp,read->queryseq);

  return;
}


static void
print_fastq_paired (FILE *fp, char *acc, T read, T mate, char *acc_ending) {
  /* Important to add the /1 and /2 endings for the remove-intragenic program, which also requires unique alignments */
  fprintf(fp,"@%s",acc);
  if (add_bounds_p == true) {
    fprintf(fp,"(%c%s:%u,%c%s:%u..%c%s:%u)",
	    read->readstrand,read->chr,read->querystart_chrpos,
	    read->readstrand,read->chr,read->queryend_chrpos,
	    mate->readstrand,mate->chr,mate->querystart_chrpos);
  }
  fprintf(fp,"%s\n",acc_ending);
  print_upper_case(fp,read->queryseq);
  if (read->quality_string != NULL) {
    fprintf(fp,"+\n");
    print_quality_string(fp,read->quality_string);
  }

  return;
}


/* Always print low coord read first, then high coord read.  This
   helps add-exons-unpaired figure out which is the donor and which is
   the acceptor, with the rule: genestrand + and + => first is donor;
   genestrand - and - => second is donor */
static void
print_unpaired_old (FILE *fp, char *acc, Bamline_T bamline_low, Bamline_T bamline_high, char *chr) {
  unsigned int flag_low, flag_high;
  char genestrand_low, genestrand_high;
  char readstrand_low, readstrand_high;

  flag_low = Bamline_flag(bamline_low);
  flag_high = Bamline_flag(bamline_high);

  if (flag_low & SECOND_READ_P) {
    /* Invert second read */
    readstrand_low = (flag_low & QUERY_MINUSP) ? '+' : '-';
  } else {
    readstrand_low = (flag_low & QUERY_MINUSP) ? '-' : '+';
  }
  if (flag_high & SECOND_READ_P) {
    /* Invert second read */
    readstrand_high = (flag_high & QUERY_MINUSP) ? '+' : '-';
  } else {
    readstrand_high = (flag_high & QUERY_MINUSP) ? '-' : '+';
  }

  genestrand_low = Bamline_splice_strand(bamline_low);
  genestrand_high = Bamline_splice_strand(bamline_high);

  if (flag_low & FIRST_READ_P) {
    fprintf(fp,"%s\t%c\t%u\t%u\t",chr,readstrand_low,Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline_low));
    if (genestrand_low == '+') {
      fprintf(fp,"+\t");
    } else if (genestrand_low == '-') {
      fprintf(fp,"-\t");
    } else {
      fprintf(fp,".\t");
    }
    Bamread_print_cigar(fp,bamline_low);
    fprintf(fp,"\t");

    fprintf(fp,"%s\t%c\t%u\t%u\t",chr,readstrand_high,Bamline_chrpos_high(bamline_high),Bamline_chrpos_low(bamline_high));
    if (genestrand_high == '+') {
      fprintf(fp,"+\t");
    } else if (genestrand_high == '-') {
      fprintf(fp,"-\t");
    } else {
      fprintf(fp,".\t");
    }
    Bamread_print_cigar(fp,bamline_high);
    fprintf(fp,"\t");

  } else {
    fprintf(fp,"%s\t%c\t%u\t%u\t",chr,readstrand_high,Bamline_chrpos_high(bamline_high),Bamline_chrpos_low(bamline_high));
    if (genestrand_high == '+') {
      fprintf(fp,"+\t");
    } else if (genestrand_high == '-') {
      fprintf(fp,"-\t");
    } else {
      fprintf(fp,".\t");
    }
    Bamread_print_cigar(fp,bamline_high);
    fprintf(fp,"\t");

    fprintf(fp,"%s\t%c\t%u\t%u\t",chr,readstrand_low,Bamline_chrpos_low(bamline_low),Bamline_chrpos_high(bamline_low));
    if (genestrand_low == '+') {
      fprintf(fp,"+\t");
    } else if (genestrand_low == '-') {
      fprintf(fp,"-\t");
    } else {
      fprintf(fp,".\t");
    }
    Bamread_print_cigar(fp,bamline_low);
    fprintf(fp,"\t");
  }

  fprintf(fp,"%s",acc);
  fprintf(fp,"\n");

  return;
}


/* Handle entire chromosomal region spanned by paired-end read */

/* Always print low coord read first, then high coord read.  This
   helps add-exons-unpaired figure out which is the donor and which is
   the acceptor, with the rule: genestrand + and + => first is donor;
   genestrand - and - => second is donor */
static void
print_unpaired (FILE *fp, char *acc, Read_T read1, Read_T read2, int divnoT) {

  if (read1->chrpos_high < read2->chrpos_low) {
    if (IIT_contains_region_with_divno(transcript_regions_iit,divnoT,
				       read1->chrpos_high,read2->chrpos_low) == true) {
      /* Entire span is within a transcript.  Ignore. */
    } else {
      fprintf(fp,"%s\t%c\t%u\t%u\t",read1->chr,read1->readstrand,read1->querystart_chrpos,read1->queryend_chrpos);
      if (read1->genestrand == '+') {
	fprintf(fp,"+\t");
      } else if (read1->genestrand == '-') {
	fprintf(fp,"-\t");
      } else {
	fprintf(fp,".\t");
      }
      fprintf(fp,"%s\t",read1->cigar_string);

      fprintf(fp,"%s\t%c\t%u\t%u\t",read2->chr,read2->readstrand,read2->querystart_chrpos,read2->queryend_chrpos);
      if (read2->genestrand == '+') {
	fprintf(fp,"+\t");
      } else if (read2->genestrand == '-') {
	fprintf(fp,"-\t");
      } else {
	fprintf(fp,".\t");
      }
      fprintf(fp,"%s\t",read2->cigar_string);
      fprintf(fp,"%s",acc);
      fprintf(fp,"\n");
    }

  } else if (read2->chrpos_high < read1->chrpos_low) {
    if (IIT_contains_region_with_divno(transcript_regions_iit,divnoT,
				       read2->chrpos_high,read1->chrpos_low) == true) {
      /* Entire span is within a transcript.  Ignore. */
    } else {
      fprintf(fp,"%s\t%c\t%u\t%u\t",read2->chr,read2->readstrand,read2->querystart_chrpos,read2->queryend_chrpos);
      if (read2->genestrand == '+') {
	fprintf(fp,"+\t");
      } else if (read2->genestrand == '-') {
	fprintf(fp,"-\t");
      } else {
	fprintf(fp,".\t");
      }
      fprintf(fp,"%s\t",read2->cigar_string);
      
      fprintf(fp,"%s\t%c\t%u\t%u\t",read1->chr,read1->readstrand,read1->querystart_chrpos,read1->queryend_chrpos);
      if (read1->genestrand == '+') {
	fprintf(fp,"+\t");
      } else if (read1->genestrand == '-') {
	fprintf(fp,"-\t");
      } else {
	fprintf(fp,".\t");
      }
      fprintf(fp,"%s\t",read1->cigar_string);
      fprintf(fp,"%s",acc);
      fprintf(fp,"\n");
    }
  }

  return;
}



static char *
string_copy (char *string) {
  char *copy;

  copy = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(copy,string);
  return copy;
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


static void
parse_bam_input (FILE *unpaired_fp, FILE *splices_fp, FILE *read_output_fp,
		 Bamreader_T bamreader, Table_T read_table1, Table_T read_table2,
		 bool *chr_omit_p) {
  Bamline_T bamline;
  char *acc, *acc_copy;
  Read_T read1, read2;
  bool found1p, found2p, complete1p, complete2p;
  bool copyp;

  char *chr = NULL, *chrtemp, *last_chr = NULL;
  int divnoT;

  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					 /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/true,
					 /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    if (chr_add_prefix != NULL) {
      FREE(chr);
      chrtemp = Bamline_chr(bamline);
      chr = (char *) CALLOC(strlen(chr_add_prefix)+strlen(chrtemp)+1,sizeof(char));
      sprintf(chr,"%s%s",chr_add_prefix,chrtemp);
    } else {
      chr = Bamline_chr(bamline);
    }

    if (last_chr == NULL) {
      divnoT = IIT_divint(transcript_regions_iit,chr);
      last_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
      strcpy(last_chr,chr);
    } else if (strcmp(chr,last_chr) != 0) {
      FREE(last_chr);
      divnoT = IIT_divint(transcript_regions_iit,chr);
      last_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
      strcpy(last_chr,chr);
    }

    if (chr_omit_p == NULL || chr_omit_p[divnoT] == false) {
      acc = Bamline_acc(bamline);
      acc = strip_illumina_acc_ending(&copyp,acc);
      debug(printf("Acc is %s",acc));
      if (Bamline_paired_read_p(bamline) == false) {
	debug(printf(" single-end\n"));
	read1 = (T) Table_get(read_table1,acc);
	found1p = (read1 ? true : false);
	read1 = Read_create_bam(&complete1p,read1,bamline,divnoT);
	if (complete1p == false) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else {
	  if (read1->splices != NULL) {
	    print_splices(splices_fp,acc,read1,read1->splices,"");
	    if (fastap == true) {
	      print_fasta_single(read_output_fp,acc,read1);
	    } else {
	      print_fastq_single(read_output_fp,acc,read1);
	    }
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
	read1 = Read_create_bam(&complete1p,read1,bamline,divnoT);
	if (complete1p == false) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
	} else if ((read2 = (T) Table_get(read_table2,acc)) == NULL) {
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
#if 0
	} else if (Read_complete(read2) == false) {
	  /* Not possible in concordant_uniq */
	  Table_put(read_table1,found1p ? acc : string_copy(acc),read1);
#endif
	} else {
	  if (read1->splices != NULL) {
	    print_splices(splices_fp,acc,read1,read1->splices,"/1");
	    if (fastap == true) {
	      print_fasta_paired(read_output_fp,acc,read1,read2,"/1");
	    } else {
	      print_fastq_paired(read_output_fp,acc,read1,read2,"/1");
	    }
	  }
	  if (read2->splices != NULL) {
	    print_splices(splices_fp,acc,read2,read2->splices,"/2");
	    if (fastap == true) {
	      print_fasta_paired(read_output_fp,acc,read2,read1,"/2");
	    } else {
	      print_fastq_paired(read_output_fp,acc,read2,read1,"/2");
	    }
	  }
	  print_unpaired(unpaired_fp,acc,read1,read2,divnoT);

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
	read2 = Read_create_bam(&complete2p,read2,bamline,divnoT);
	if (complete2p == false) {
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
	} else if ((read1 = Table_get(read_table1,acc)) == NULL) {
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
#if 0
	} else if (Read_complete(read1) == false) {
	  /* Not possible in concordant_uniq */
	  Table_put(read_table2,found2p ? acc : string_copy(acc),read2);
#endif
	} else {
	  if (read1->splices != NULL) {
	    print_splices(splices_fp,acc,read1,read1->splices,"/1");
	    if (fastap == true) {
	      print_fasta_paired(read_output_fp,acc,read1,read2,"/1");
	    } else {
	      print_fastq_paired(read_output_fp,acc,read1,read2,"/1");
	    }
	  }
	  if (read2->splices != NULL) {
	    print_splices(splices_fp,acc,read2,read2->splices,"/2");
	    if (fastap == true) {
	      print_fasta_paired(read_output_fp,acc,read2,read1,"/2");
	    } else {
	      print_fastq_paired(read_output_fp,acc,read2,read1,"/2");
	    }
	  }
	  print_unpaired(unpaired_fp,acc,read1,read2,divnoT);

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
    }

    if (copyp) {
      FREE(acc);
    }

    Bamline_free(&bamline);
  }

  if (chr_add_prefix != NULL) {
    FREE(chr);
  }

  if (last_chr != NULL) {
    FREE(last_chr);
  }

  return;
}



int
main (int argc, char *argv[]) {
  FILE *unpaired_fp = NULL, *splices_fp = NULL, *read_output_fp = NULL;
  char *transcript_regions_iitfile = NULL, *chr;
  Bamreader_T bamreader;
  Table_T read_table1, read_table2;
  List_T p;
  bool *chr_omit_p;
  int divno;
  int n, i;

  int opt;
  extern int optind;
  /* extern char *optarg; */
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"u:s:f:t:x:P:",
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
      } else if (!strcmp(long_name,"add-bounds")) {
	add_bounds_p = true;
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'bam_fasta --help'",long_name);
	exit(9);
      }
      break;

    case 'u': unpaired_filename = optarg; break;
    case 's': splices_filename = optarg; break;
    case 'f': read_output_filename = optarg; break;
    case 't': transcript_regions_iitfile = optarg; break;
    case 'x': chr_omit = List_from_string(optarg); break;
    case 'P': chr_add_prefix = optarg; break;

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      

  if (transcript_regions_iitfile == NULL) {
    chr_omit_p = (bool *) NULL;
  } else {
    transcript_regions_iit = IIT_read(transcript_regions_iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true);
    chr_omit_p = (bool *) CALLOC(IIT_ndivs(transcript_regions_iit)+1,sizeof(bool));
    for (p = chr_omit; p != NULL; p = List_next(p)) {
      chr = (char *) List_head(p);
      if ((divno = IIT_divint(transcript_regions_iit,chr)) < 0) {
	fprintf(stderr,"Chromosome %s not found in IIT file\n",chr);
      } else {
	chr_omit_p[divno] = true;
      }
      FREE(chr);
    }
    List_free(&chr_omit);
  }

  if (read_output_filename == NULL) {
    fprintf(stderr,"Need to specify fasta file with -f flag\n");
    exit(9);
  } else {
    read_output_fp = fopen(read_output_filename,"a");
  }

  if (unpaired_filename == NULL) {
    fprintf(stderr,"Need to specify unpaired output file with -u flag\n");
    exit(9);
  } else {
    unpaired_fp = fopen(unpaired_filename,"a");
  }

  if (splices_filename == NULL) {
    fprintf(stderr,"Need to specify splices output file with -s flag\n");
    exit(9);
  } else {
    splices_fp = fopen(splices_filename,"a");
  }

  for (i = 0; i < argc; i++) {
    fprintf(stderr,"Processing file %s\n",argv[i]);
    read_table1 = Table_new(1000000,Table_string_compare,Table_string_hash);
    read_table2 = Table_new(1000000,Table_string_compare,Table_string_hash);

    bamreader = Bamread_new(argv[i]);
    parse_bam_input(unpaired_fp,splices_fp,read_output_fp,bamreader,read_table1,read_table2,chr_omit_p);
    Bamread_free(&bamreader);

    if ((n = Table_length(read_table2)) > 0) {
      fprintf(stderr,"Warning: for file %s, %d high-end reads were not completed or paired\n",argv[i],n);
      /* empty_table(stdout,read_table1,"/1"); */
    }
    if ((n = Table_length(read_table1)) > 0) {
      fprintf(stderr,"Warning: for file %s, %d low-end reads were not completed or paired\n",argv[i],n);
      /* empty_table(stdout,read_table1,"/1"); */
    }

    Table_free(&read_table2);
    Table_free(&read_table1);
  }

  fclose(splices_fp);
  fclose(unpaired_fp);
  fclose(read_output_fp);

  if (chr_omit_p != NULL) {
    FREE(chr_omit_p);
  }
  if (transcript_regions_iitfile != NULL) {
    IIT_free(&transcript_regions_iit);
  }

  return 0;
}

