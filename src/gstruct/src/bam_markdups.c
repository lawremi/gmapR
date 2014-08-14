static char rcsid[] = "$Id: bam_markdups.c 108654 2013-09-19 23:11:00Z twu $";
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


static bool coordonlyp = false;
static bool inplacep = false;


static struct option long_options[] = {
  /* Input options */
  {"coordonly", no_argument, 0, 0}, /* coordonlyp */
  {"inplace", no_argument, 0, 0}, /* inplacep */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_MARKDUPS -- Marks duplicate reads in BAM files\n");
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
Usage: bam_markdups [options] [-o <output>] <bam file...>\n\
\n\
Options:\n\
");

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




/************************************************************************
 *   Read_T (simplified from bam_fastq)
 ************************************************************************/


#define T Read_T
typedef struct T *T;
struct T {
  int flags;
  char *queryseq;
  char *quality_string;
  Genomicpos_T chrpos_low;

  char querystart_strand;
  char *querystart_chr;
  Genomicpos_T querystart_chrpos;

  char queryend_strand;
  char *queryend_chr;
  Genomicpos_T queryend_chrpos;
};


static void
Read_free (T *old) {
  FREE((*old)->queryend_chr);
  FREE((*old)->querystart_chr);
  FREE((*old)->queryseq);
  if ((*old)->quality_string != NULL) {
    FREE((*old)->quality_string);
  }
  FREE(*old);
  return;
}


static int
Read_quality_score (T read) {
  int score = 0;
  int i;

  if (read->quality_string != NULL) {
    for (i = 0; i < (int) strlen(read->quality_string); i++) {
      score += read->quality_string[i];
    }
  }

  return score;
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


static void
Read_add_chrinfo_querystart (Read_T read, int flag, char *chr, Genomicpos_T chrpos_low, Genomicpos_T chrpos_high) {

  if (read->querystart_chr != NULL) {
    fprintf(stderr,"Read already has querystart chr %s\n",read->querystart_chr);
    FREE(read->querystart_chr);
  }

  if (chr == NULL) {
    read->querystart_chr = (char *) CALLOC(2,sizeof(char));
    strcpy(read->querystart_chr,"*");
    read->querystart_chrpos = 0;
  } else {
    read->querystart_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(read->querystart_chr,chr);
    if (flag & QUERY_MINUSP) {
      read->querystart_strand = '-';
      read->querystart_chrpos = chrpos_high;
    } else {
      read->querystart_strand = '+';
      read->querystart_chrpos = chrpos_low;
    }
  }

  return;
}


static void
Read_add_chrinfo_queryend (Read_T read, int flag, char *chr, Genomicpos_T chrpos_low, Genomicpos_T chrpos_high) {

  if (read->queryend_chr != NULL) {
    fprintf(stderr,"Read already has queryend chr %s\n",read->queryend_chr);
    FREE(read->queryend_chr);
  }

  if (chr == NULL) {
    read->queryend_chr = (char *) CALLOC(2,sizeof(char));
    strcpy(read->queryend_chr,"*");
    read->queryend_chrpos = 0;
  } else {
    read->queryend_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(read->queryend_chr,chr);
    if (flag & QUERY_MINUSP) {
      read->queryend_strand = '-';
      read->queryend_chrpos = chrpos_low;
    } else {
      read->queryend_strand = '+';
      read->queryend_chrpos = chrpos_high;
    }
  }

  return;
}


static T
Read_create_bam (bool *completep, T read, Bamline_T bamline) {
  unsigned int flag;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  Genomicpos_T chrpos;
  char *shortread, *quality_string, *p, *q;
  char strand;
  int type;
  int cigar_querylength, mlength, querypos, i;
  /* bool insidep; */

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
    read->chrpos_low = 0U;
    read->querystart_chr = (char *) NULL;
    read->queryend_chr = (char *) NULL;
    Read_add_chrinfo_querystart(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    Read_add_chrinfo_queryend(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    *completep = true;
    return read;
  }

  types = Bamline_cigar_types(bamline);
  npositions = Bamline_cigar_npositions(bamline);
  cigar_querylength = Bamline_cigar_querylength(bamline);

  if (read == NULL) {
    /* Initialize */
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
    read->chrpos_low = Bamline_chrpos_low(bamline);
    read->querystart_chr = (char *) NULL;
    read->queryend_chr = (char *) NULL;
  }

  if (flag & QUERY_MINUSP) {
    strand = '-';
    querypos = cigar_querylength;
    if (Intlist_head(types) != 'H') {
      Read_add_chrinfo_queryend(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    }
  } else {
    strand = '+';
    querypos = 0;
    if (Intlist_head(types) != 'H') {
      Read_add_chrinfo_querystart(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    }
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

  if (type != 'H') {
    if (flag & QUERY_MINUSP) {
      Read_add_chrinfo_querystart(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
    } else {
      Read_add_chrinfo_queryend(read,flag,Bamline_chr(bamline),Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
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


/************************************************************************
 *   Readpair_T
 ************************************************************************/


typedef struct Readpair_T *Readpair_T;
struct Readpair_T {
  char *acc;
  Read_T lowread;
  Read_T highread;
};


static Readpair_T
Readpair_new (char *acc, Read_T lowread, Read_T highread) {
  Readpair_T new = (Readpair_T) MALLOC(sizeof(*new));

  debug(printf("Creating read pair from %u..%u and %u..%u\n",
	       lowread->querystart_chrpos,lowread->queryend_chrpos,
	       highread->querystart_chrpos,highread->queryend_chrpos));
  new->acc = (char *) CALLOC(strlen(acc)+1,sizeof(char));
  strcpy(new->acc,acc);
  new->lowread = lowread;
  new->highread = highread;
  return new;
}

static void
Readpair_free (Readpair_T *old) {
  FREE((*old)->acc);
  Read_free(&(*old)->lowread);
  Read_free(&(*old)->highread);
  FREE(*old);
  return;
}

static int
Readpair_quality_score (Readpair_T this) {
  return Read_quality_score(this->lowread) + Read_quality_score(this->highread);
}


static int
sequence_cmp (char *queryseq1, char *queryseq2) {
  int i;

  for (i = 0; i < (int) strlen(queryseq1); i++) {
    if (toupper(queryseq1[i]) < toupper(queryseq2[i])) {
      return -1;
    } else if (toupper(queryseq1[i]) > toupper(queryseq2[i])) {
      return +1;
    }
  }

  return 0;
}


static int
Readpair_cmp_coordonly (const void *a, const void *b) {
  Readpair_T x = * (Readpair_T *) a;
  Readpair_T y = * (Readpair_T *) b;
  int cmp;

  debug(printf("Comparing %u..%u..%u..%u with %u..%u..%u..%u\n",
	       x->lowread->querystart_chrpos,x->lowread->queryend_chrpos,
	       x->highread->querystart_chrpos,x->highread->queryend_chrpos,
	       y->lowread->querystart_chrpos,y->lowread->queryend_chrpos,
	       y->highread->querystart_chrpos,y->highread->queryend_chrpos));

  if ((cmp = strcmp(x->lowread->querystart_chr,y->lowread->querystart_chr)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  } else if (x->lowread->querystart_chrpos < y->lowread->querystart_chrpos) {
    return -1;
  } else if (y->lowread->querystart_chrpos < x->lowread->querystart_chrpos) {
    return +1;

  } else if ((cmp = strcmp(x->lowread->queryend_chr,y->lowread->queryend_chr)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  } else if (x->lowread->queryend_chrpos < y->lowread->queryend_chrpos) {
    return -1;
  } else if (y->lowread->queryend_chrpos < x->lowread->queryend_chrpos) {
    return +1;

  } else if ((cmp = strcmp(x->highread->querystart_chr,y->highread->querystart_chr)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  } else if (x->highread->querystart_chrpos < y->highread->querystart_chrpos) {
    return -1;
  } else if (y->highread->querystart_chrpos < x->highread->querystart_chrpos) {
    return +1;

  } else if ((cmp = strcmp(x->highread->queryend_chr,y->highread->queryend_chr)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;

  } else if (x->highread->queryend_chrpos < y->highread->queryend_chrpos) {
    return -1;
  } else if (y->highread->queryend_chrpos < x->highread->queryend_chrpos) {
    return +1;

  } else {
    return 0;
  }
}


static int
Readpair_cmp_coord_seq (const void *a, const void *b) {
  Readpair_T x, y;
  int cmp;

  if ((cmp = Readpair_cmp_coordonly(a,b)) != 0) {
    return cmp;

  } else {
    x = * (Readpair_T *) a;
    y = * (Readpair_T *) b;

    if ((cmp = sequence_cmp(x->lowread->queryseq,y->lowread->queryseq)) != 0) {
      return cmp;

    } else if ((cmp = sequence_cmp(x->highread->queryseq,y->highread->queryseq)) != 0) {
      return cmp;

    } else {
      return 0;
    }
  }
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


#if 0
static Genomicpos_T
process_dups_single (Genomicpos_T last_chrpos, Read_T read) {
  List_T low_dups;

  if (last_chrpos == 0U) {
    /* Store read for future comparison */
    Uinttable_put(chrpos_low_table,read->chrpos_low,List_push(NULL,read));
    last_chrpos = read->chrpos_low;
  } else if (read->chrpos_low != last_chrpos) {
    /* Process stored reads */
    low_dups = (List_T) Uinttable_get(chrpos_low_table,read->chrpos_low);
    /* find_duplicates(low_dups); */

    /* Clear stored reads */
    /* Uinttable_clear(chrpos_low_table); */
    Uinttable_put(chrpos_low_table,read->chrpos_low,read);
    last_chrpos = read->chrpos_low;
  } else {
    low_dups = (List_T) Uinttable_get(chrpos_low_table,read->chrpos_low);
    Uinttable_put(chrpos_low_table,read->chrpos_low,List_push(low_dups,read));
  }

  return last_chrpos;
}
#endif


static int
count_duplicates_coordonly (Readpair_T *array, int n) {
  int ndups = 0;
  int i, j;

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && Readpair_cmp_coordonly(&(array[j]),&(array[i])) == 0) {
      j++;
    }
    ndups += (j - 1) - i;
    i = j;
  }

  return ndups;
}

static int
count_duplicates_coord_seq (Readpair_T *array, int n) {
  int ndups = 0;
  int i, j;

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && Readpair_cmp_coord_seq(&(array[j]),&(array[i])) == 0) {
      j++;
    }
    ndups += (j - 1) - i;
    i = j;
  }

  return ndups;
}




static void
find_duplicates_paired (int *ndups_by_coordonly, int *ndups_by_coord_seq,
			Table_T dup_table, List_T samecoords) {
  Readpair_T *array;
  int n, i, j, k, kbest;
  int bestscore, score;
  char *acc;

  array = (Readpair_T *) List_to_array_n(&n,samecoords);
  qsort(array,n,sizeof(Readpair_T),Readpair_cmp_coord_seq);
  *ndups_by_coordonly += count_duplicates_coordonly(array,n);
  *ndups_by_coord_seq += count_duplicates_coord_seq(array,n);

  i = 0;
  while (i < n) {
    j = i+1;
    if (coordonlyp == true) {
      while (j < n && Readpair_cmp_coordonly(&(array[j]),&(array[i])) == 0) {
	j++;
      }
    } else {
      while (j < n && Readpair_cmp_coord_seq(&(array[j]),&(array[i])) == 0) {
	j++;
      }
    }

    if (j - 1 > i) {
      kbest = i;
      bestscore = Readpair_quality_score(array[i]);
      for (k = i + 1; k < j; k++) {
	if ((score = Readpair_quality_score(array[k])) > bestscore) {
	  kbest = k;
	  bestscore = score;
	}
      }
      for (k = i; k < j; k++) {
	if (k != kbest) {
	  acc = (char *) CALLOC(strlen(array[k]->acc)+1,sizeof(char));
	  strcpy(acc,array[k]->acc);
	  if (Table_get(dup_table,acc) != NULL) {
	    fprintf(stderr,"%s already in dup_table\n",acc);
	    FREE(acc);
	  } else {
	    Table_put(dup_table,acc,(void *) 1);
	  }
	}
      }
      /* printf("\n"); */
    }
    i = j;
  }
  FREE(array);

  return;
}



/* Stored based on the low coordinate of the highread */
static List_T
process_samecoords_paired (int *ndups_by_coordonly, int *ndups_by_coord_seq,
			   Genomicpos_T *last_chrpos, List_T samecoords, Table_T dup_table,
			   char *acc, Read_T lowread, Read_T highread) {
  List_T p;
  Readpair_T readpair;

  if (*last_chrpos == 0U) {
    /* Store read for future comparison */
    debug(printf("  last_chrpos == 0, so storing read and setting last_chrpos = %u\n",highread->chrpos_low));
    /* Clear stored reads */
    for (p = samecoords; p != NULL; p = List_next(p)) {
      readpair = (Readpair_T) List_head(p);
      Readpair_free(&readpair);
    }
    List_free(&samecoords);

    samecoords = List_push(NULL,Readpair_new(acc,lowread,highread));
    *last_chrpos = highread->chrpos_low;

  } else if (highread->chrpos_low != *last_chrpos) {
    /* Process stored reads */
    debug(printf("  highread->chrpos_low %u != last_chrpos %u, so finding duplicates\n",highread->chrpos_low,*last_chrpos));
    find_duplicates_paired(&(*ndups_by_coordonly),&(*ndups_by_coord_seq),dup_table,samecoords);

    /* Clear stored reads */
    for (p = samecoords; p != NULL; p = List_next(p)) {
      readpair = (Readpair_T) List_head(p);
      Readpair_free(&readpair);
    }
    List_free(&samecoords);
    samecoords = List_push(NULL,Readpair_new(acc,lowread,highread));

    *last_chrpos = highread->chrpos_low;

  } else {
    debug(printf("  Have a possible duplicate at position %u, so pushing paired-end read\n",*last_chrpos));
    samecoords = List_push(samecoords,Readpair_new(acc,lowread,highread));
  }

  return samecoords;
}


static Table_T
parse_bam_input (int *nprimary, int *npaired, int *ndups_by_coordonly, int *ndups_by_coord_seq,
		 Bamreader_T bamreader, Table_T read_table1, Table_T read_table2) {
  Table_T dup_table;
  Bamline_T bamline;
  char *acc, *acc_copy;
  Read_T read1, read2;
  bool found1p, found2p, complete1p, complete2p;
  bool copyp;
  Readpair_T readpair;
  List_T samecoords = NULL, p;

  char *chr = NULL, *last_chr = NULL;
  Genomicpos_T last_chrpos = 0U;

  *nprimary = *npaired = *ndups_by_coordonly = *ndups_by_coord_seq = 0;
  dup_table = Table_new(1000000,Table_string_compare,Table_string_hash);

  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					 /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/true,
					 /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    *nprimary += 1;
    if (chr == NULL) {
      if (Bamline_chr(bamline) == NULL) {
	chr = (char *) CALLOC(2,sizeof(char));
	strcpy(chr,"*");
      } else {
	chr = (char *) CALLOC(strlen(Bamline_chr(bamline)) + 1,sizeof(char));
	strcpy(chr,Bamline_chr(bamline));
      }
    }

    if (last_chr == NULL) {
      last_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
      strcpy(last_chr,chr);

    } else if (strcmp(chr,last_chr) != 0) {
      /* Clear stored reads */
      for (p = samecoords; p != NULL; p = List_next(p)) {
	readpair = (Readpair_T) List_head(p);
	Readpair_free(&readpair);
      }
      List_free(&samecoords);

      FREE(last_chr);
      last_chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
      strcpy(last_chr,chr);
      last_chrpos = 0U;

      FREE(chr);
      chr = NULL;
    }

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
	/* Have a complete read */
	/* last_chrpos = process_dups_single(last_chrpos,read1); */

	if (found1p == true) {
	  acc_copy = Table_remove(read_table1,acc);
	  FREE(acc_copy);
	}
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
	/* Have complete paired-end read */
	debug(printf("  Have a complete paired-end read\n"));
	*npaired += 1;
	samecoords = process_samecoords_paired(&(*ndups_by_coordonly),&(*ndups_by_coord_seq),
					       &last_chrpos,samecoords,dup_table,acc,
					       /*lowread*/read2,/*highread*/read1);

	if (found1p == true) {
	  acc_copy = Table_remove(read_table1,acc);
	  FREE(acc_copy);
	}
	acc_copy = Table_remove(read_table2,acc);
	FREE(acc_copy);
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
	/* Have complete paired-end read */
	debug(printf("  Have a complete paired-end read\n"));
	*npaired += 1;
	samecoords = process_samecoords_paired(&(*ndups_by_coordonly),&(*ndups_by_coord_seq),
					       &last_chrpos,samecoords,dup_table,acc,
					       /*lowread*/read2,/*highread*/read1);

	if (found2p == true) {
	  acc_copy = Table_remove(read_table2,acc);
	  FREE(acc_copy);
	}
	acc_copy = Table_remove(read_table1,acc);
	FREE(acc_copy);
      }
    }
    
    if (copyp) {
      FREE(acc);
    }

    Bamline_free(&bamline);
  }

  /* Clear stored reads */
  for (p = samecoords; p != NULL; p = List_next(p)) {
    readpair = (Readpair_T) List_head(p);
    Readpair_free(&readpair);
  }
  List_free(&samecoords);

  if (last_chr != NULL) {
    FREE(last_chr);
  }

  FREE(chr);

  return dup_table;
}


static void
mark_dups (int *nsam, int *nmarked, Bamreader_T bamreader, Table_T dup_table) {
  Bamline_T bamline;
  char *acc;
  bool copyp;

  *nsam = *nmarked = 0;
  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					 /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/false,
					 /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    (*nsam) += 1;
    acc = Bamline_acc(bamline);
    acc = strip_illumina_acc_ending(&copyp,acc);
    if (Table_get(dup_table,acc) != NULL) {
      (*nmarked) += 1;
      Bamline_print(stdout,bamline,/*newflag*/Bamline_flag(bamline) | DUPLICATE_READ,/*quality_score_adj*/0);
    } else {
      Bamline_print(stdout,bamline,/*newflag*/Bamline_flag(bamline),/*quality_score_adj*/0);
    }
    if (copyp) {
      FREE(acc);
    }
    Bamline_free(&bamline);
  }

  return;
}


#if 0
static void
modify_dups (Bamreader_T bamreader, Table_T dup_table) {
  Bamline_T bamline;
  char *acc;
  bool copyp;

  while ((bamline = Bamread_next_bamline(bamreader,/*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					 /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/true,
					 /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    acc = Bamline_acc(bamline);
    acc = strip_illumina_acc_ending(&copyp,acc);
    if (Table_get(dup_table,acc) != NULL) {
      Bamread_modify_flag(bamreader,/*newflag*/Bamline_flag(bamline) | DUPLICATE_READ);
    }
    if (copyp) {
      FREE(acc);
    }
    Bamline_free(&bamline);
  }

  return;
}
#endif



static void
empty_read_table (Table_T read_table) {
  char *acc, *acc_copy;
  char **accessions;
  Read_T read;
  int n, i;

  n = Table_length(read_table);
  accessions = (char **) Table_keys(read_table,(void *) NULL);
  for (i = 0; i < n; i++) {
    acc = accessions[i];
    read = (T) Table_get(read_table,acc);
    acc_copy = Table_remove(read_table,acc);
    FREE(acc_copy);
    Read_free(&read);
  }
  FREE(accessions);

  return;
}


static void
empty_dup_table (Table_T dup_table) {
  char *acc;
  char **accessions;
  int n, i;

  n = Table_length(dup_table);
  accessions = (char **) Table_keys(dup_table,(void *) NULL);
  for (i = 0; i < n; i++) {
    acc = accessions[i];
    acc = Table_remove(dup_table,acc);
    FREE(acc);
  }
  FREE(accessions);

  return;
}



int
main (int argc, char *argv[]) {
  Bamreader_T bamreader;
  Table_T dup_table, read_table1, read_table2;
  int nsam, nmarked, nprimary, npaired, ndups_by_coordonly, ndups_by_coord_seq;
  int n, i;

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
      if (!strcmp(long_name,"coordonly")) {
	coordonlyp = true;
      } else if (!strcmp(long_name,"inplace")) {
	inplacep = true;
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

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  for (i = 0; i < argc; i++) {
    fprintf(stderr,"Processing file %s\n",argv[i]);
    read_table1 = Table_new(1000000,Table_string_compare,Table_string_hash);
    read_table2 = Table_new(1000000,Table_string_compare,Table_string_hash);

    bamreader = Bamread_new(argv[i]);
    dup_table = parse_bam_input(&nprimary,&npaired,&ndups_by_coordonly,&ndups_by_coord_seq,bamreader,read_table1,read_table2);
    Bamread_free(&bamreader);

    if ((n = Table_length(read_table2)) > 0) {
      fprintf(stderr,"Warning: for file %s, %d high-end reads were not completed or paired\n",argv[i],n);
      empty_read_table(read_table2);
    }
    if ((n = Table_length(read_table1)) > 0) {
      fprintf(stderr,"Warning: for file %s, %d low-end reads were not completed or paired\n",argv[i],n);
      empty_read_table(read_table1);
    }

    Table_free(&read_table2);
    Table_free(&read_table1);

    /* Process BAM file for dups */
    if (inplacep == true) {
      fprintf(stderr,"Not able to modify BAM in place\n");
      exit(9);
#if 0
      bamreader = Bamread_modify(argv[i]);
      modify_dups(bamreader,dup_table);
      Bamread_free(&bamreader);
#endif
    } else {
      bamreader = Bamread_new(argv[i]);
      Bamread_write_header(bamreader);
      mark_dups(&nsam,&nmarked,bamreader,dup_table);
      Bamread_free(&bamreader);
    }

    fprintf(stderr,"For file %s, %d SAM lines, %d primary, %d (%.1f%%) marked as dups; ",
	    argv[i],nsam,nprimary,nmarked,100.0*(double) nmarked/(double) nsam);
    fprintf(stderr,"  %d accessions, %d dups (%.1f%%) by coord only, and %d (%.1f%%) by coord/sequence\n",
	      npaired,ndups_by_coordonly,100.0*(double) ndups_by_coordonly/(double) npaired,
	      ndups_by_coord_seq,100.0*(double) ndups_by_coord_seq/(double) npaired);

    empty_dup_table(dup_table);
    Table_free(&dup_table);
  }

  return 0;
}

