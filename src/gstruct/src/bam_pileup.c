static char rcsid[] = "$Id: bam_pileup.c 132964 2014-04-10 20:25:16Z twu $";
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
#include <math.h>


#include "except.h"
#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "complement.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "genome.h"

#include "datadir.h"
#include "getopt.h"

/* Specific to BAM */
#include "samflags.h"		/* For flags */
#include "samread.h"
#include "bamread.h"
#include "parserange.h"


/* Alloc and block control structure */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Revise per read */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

/* Filters */
static int min_depth = 1;
static int variant_strands = 0;

/* Output options */
static bool blockp = true;
static bool print_cycles_p = false;
static bool print_quality_scores_p = false;
static bool want_genotypes_p = false;

static char *chromosome = NULL;

static bool print_chrpos_p = false;
static bool print_totals_p = false;

static bool signed_counts_p = false;

static int blocksize = 1000;



#if 0
static bool need_concordant_p = false;
#endif
static bool uniquep = false;

static int minimum_mapq = 0;

/* static int min_mlength = 0; */

#define DEFAULT_QUALITY 40  /* quality_score_adj + 40 */

static Genomicpos_T alloclength = 200000;

static int *block_max_column;
static int *alloc_max_column;

static char **block_pileup;
static char **alloc_pileup;

static List_T *block_start_accessions;
static List_T *block_end_accessions;

static List_T *alloc_start_accessions;
static List_T *alloc_end_accessions;



typedef enum {FIRST, SECOND} Pairend_T;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Filters */
  {"min-depth", required_argument, 0, 'n'}, /* min_depth */
  {"variants", required_argument, 0, 'X'}, /* variant_strands */

  /* Output options */
  {"block-format", required_argument, 0, 'B'}, /* blockp */
  {"signed-counts", no_argument, 0, 'S'}, /* signed_counts_p */
  {"positions", no_argument, 0, 'P'}, /* print_chrpos_p */
  {"totals", no_argument, 0, 'T'}, /* print_totals_p */
  {"cycles", no_argument, 0, 'C'}, /* print_cycles_p */
  {"quality-scores", no_argument, 0, 'Q'}, /* print_quality_scores_p */
  {"genotypes", no_argument, 0, 'G'}, /* want_genotypes_p */

  {"unique", required_argument, 0, 'U'}, /* uniquep */

  {"pairmax", required_argument, 0, 'p'}, /* alloclength */
  {"blocksize", required_argument, 0, 'b'}, /* blocksize */

  {"minimum-mapq", required_argument, 0, 'q'}, /* minimum_mapq */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_TALLY\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}



static void
print_program_usage ();


static void
accessions_gc (List_T *old) {
  List_T p;
  char *accession;

  for (p = *old; p != NULL; p = List_next(p)) {
    accession = (char *) List_head(p);
    FREE(accession);
  }
  List_free(&(*old));
  return;
}


static void
print_block (Genomicpos_T blockstart, Genomicpos_T blockptr, 
	     Genome_T genome, char *chr, Genomicpos_T chroffset) {
  Genomicpos_T chrpos;
  int blocki, lasti;
  bool foundp;
  int c;
  List_T p;

  lasti = blockptr - blockstart;
  debug1(printf("Printing blocki 0 to %d\n",lasti));

  /* Check for foundp */
  foundp = false;
  for (blocki = 0; blocki < lasti; blocki++) {
    if (block_max_column[blocki] >= 0) {
      foundp = true;
    }
  }

  
  /* Total could be 0 if the block is outside chrstart..chrend */
  if (foundp == true) {
    if (blockp == true) {
      printf(">0 %s:%u..%u\n",chr,blockstart,blockptr-1U);
    }

    for (blocki = 0; blocki < lasti; blocki++) {
      chrpos = blockstart + blocki;
      printf("%s\t%u\t",chr,chrpos);
      for (c = 0; c <= block_max_column[blocki]; c++) {
	if (block_pileup[blocki][c] == '\0') {
	  printf(" ");
	} else {
	  printf("%c",block_pileup[blocki][c]);
	}
      }

      printf("\t");
      for (p = block_start_accessions[blocki]; p != NULL; p = List_next(p)) {
	printf("%s ",(char *) List_head(p));
      }

      printf("\t");
      for (p = block_end_accessions[blocki]; p != NULL; p = List_next(p)) {
	printf("%s ",(char *) List_head(p));
      }

      printf("\n");

      /* Reset */
      block_max_column[blocki] = -1;
      FREE(block_pileup[blocki]);
      accessions_gc(&(block_start_accessions[blocki]));
      accessions_gc(&(block_end_accessions[blocki]));
    }

  }
    

  return;
}



static void
transfer_position (Genomicpos_T *blockptr, Genomicpos_T *blockstart, Genomicpos_T *blockend,
		   Genomicpos_T chrpos, int alloci,
		   Genome_T genome, char *chr, Genomicpos_T chroffset,
		   Genomicpos_T chrstart, Genomicpos_T chrend) {
  int blocki;

  if (chrpos < chrstart) {
    debug0(printf("    excluding chrpos %u < chrstart %u\n",chrpos,chrstart));
    alloc_max_column[alloci] = -1;
    alloc_pileup[alloci] = (char *) NULL;
    alloc_start_accessions[alloci] = (List_T) NULL;
    alloc_end_accessions[alloci] = (List_T) NULL;
  } else if (chrpos > chrend) {
    debug0(printf("    excluding chrpos %u > chrend %u\n",chrpos,chrend));
    alloc_max_column[alloci] = -1;
    alloc_pileup[alloci] = (char *) NULL;
    alloc_start_accessions[alloci] = (List_T) NULL;
    alloc_end_accessions[alloci] = (List_T) NULL;
  } else {
    if (chrpos >= *blockend) {
      debug0(printf("    chrpos %u >= blockend %u\n",chrpos,*blockend));
      debug0(printf("      print block from blockstart %u to blockptr %u\n",*blockstart,*blockptr));
      
      print_block(*blockstart,*blockptr,genome,chr,chroffset);
      debug0(printf("      set blockstart to chrpos, blockend to chrpos + blocksize\n"));
      *blockstart = chrpos;
      *blockend = chrpos + blocksize;
    }

    blocki = chrpos - (*blockstart);
#if 0
    debug0(printf("      transfer position from alloci %d to blocki %d\n",alloci,blocki));
#endif

    block_max_column[blocki] = alloc_max_column[alloci];
    alloc_max_column[alloci] = -1;
    
    block_pileup[blocki] = alloc_pileup[alloci];
    alloc_pileup[alloci] = (char *) NULL;

    block_start_accessions[blocki] = alloc_start_accessions[alloci];
    alloc_start_accessions[blocki] = (List_T) NULL;

    block_end_accessions[blocki] = alloc_end_accessions[alloci];
    alloc_end_accessions[blocki] = (List_T) NULL;

    /* Points just beyond maximum chrpos observed */
    if (chrpos + 1U > *blockptr) {
      *blockptr = chrpos + 1U;
#if 0
      debug0(printf("    advancing blockptr to %u\n",*blockptr));
#endif
    }
  }

  return;
}


static void
revise_read (Genomicpos_T chrpos_low, Genomicpos_T chrpos_high,
	     unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
	     char *acc, char *shortread, char *quality_string, Genomicpos_T alloc_low,
	     Genome_T genome, char *chr, Genomicpos_T chroffset,
	     Genomicpos_T chrstart, Genomicpos_T chrend) {
  bool availp[1000000];

  int alloci, firsti = -1, lasti = -1;
  char *copy;

  int shift;
  char strand;
  Genomicpos_T pos;
  char genomic[1024], *p, *q, *r;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;

  int c;

  debug1(printf("Revising read at %s:%u..%u\n",chr,chrpos_low,chrpos_high));

  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
  }


  /* Pass 1: Get best available column */
  for (c = 0; c < 1000000; c++) {
    availp[c] = true;
  }

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */
  p = shortread;
  r = quality_string;
  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
      /* r += mlength; -- hard clip means quality string is absent */
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'N') {
      pos += mlength;
    } else if (type == 'P') {
      /* Do nothing */
    } else if (type == 'I') {
      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'D') {
      pos += mlength;
    } else if (type == 'M') {
      if (0 /* mlength < min_mlength */) {
#if 0
	debug1(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
#endif
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	/* printf("After (+): %s\n",genomic); */
	q = genomic;

	/* psave = p; qsave = q; */
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));

	  /* signed_shift = (strand == '+') ? shift : -shift; */

	  if (firsti == -1) {
	    firsti = alloci;
	  }
	  lasti = alloci;

	  /******************************/
	  for (c = 0; c <= alloc_max_column[alloci]; c++) {
	    if (alloc_pileup[alloci][c] != '\0') {
	      availp[c] = false;
	    }
	  }
	  /******************************/

	  p++;
	  q++;
	  r++;
	  pos++;
	  shift += ((strand == '+') ? +1 : -1);
	}
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }
  
  c = 0;
  while (c < 1000000 && availp[c] == false) {
    c++;
  }
  if (c >= 1000000) {
    return;
  }


  if (firsti >= 0) {
    copy = (char *) CALLOC(strlen(acc)+1,sizeof(char));
    strcpy(copy,acc);
    alloc_start_accessions[firsti] = List_push(alloc_start_accessions[firsti],(void *) copy);
  }

  if (lasti >= 0) {
    copy = (char *) CALLOC(strlen(acc)+1,sizeof(char));
    strcpy(copy,acc);
    alloc_end_accessions[lasti] = List_push(alloc_end_accessions[lasti],(void *) copy);
  }


  /* Pass 2: Store in column */
  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */
  p = shortread;
  r = quality_string;

  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
      /* r += mlength; -- hard clip means quality string is absent */
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'P') {
      /* Do nothing */
    } else if (type == 'N') {
      pos += mlength;
    } else if (type == 'I') {
      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);
    } else if (type == 'D') {
      pos += mlength;
    } else if (type == 'M') {
      if (0 /* mlength < min_mlength */) {
#if 0
	debug1(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
#endif
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	/* printf("After (+): %s\n",genomic); */
	q = genomic;

	/* psave = p; qsave = q; */
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));

	  /* signed_shift = (strand == '+') ? shift : -shift; */

	  /******************************/
	  if (alloc_max_column[alloci] == -1) {
	    alloc_pileup[alloci] = (char *) CALLOC(1000000,sizeof(char));
	  }
	  alloc_pileup[alloci][c] = /*querynt*/ *p;
	  if (c > alloc_max_column[alloci]) {
	    alloc_max_column[alloci] = c;
	  }
	  /******************************/

	  p++;
	  q++;
	  r++;
	  pos++;
	  shift += ((strand == '+') ? +1 : -1);
	}
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }

  return;
}


static void
parse_bam (Bamreader_T bamreader, char *chromosome, Genome_T genome,
	   Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend) {
  char *acc, *chr, *mate_chr, *shortread, *quality_string, *read_group;
  unsigned int flag;
  Genomicpos_T alloc_ptr, alloc_low, alloc_high, chrpos_low, chrpos_high, mate_chrpos_low, chrpos;
  Genomicpos_T blockptr = 0U, blockstart = 0U, blockend = 0U;
  Intlist_T types;
  Uintlist_T npositions;
  int readlength, cigar_querylength, mapq;
  int delta, alloci;
  bool goodp;
  bool terminalp;

  alloc_ptr = 0U;
  alloc_low = 0U;
  alloc_high = alloc_low + alloclength;
  goodp = false;

  while (goodp == false && Bamread_next_line(bamreader,&acc,&flag,&mapq,&chr,&chrpos_low,
					     &mate_chr,&mate_chrpos_low,
					     &types,&npositions,&cigar_querylength,
					     &readlength,&shortread,&quality_string,&read_group,
					     &terminalp) > 0) {
    chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low);
    debug0(printf(">%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Samread_print_cigar(types,npositions));
    debug0(printf("\n"));


    if (mapq < minimum_mapq) {
      /* Skip */
    } else if (uniquep == true && (flag & NOT_PRIMARY)) {
      /* Skip */
    } else {
      alloc_low = chrpos_low;
      alloc_high = alloc_low + alloclength;
      alloc_ptr = chrpos_low;

      blockstart = chrpos_low;
      blockend = blockstart + blocksize;
      blockptr = chrpos_low;

      debug0(printf("    initialize alloc_low %u, alloc_high = %u, blockstart = %u, blockend = %u\n",
		   alloc_low,alloc_high,blockstart,blockend));

      if (chrpos_high > alloc_high) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		acc,chr,chrpos_low,chrpos_high,alloclength);
      } else {
	if (chrpos_high + 1U > alloc_ptr) {
	  alloc_ptr = chrpos_high + 1U;
	  debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	}

	revise_read(chrpos_low,chrpos_high,flag,types,npositions,cigar_querylength,
		    acc,shortread,quality_string,alloc_low,
		    genome,chr,chroffset,chrstart,chrend);

	goodp = true;
      }
    }

    Intlist_free(&types);
    Uintlist_free(&npositions);
    FREE(shortread);
  }

  while (Bamread_next_line(bamreader,&acc,&flag,&mapq,&chr,&chrpos_low,
			   &mate_chr,&mate_chrpos_low,&types,&npositions,&cigar_querylength,
			   &readlength,&shortread,&quality_string,&read_group,
			   &terminalp) > 0) {
    chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low);

    debug0(printf("*  alloc: %u..%u..%u  block: %u..%u..%u  ",
		  alloc_low,alloc_ptr,alloc_high,blockstart,blockptr,blockend));
    debug0(printf("%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Samread_print_cigar(types,npositions));
    debug0(printf("\n"));
    
    if (mapq < minimum_mapq) {
      /* Skip */
    } else if (uniquep == true && (flag & NOT_PRIMARY)) {
      /* Skip */
    } else {
      if (chrpos_low > alloc_ptr) {
	debug0(printf("Case 1: chrpos_low %u > alloc_ptr %u\n",chrpos_low,alloc_ptr));
	debug0(printf("    transfer from alloc_low %u to alloc_ptr %u\n",alloc_low,alloc_ptr));
	for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
	  alloci = chrpos - alloc_low;
	  if (alloc_max_column[alloci] >= 0) {
	    transfer_position(&blockptr,&blockstart,&blockend,chrpos,alloci,
			      genome,chr,chroffset,chrstart,chrend);
	  }
	}

	debug0(printf("    reset alloc_low = chrpos_low\n"));
	alloc_low = chrpos_low;
	alloc_high = alloc_low + alloclength;
	
      } else if (chrpos_low > alloc_low) {
	debug0(printf("Case 2: chrpos_low %u > alloc_low %u\n",chrpos_low,alloc_low));
	debug0(printf("    transfer from alloc_low %u up to chrpos_low %u\n",alloc_low,chrpos_low));
	for (chrpos = alloc_low; chrpos < chrpos_low; chrpos++) {
	  alloci = chrpos - alloc_low;
	  if (alloc_max_column[alloci] >= 0) {
	    transfer_position(&blockptr,&blockstart,&blockend,chrpos,alloci,
			      genome,chr,chroffset,chrstart,chrend);
	  }
	}

	delta = chrpos_low - alloc_low;
	debug0(printf("    shift alloc downward by %d so alloc_low = chrpos_low\n",delta));
	for (alloci = 0; alloci < alloc_ptr - alloc_low - delta; alloci++) {
	  alloc_max_column[alloci] = alloc_max_column[alloci+delta];
	  alloc_pileup[alloci] = alloc_pileup[alloci+delta];
	  alloc_start_accessions[alloci] = alloc_start_accessions[alloci+delta];
	  alloc_end_accessions[alloci] = alloc_end_accessions[alloci+delta];
	}
	for ( ; alloci < alloc_ptr - alloc_low; alloci++) {
	  debug1(printf("Resetting alloci %d\n",alloci));
	  alloc_max_column[alloci] = -1;
	  alloc_pileup[alloci] = (char *) NULL;
	  alloc_start_accessions[alloci] = (List_T) NULL;
	  alloc_end_accessions[alloci] = (List_T) NULL;
	}
	alloc_low = chrpos_low;
	alloc_high = alloc_low + alloclength;
	
      } else if (chrpos_low == alloc_low) {
	debug0(printf("Case 3: chrpos_low %u == alloc_low %u\n",chrpos_low,alloc_low));

      } else {
	fprintf(stderr,"Sequences are not in order.  Got chrpos_low %u < alloc_low %u\n",
		chrpos_low,alloc_low);
	abort();
      }

      if (chrpos_high > alloc_high) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		acc,chr,chrpos_low,chrpos_high,alloc_high);
      } else {
	if (chrpos_high + 1U > alloc_ptr) {
	  alloc_ptr = chrpos_high + 1U;
	  debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	}
	revise_read(chrpos_low,chrpos_high,flag,types,npositions,cigar_querylength,
		    acc,shortread,quality_string,alloc_low,
		    genome,chr,chroffset,chrstart,chrend);
      }
    }

    Intlist_free(&types);
    Uintlist_free(&npositions);
    FREE(shortread);
  }

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));
  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
    alloci = chrpos - alloc_low;
    if (alloc_max_column[alloci] >= 0) {
      transfer_position(&blockptr,&blockstart,&blockend,chrpos,alloci,
			genome,chr,chroffset,chrstart,chrend);
    }
  }

  debug0(printf("print block from blockstart %u to blockptr %u\n",blockstart,blockptr));
  print_block(blockstart,blockptr,genome,chr,chroffset);

  return;
}





int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL;
  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile;
  int index;
  IIT_T chromosome_iit;
  bool allocp;

  Genome_T genome = NULL;
  Bamreader_T bamreader;

  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;

  int i;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:B:p:n:X:SPTGCQq:U:b:",
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

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'B':
      if (!strcmp(optarg,"0")) {
	blockp = false;
      } else if (!strcmp(optarg,"1")) {
	blockp = true;
      } else {
	fprintf(stderr,"Argument to -B flag must be 0 or 1\n");
	exit(9);
      }
      break;

    case 'n': min_depth = atoi(optarg); break;
    case 'X': variant_strands = atoi(optarg); break;

    case 'p': alloclength = strtoul(optarg,NULL,10); break;

    case 'S': signed_counts_p = true; break;
    case 'P': print_chrpos_p = true; break;
    case 'T': print_totals_p = true; break;
    case 'G': want_genotypes_p = true; break;
    case 'C': print_cycles_p = true; break;
    case 'Q': print_quality_scores_p = true; break;

    case 'q': minimum_mapq = atoi(optarg); break;

#if 0
    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;
#endif

    case 'U':
      switch (atoi(optarg)) {
      case 0: uniquep = false; break;
      case 1: uniquep = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'b': blocksize = atoi(optarg); break;

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

  fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);

  alloc_max_column = (int *) CALLOC(alloclength,sizeof(int));
  alloc_pileup = (char **) CALLOC(alloclength,sizeof(char *));
  alloc_start_accessions = (List_T *) CALLOC(alloclength,sizeof(List_T));
  alloc_end_accessions = (List_T *) CALLOC(alloclength,sizeof(List_T));
  for (i = 0; i < alloclength; i++) {
    alloc_max_column[i] = -1;
  }

  block_max_column = (int *) CALLOC(blocksize,sizeof(int));
  block_pileup = (char **) CALLOC(blocksize,sizeof(char *));
  block_start_accessions = (List_T *) CALLOC(blocksize,sizeof(List_T));
  block_end_accessions = (List_T *) CALLOC(blocksize,sizeof(List_T));
  for (i = 0; i < blocksize; i++) {
    block_max_column[i] = -1;
  }


  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);

  bamreader = Bamread_new(argv[0]);
  if (argc > 1) {
    Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			 &chroffset,&chrlength,argv[1],genomesubdir,fileroot);
    Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
    parse_bam(bamreader,chromosome,genome,chroffset,chrstart,chrend);
    Bamread_unlimit_region(bamreader);

  } else {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chromosome = IIT_label(chromosome_iit,index,&allocp);
      chrstart = 1;
      chrend = Interval_length(IIT_interval(chromosome_iit,index));
      chroffset = Interval_low(IIT_interval(chromosome_iit,index));

      Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
      parse_bam(bamreader,chromosome,genome,chroffset,chrstart,chrend);
      Bamread_unlimit_region(bamreader);

      if (allocp == true) {
	FREE(chromosome);
      }
    }

    IIT_free(&chromosome_iit);
  }

  Bamread_free(&bamreader);

  if (genome != NULL) {
    Genome_free(&genome);
  }

  FREE(fileroot);
  FREE(genomesubdir);

  
  FREE(alloc_pileup);
  FREE(alloc_max_column);

  FREE(block_pileup);
  FREE(block_max_column);

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_pileup [OPTIONS...] bamfile chromosome:range, or\n\
       bam_pileup [OPTIONS...] bamfile\n\
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
Input options (must include -d)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
\n\
Compute options\n\
  -q, --min-mapq=INT             Use only alignments with this mapping quality and higher\n\
  -p, --pairmax=INT              Expected insert length (discards alignments longer than\n\
                                     this value) [default=200000]\n\
\n\
Filtering of output (options may be combined)\n\
  -n, --min-depth=INT            Print only positions with this depth or more\n\
  -X, --variants=INT             Print only positions showing a variant allele\n\
                                   Argument of 0 means all positions (default)\n\
                                   Argument of 1 means an observation of 1 strand is required\n\
                                   Argument of 2 means an observation of 2 strands is required\n\
\n\
Output options\n\
  -B, --block-format=INT         Print in block format (0=no, 1=yes (default))\n\
  -b, --blocksize=INT            Block size for printing in block format [default 1000]\n\
  -S, --signed-counts            Print signed allele counts (as plus_count|minus_count)\n\
  -P, --positions                Print chromosomal position on each line\n\
  -T, --totals                   Print total count (i.e., coverage or depth)\n\
  -G, --genotypes                Print genotype information (probability and likelihood)\n\
  -C, --cycles                   Include details about cycles\n\
  -Q, --quality-scores           Include details about quality scores\n\
\n\
");
    return;
}
