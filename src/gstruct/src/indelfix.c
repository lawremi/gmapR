static char rcsid[] = "$Id: indelfix.c 143438 2014-08-05 21:46:08Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "indelfix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mem.h"
#include "assert.h"
#include "genomicpos.h"

#include "samflags.h"
#include "bamread.h"
#include "samread.h"

#include "dynprog.h"
#include "dynprog_single.h"


#define ARRAY_THRESHOLD 20
#define DEFAULT_QUALITY 40  /* quality_score_adj + 40 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* tally_reads_multiple */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* realign_reads_two */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* realign_reads_multiple */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* read_nmismatches */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif


/* Processing neighborhoods */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


/* Parsing variants */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



typedef struct Variant_T *Variant_T;
struct Variant_T {
  int count;			/* When set to 0, indicates that it has been processed already */
  int delta;
  Intlist_T types;
  Uintlist_T npositions;
  Uintlist_T chrpos;
  Genomicpos_T first_chrpos;
  Genomicpos_T last_chrpos;
  List_T segments;
};


static void
Variant_print (Variant_T this) {
  Uintlist_T p, q;
  Intlist_T t;

  printf("Variant at %u..%u (count %d):",this->first_chrpos,this->last_chrpos,this->count);
  for (t = this->types, q = this->npositions, p = this->chrpos; t != NULL;
       t = Intlist_next(t), q = Uintlist_next(q), p = Uintlist_next(p)) {
    printf(" %u%c at %u",Uintlist_head(q),Intlist_head(t),Uintlist_head(p));
  }
  printf("\n");

  return;
}


static Variant_T
Variant_new (Intlist_T types, Uintlist_T npositions, Uintlist_T chrpos, List_T segments) {
  Variant_T new = (Variant_T) MALLOC(sizeof(*new));
  Genomicpos_T position;
  Uintlist_T p, q;

  new->count = 1;
  new->types = types;
  new->npositions = npositions;
  new->chrpos = chrpos;
  new->segments = segments;

  p = chrpos;
  q = npositions;
  new->first_chrpos = Uintlist_head(p) - Uintlist_head(q);
  new->last_chrpos = Uintlist_head(p) + Uintlist_head(q);
  for (p = Uintlist_next(p), q = Uintlist_next(q); p != NULL; p = Uintlist_next(p), q = Uintlist_next(q)) {
    if ((position = Uintlist_head(p) - Uintlist_head(q)) < new->first_chrpos) {
      new->first_chrpos = position;
    }
    if ((position = Uintlist_head(p) + Uintlist_head(q)) > new->last_chrpos) {
      new->last_chrpos = position;
    }
  }

  return new;
}

static void
Variant_free (Variant_T *old) {
  List_T p;
  char *segment;

  Intlist_free(&(*old)->types);
  Uintlist_free(&(*old)->npositions);
  Uintlist_free(&(*old)->chrpos);

  for (p = (*old)->segments; p != NULL; p = List_next(p)) {
    segment = (char *) List_head(p);
    FREE(segment);
  }
  List_free(&(*old)->segments);

  FREE(*old);
  return;
}


static int
Variant_count_cmp (const void *a, const void *b) {
  Variant_T x = * (Variant_T *) a;
  Variant_T y = * (Variant_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (y->count > x->count) {
    return +1;
  } else if (x->first_chrpos < y->first_chrpos) {
    return -1;
  } else if (y->first_chrpos < x->first_chrpos) {
    return +1;
  } else if (x->last_chrpos < y->last_chrpos) {
    return -1;
  } else if (y->last_chrpos < x->last_chrpos) {
    return +1;
  } else {
    return 0;
  }
}



static Variant_T
Variant_find (List_T variants, Intlist_T types, Uintlist_T npositions, Uintlist_T chrpos, List_T segments) {
  List_T p;
  Variant_T var;

  for (p = variants; p != NULL; p = List_next(p)) {
    var = (Variant_T) List_head(p);
    if (Intlist_equal(types,var->types) == true &&
	Uintlist_equal(npositions,var->npositions) == true &&
	Uintlist_equal(chrpos,var->chrpos) == true &&
	List_equal_strings(segments,var->segments) == true) {
      return var;
    }
  }

  return (Variant_T) NULL;
}


#if 0
static int **
Matrix_new (int length1, int length2) {
  int **matrix;
  int i;

  matrix = (int **) MALLOC(length2 * sizeof(int *));
  matrix[0] = (int *) CALLOC(length1*length2,sizeof(int));
  for (i = 1; i < length2; i++) {
    matrix[i] = &(matrix[i-1][length1]);
  }

  return matrix;
}
#endif


static void
find_segmentation (int *equiv_x1, int *equiv_x2, int *data, int n) {
  int *cumdata;
  int x1, x2, i;
  double lowest_chisq, chisq;
  bool firstp = true;
  int sum_inside, sum_outside;
  int n_inside, n_outside;
  double mean_inside, mean_outside, delta;


  cumdata = (int *) CALLOC(n+1,sizeof(int));
  cumdata[0] = 0;
  for (i = 1; i <= n; i++) {
    cumdata[i] = cumdata[i-1] + data[i-1];
  }

  /* equivalence class is [x1, x2] inclusive */
  for (x1 = 0; x1 < n - 1; x1++) {
    for (x2 = x1; x2 < n; x2++) {
      sum_inside = cumdata[x2 + 1] - cumdata[x1];
      sum_outside = cumdata[n] - sum_inside;
      n_inside = (x2 + 1) - x1;
      n_outside = n - n_inside;
      mean_inside = (double) (sum_inside + 0.5) / (double) (n_inside + 0.5);
      mean_outside = (double) (sum_outside + 0.5) / (double) (n_outside + 0.5);
      debug8(printf("%d %d: inside %d (%d) => %f, outside %d (%d) => %f",
		    x1,x2,sum_inside,n_inside,mean_inside,sum_outside,n_outside,mean_outside));

      chisq = 0.0;
      i = 0;
      while (i < x1) {
	delta = (data[i] - mean_outside);
	chisq += (delta*delta)/mean_outside;
	i++;
      }
      while (i <= x2) {
	delta = (data[i] - mean_inside);
	chisq += (delta*delta)/mean_inside;
	i++;
      }
      while (i < n) {
	delta = (data[i] - mean_outside);
	chisq += (delta*delta)/mean_outside;
	i++;
      }
      debug8(printf("  chisq = %.1f",chisq));
      
      if (firstp == true) {
	lowest_chisq = chisq;
	firstp = false;
	*equiv_x1 = x1;
	*equiv_x2 = x2;
	debug8(printf("**"));
      } else if (chisq < lowest_chisq) {
	lowest_chisq = chisq;
	*equiv_x1 = x1;
	*equiv_x2 = x2;
	debug8(printf("**"));
      }

      debug8(printf("\n"));
    }
  }

  FREE(cumdata);

  return;
}



static int
read_nmismatches (Genomicpos_T chrpos_low, Intlist_T types, Uintlist_T npositions,
		  char *shortread, Genome_T genome, Genomicpos_T chroffset) {
  int nmismatches = 0;
  Genomicpos_T pos;
  char *genomic = NULL, *p, *q;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;


  debug5(printf("Count read mismatches at %u\n",chrpos_low));

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */

  p = shortread;

  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;

    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
      /* r += mlength; -- hard clip means quality string is absent */

    } else if (type == 'N') {
      pos += mlength;

    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'I') {
      p += mlength;

    } else if (type == 'D') {
      pos += mlength;

    } else if (type == 'M') {
      if (0 /* mlength < min_mlength */) {
	p += mlength;
	pos += mlength;

      } else {
	debug5(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  q = &(*p);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	/* psave = p; qsave = q; */
	debug5(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  debug5(printf("Processing %c and %c at pos %u, mlength %u\n",*p,*q,pos+1U,mlength));
	  if (*p != *q) {
	    nmismatches += 1;
	  }

	  p++;
	  q++;
	  pos++;
	}
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }
  
  FREE(genomic);
  return nmismatches;
}




static int
find_max_readlength (Bamreader_T realign_bamreader, char *chr, Genomicpos_T indelpos_start, Genomicpos_T indelpos_end) {
  int max_readlength = 0, readlength;
  Bamline_T bamline;

  Bamread_limit_region(realign_bamreader,chr,indelpos_start,indelpos_end);
  while ((bamline = Bamread_next_bamline(realign_bamreader,/*desired_read_group*/NULL,
					 /*minimum_mapq*/0,/*good_unique_mapq*/0,/*maximum_nhits*/1000000,
					 /*need_unique_p*/false,/*need_primary_p*/false,/*ignore_duplicates_p*/false,
					 /*need_concordant_p*/false)) != NULL) {
    if ((readlength = Bamline_readlength(bamline)) > max_readlength) {
      max_readlength = readlength;
    }
    Bamline_free(&bamline);
  }
   
  Bamread_unlimit_region(realign_bamreader);
  return max_readlength;
}


/* Selects among multiple candidates */
static int
tally_reads_multiple (Bamreader_T realign_bamreader, char *chr, Genomicpos_T indelpos_start, Genomicpos_T indelpos_end,
		      Genomicpos_T chrpos_start, Genomicpos_T chrpos_end,
		      char *genomicseg0, char **genomicsegs_alt, int *glengths_alt, Genomicpos_T *initpos_alt,
		      int *nindels0, int **nindels_alt, char ***delstrings_alt,
		      Variant_T *var_array, int nvariants, int keyi, int max_varlength, Dynprog_T dynprog) {
  int *support;
  int *equivalent;
  int equiv_x1, equiv_x2;
  Bamline_T *bamlines, bamline, bamline_rep, prev_bamline = NULL;
  int nlines, linei_start, linei_end, linei, nreps;
  char *rsequence;
  int readlength;

  int glength0;
  int best_finalc, finalc0, finalc1;
  int besti, best_score, score0, *scores;
  char *best_cigar, *orig_cigar, *cigar0, *cigar1;
  char *best_md_string, *md_string_0, *md_string_1;
  int i, j;

  int best_nmismatches, orig_nmismatches, nmismatches0, nmismatches1;

  Genomicpos_T readalign_start, readalign_end, total_ins;
  int readalign_length;

	       
  support = (int *) CALLOC(nvariants,sizeof(int));
  scores = (int *) CALLOC(nvariants,sizeof(int));
  equivalent = (int *) CALLOC(nvariants,sizeof(int));

  glength0 = chrpos_end - chrpos_start + 1;
  Bamread_limit_region(realign_bamreader,chr,indelpos_start,indelpos_end);

  while ((bamlines = Bamread_next_bamline_set(&nlines,&prev_bamline,realign_bamreader,
					      /*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					      /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/false,
					      /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {
    fprintf(stderr,"Have %d lines at chrpos %u\n",nlines,Bamline_chrpos_low(prev_bamline));

    for (linei_start = 0; linei_start < nlines; linei_start = linei_end) {
      linei_end = linei_start + 1;
      while (linei_end < nlines && !strcmp(Bamline_read(bamlines[linei_start]),Bamline_read(bamlines[linei_end]))) {
	linei_end++;
      }

      bamline_rep = bamlines[linei_start];
      nreps = linei_end - linei_start;
      orig_nmismatches = Bamline_nmismatches(bamline_rep);
      
      debug1(printf("\n"));
      debug1(printf(">%s:%u..%u %s ",chr,Bamline_chrpos_low(bamline_rep),Bamline_chrpos_high(bamline_rep),Bamline_acc(bamline_rep)));
      debug1(Bamread_print_cigar(stdout,bamline_rep));
      debug1(printf(" %d mismatches\n",orig_nmismatches));

      rsequence = Bamline_read(bamline_rep);
      readlength = Bamline_readlength(bamline_rep);
      debug1(printf("%s\n",rsequence));

#ifdef WHOLE_ALIGNMENT
      readalign_start = chrpos_start;
      readalign_length = glength0;
#else
      total_ins = Bamline_total_ins(bamline_rep);
      if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	readalign_start = chrpos_start;
      }
      if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_end) {
	readalign_end = chrpos_end;
      }
      readalign_length = readalign_end - readalign_start + 1;
#endif

      /* Use nogap here */
      cigar0 = Dynprog_single_nogap(&score0,&nmismatches0,&finalc0,&md_string_0,dynprog,rsequence,
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(nindels0[readalign_start - chrpos_start]),/*deletion_string*/NULL,
				    /*rlength*/readlength,/*glength*/readalign_length,
				    /*extraband_single*/total_ins + max_varlength);

      debug1(printf(" wildtype:%d with %d mismatches at %d (%u) %s",
		    score0,nmismatches0,finalc0,readalign_start+finalc0,cigar0));

      besti = -1;
      best_score = score0;
      best_nmismatches = nmismatches0;
      best_finalc = finalc0;
      best_cigar = cigar0;
      best_md_string = md_string_0;
      orig_cigar = Bamline_cigar_string(bamline_rep);

      for (i = 0; i < nvariants; i++) {
	debug1(printf("  variant-%d[%u..%u]:",i,var_array[i]->first_chrpos,var_array[i]->last_chrpos));

#ifdef WHOLE_ALIGNMENT
	readalign_start = chrpos_start;
	readalign_length = glengths_alt[i];
#else
	/* if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	   readalign_start = chrpos_start;
	   } */
	if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_start + glengths_alt[i]) {
	  readalign_end = chrpos_start + glengths_alt[i];
	}
#endif

	if (var_array[i]->last_chrpos < readalign_start) {
	  /* Skip, because variants are all to the left of the alignment */
	  debug1(printf("all-left"));
	  scores[i] = -1;
	} else if (var_array[i]->first_chrpos > readalign_end) {
	  /* Skip, because variants are all to the right of the alignment */
	  debug1(printf("all-right"));
	  scores[i] = -1;
	} else {
	  /* Align */
	  readalign_length = readalign_end - readalign_start + 1;

	  /* Use nogap here */
	  cigar1 = Dynprog_single_nogap(&(scores[i]),&nmismatches1,&finalc1,&md_string_1,dynprog,rsequence,
					&(genomicsegs_alt[i][readalign_start - chrpos_start]),
					&(genomicsegs_alt[i][readalign_start - chrpos_start]),
					&(nindels_alt[i][readalign_start - chrpos_start]),
					&(delstrings_alt[i][readalign_start - chrpos_start]),
					/*rlength*/readlength,/*glength*/readalign_length,
					/*extraband_single*/total_ins + max_varlength);
	  debug1(printf("%d with %d mismatches at %d (%u) %s",scores[i],nmismatches1,finalc1,readalign_start+finalc1,cigar1));


	  if (scores[i] < best_score) {
	    FREE(cigar1);
	    FREE(md_string_1);

	  } else if (scores[i] == best_score) {
	    FREE(cigar1);
	    FREE(md_string_1);
	    besti = -1;

	  } else {
	    debug1(printf("** "));
	    FREE(best_cigar);
	    FREE(best_md_string);

	    best_finalc = finalc1;
	    best_cigar = cigar1;
	    best_md_string = md_string_1;
	    best_score = scores[i];
	    best_nmismatches = nmismatches1;
	    besti = i;
	  }
	}
      }
      debug1(printf("\n"));

      if (best_nmismatches <= orig_nmismatches + 1) {
	if (scores[keyi] == best_score) {
	  for (j = 0; j < nvariants; j++) {
	    if (scores[j] < 0) {
	      /* Skip since not evaluated */
	    } else if (scores[j] == best_score) {
	      equivalent[j] += 1;
	    } else if (scores[j] < best_score) {
	      /* equivalent means equal or dominated by key_var */
	      equivalent[j] += 1;
	    }
	  }
	}

	for (i = 0; i < nvariants; i++) {
	  if (scores[i] < 0) {
	    /* Skip, since not evaluated */
	  } else if (scores[i] == best_score) {
	    support[i] += nreps;
	  }
	}
	  
	debug1(printf("Support for candidate %d is %d\n",0,support[0]));
	for (i = 1; i < nvariants; i++) {
	  debug1(printf("Support for candidate %d is %d\n",i,support[i]));
	}
      }
      
      FREE(orig_cigar);
      FREE(best_cigar);
      FREE(best_md_string);
    }

    for (linei = 0; linei < nlines; linei++) {
      bamline = bamlines[linei];
      Bamline_free(&bamline);
    }
    FREE(bamlines);
  }

  Bamread_unlimit_region(realign_bamreader);

#ifdef DEBUG8
  printf("Equivalence:");
  for (j = 0; j < nvariants; j++) {
    printf(" %06d",equivalent[j]);
  }
  printf("\n");
#endif

  find_segmentation(&equiv_x1,&equiv_x2,equivalent,nvariants);
  

  FREE(equivalent);

#ifdef DEBUG8
  for (i = 0; i < nvariants; i++) {
    printf("Support for candidate %d (initpos %u) is %d",i,initpos_alt[i],support[i]);
    if (i >= equiv_x1 && i <= equiv_x2) {
      printf(" **");
    }
    printf("\n");
  }
#endif

  besti = equiv_x1;
  for (i = equiv_x1; i <= equiv_x2; i++) {
    if (support[i] > support[besti]) {
      besti = i;
    } else if (support[i] == support[besti] && initpos_alt[i] < initpos_alt[besti]) {
      besti = i;
    }
    var_array[i]->count = 0;	/* Mark so it is not evaluated again */
  }
    
  FREE(scores);
  FREE(support);

  debug1(printf("Candidate in equivalence class with best support is %d\n",besti));
  return besti;
}


/* Allows for multiple candidates */
static void
realign_reads_multiple (Bamreader_T realign_bamreader, char *chr, Genomicpos_T indelpos_start, Genomicpos_T indelpos_end,
			Genomicpos_T chrpos_start, Genomicpos_T chrpos_end, Genome_T genome, Genomicpos_T chroffset,
			char *genomicseg0, char **genomicsegs_alt, int *glengths_alt,
			int *nindels0, int **nindels_alt, char ***delstrings_alt, Genomicpos_T **wildtype_chrpos,
			int nvariants, int max_varlength, Dynprog_T dynprog) {
  Bamline_T bamline, bamline_rep, *bamlines, prev_bamline = NULL;
  int nlines, linei_start, linei_end, linei;
  char *rsequence;
  int readlength, new_readlength;
  int orig_nmismatches, new_nmismatches;
  Intlist_T new_cigar_types;
  Uintlist_T new_cigar_npositions;

  int glength0;
  Genomicpos_T best_chrpos;
  int finalc0, finalc1;
  int best_score, score0, score1;
  char *best_cigar, *orig_cigar, *cigar0, *cigar1;
  char *best_md_string, *md_string_0, *md_string_1;
  int i;

  Genomicpos_T readalign_start, readalign_end, total_ins;
  int readalign_length;

	       
  glength0 = chrpos_end - chrpos_start + 1;
  Bamread_limit_region(realign_bamreader,chr,indelpos_start,indelpos_end);

  while ((bamlines = Bamread_next_bamline_set(&nlines,&prev_bamline,realign_bamreader,
					      /*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					      /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/false,
					      /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {

    for (linei_start = 0; linei_start < nlines; linei_start = linei_end) {
      linei_end = linei_start + 1;
      while (linei_end < nlines && !strcmp(Bamline_read(bamlines[linei_start]),Bamline_read(bamlines[linei_end]))) {
	linei_end++;
      }

      bamline_rep = bamlines[linei_start];

      debug3(printf("\n"));
      debug3(printf(">%s:%u..%u ",chr,Bamline_chrpos_low(bamline_rep),Bamline_chrpos_high(bamline_rep)));
      debug3(Bamread_print_cigar(stdout,bamline_rep));

      rsequence = Bamline_read(bamline_rep);
      readlength = Bamline_readlength(bamline_rep);

#ifdef WHOLE_ALIGNMENT
      readalign_start = chrpos_start;
      readalign_length = glength0;
#else
      total_ins = Bamline_total_ins(bamline_rep);
      if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	readalign_start = chrpos_start;
      }
      if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_end) {
	readalign_end = chrpos_end;
      }
      readalign_length = readalign_end - readalign_start + 1;
#endif
    

#ifdef NOGAP
      cigar0 = Dynprog_single_nogap(&score0,&finalc0,&md_string_0,dynprog,rsequence,
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(nindels0[readalign_start - chrpos_start]),/*deletion_string*/NULL,
				    /*rlength*/readlength,/*glength*/readalign_length,
				    /*extraband_single*/total_ins + max_varlength);
#else
      /* Use gap here to account for sequence errors */
      cigar0 = Dynprog_single_gap(&score0,&finalc0,&md_string_0,dynprog,rsequence,
				  &(genomicseg0[readalign_start - chrpos_start]),
				  &(genomicseg0[readalign_start - chrpos_start]),
				  &(nindels0[readalign_start - chrpos_start]),/*deletion_string*/NULL,
				  /*rlength*/readlength,/*glength*/readalign_length,
				  /*extraband_single*/total_ins + max_varlength);
#endif

      debug3(printf(" wildtype:%d at %d (%u) %s",score0,finalc0,readalign_start+finalc0,cigar0));

      best_score = score0;
      best_chrpos = readalign_start + finalc0;
      best_cigar = cigar0;
      best_md_string = md_string_0;

      for (i = 0; i < nvariants; i++) {
#ifdef WHOLE_ALIGNMENT
	readalign_start = chrpos_start;
	readalign_length = glengths_alt[i];
#else
	/* if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	   readalign_start = chrpos_start;
	   } */
	if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_start + glengths_alt[i]) {
	  readalign_end = chrpos_start + glengths_alt[i];
	}
	readalign_length = readalign_end - readalign_start + 1;
#endif

#ifdef NOGAP
	cigar1 = Dynprog_single_nogap(&score1,&finalc1,&md_string_1,dynprog,rsequence,
				      &(genomicsegs_alt[i][readalign_start - chrpos_start]),
				      &(genomicsegs_alt[i][readalign_start - chrpos_start]),
				      &(nindels_alt[i][readalign_start - chrpos_start]),
				      &(delstrings_alt[i][readalign_start - chrpos_start]),
				      /*rlength*/readlength,/*glength*/readalign_length,
				      /*extraband_single*/total_ins + max_varlength);
#else
	/* Use gap here to account for sequence errors */
	cigar1 = Dynprog_single_gap(&score1,&finalc1,&md_string_1,dynprog,rsequence,
				    &(genomicsegs_alt[i][readalign_start - chrpos_start]),
				    &(genomicsegs_alt[i][readalign_start - chrpos_start]),
				    &(nindels_alt[i][readalign_start - chrpos_start]),
				    &(delstrings_alt[i][readalign_start - chrpos_start]),
				    /*rlength*/readlength,/*glength*/readalign_length,
				    /*extraband_single*/total_ins + max_varlength);
#endif
	debug3(printf("  variant-%d:%d at %d (%u) %s",i,score1,finalc1,readalign_start+finalc1,cigar1));
	if (score1 <= best_score) {
	  FREE(cigar1);
	  FREE(md_string_1);
	} else {
	  debug3(printf("** "));
	  FREE(best_cigar);
	  FREE(best_md_string);

	  /* ? Is this right for best_chrpos */
	  best_chrpos = wildtype_chrpos[i][readalign_start - chrpos_start + finalc1];
	  best_cigar = cigar1;
	  best_md_string = md_string_1;
	  best_score = score1;
	}
      }
      debug3(printf("\n"));

      if (best_score > score0) {
	for (linei = linei_start; linei < linei_end; linei++) {
	  bamline = bamlines[linei];

	  orig_cigar = Bamline_cigar_string(bamline);
	  if (strcmp(best_cigar,orig_cigar)) {
	    /* Aligns better to genome with indel and changes CIGAR string */
	    orig_nmismatches = read_nmismatches(Bamline_chrpos_low(bamline),Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
						Bamline_read(bamline),genome,chroffset);

	    new_cigar_types = Samread_parse_cigar(&new_cigar_npositions,&new_readlength,best_cigar);
	    new_nmismatches = read_nmismatches(best_chrpos,new_cigar_types,new_cigar_npositions,
					       Bamline_read(bamline),genome,chroffset);
	    debug2(printf(" => orig nmismatches:%d vs new_nmismatches:%d\n",orig_nmismatches,new_nmismatches));

	    if (new_nmismatches <= orig_nmismatches) {
	      Bamline_print_new_cigar(stdout,bamline,best_chrpos,best_cigar,best_md_string,
				      /*quality_score_adj*/33);
	    }
	    Uintlist_free(&new_cigar_npositions);
	    Intlist_free(&new_cigar_types);
	  }
	  FREE(orig_cigar);
	}
      }

      FREE(best_md_string);
      FREE(best_cigar);
    }

    for (linei = 0; linei < nlines; linei++) {
      bamline = bamlines[linei];
      Bamline_free(&bamline);
    }
    FREE(bamlines);
  }

  Bamread_unlimit_region(realign_bamreader);
  return;
}


/* Assume we are down to our best two candidates */
static void
realign_reads_two (Bamreader_T realign_bamreader, char *chr, Genomicpos_T indelpos_start, Genomicpos_T indelpos_end,
		   Genomicpos_T chrpos_start, Genomicpos_T chrpos_end, Genome_T genome, Genomicpos_T chroffset,
		   char *genomicseg0, char *genomicseg_alt, int glength_alt,
		   int *nindels0, int *nindels_alt, char **delstrings_alt, Genomicpos_T *wildtype_chrpos,
		   Genomicpos_T var_first_chrpos, Genomicpos_T var_last_chrpos, int max_varlength, Dynprog_T dynprog) {
  Bamline_T bamline, bamline_rep, *bamlines, prev_bamline = NULL;
  int nlines, linei_start, linei_end, linei;
  char *rsequence;
  int readlength, new_readlength;
  int orig_nmismatches, new_nmismatches;
  Intlist_T new_cigar_types;
  Uintlist_T new_cigar_npositions;

  int glength0;
  Genomicpos_T chrpos1;
  int finalc0, finalc1;
  int score0, score1;
  char *orig_cigar, *cigar0, *cigar1;
  char *md_string_0, *md_string_1;

  Genomicpos_T readalign_start, readalign_end, total_ins;
  int readalign_length;

	       
  glength0 = chrpos_end - chrpos_start + 1;
  Bamread_limit_region(realign_bamreader,chr,indelpos_start,indelpos_end);

  while ((bamlines = Bamread_next_bamline_set(&nlines,&prev_bamline,realign_bamreader,
					      /*desired_read_group*/NULL,/*minimum_mapq*/0,/*good_unique_mapq*/0,
					      /*maximum_nhits*/1000000,/*need_unique_p*/false,/*need_primary_p*/false,
					      /*ignore_duplicates_p*/false,/*need_concordant_p*/false)) != NULL) {

    debug2(printf("Got a set of %d lines\n",nlines));
    for (linei_start = 0; linei_start < nlines; linei_start = linei_end) {
      linei_end = linei_start + 1;
      while (linei_end < nlines && !strcmp(Bamline_read(bamlines[linei_start]),Bamline_read(bamlines[linei_end]))) {
	linei_end++;
      }
      debug2(printf("Analyzing %d to %d\n",linei_start,linei_end - 1));

      bamline_rep = bamlines[linei_start];

      debug2(printf("\n"));
      debug2(printf(">%s:%u..%u ",chr,Bamline_chrpos_low(bamline_rep),Bamline_chrpos_high(bamline_rep)));
      debug2(Bamread_print_cigar(stdout,bamline_rep));

      rsequence = Bamline_read(bamline_rep);
      readlength = Bamline_readlength(bamline_rep);

#ifdef WHOLE_ALIGNMENT
      readalign_start = chrpos_start;
      readalign_length = glength0;
#else
      total_ins = Bamline_total_ins(bamline_rep);
      if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	readalign_start = chrpos_start;
      }
      if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_end) {
	readalign_end = chrpos_end;
      }
#endif

      if (var_last_chrpos < readalign_start) {
	/* Skip, because variants are all to the left of the alignment */
      } else if (var_first_chrpos > readalign_end) {
	/* Skip, because variants are all to the right of the alignment */
      } else {
	/* Align */
	readalign_length = readalign_end - readalign_start + 1;

#ifdef NOGAP
	cigar0 = Dynprog_single_nogap(&score0,&finalc0,&md_string_0,dynprog,rsequence,
				      &(genomicseg0[readalign_start - chrpos_start]),
				      &(genomicseg0[readalign_start - chrpos_start]),
				      &(nindels0[readalign_start - chrpos_start]),/*deletion_string*/NULL,
				      /*rlength*/readlength,/*glength*/readalign_length,
				      /*extraband_single*/total_ins + max_varlength);
#else
	/* Use gap here to account for sequence errors */
	cigar0 = Dynprog_single_gap(&score0,&finalc0,&md_string_0,dynprog,rsequence,
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(genomicseg0[readalign_start - chrpos_start]),
				    &(nindels0[readalign_start - chrpos_start]),/*deletion_string*/NULL,
				    /*rlength*/readlength,/*glength*/readalign_length,
				    /*extraband_single*/total_ins + max_varlength);
#endif
	/* chrpos0 = chrpos_start + finalc0; */
	debug2(printf(" wildtype:%d at %d (%u) %s",score0,finalc0,readalign_start+finalc0,cigar0));


#ifdef WHOLE_ALIGNMENT
	readalign_start = chrpos_start;
	readalign_length = glength_alt;
#else
	/* if ((readalign_start = Bamline_chrpos_low_noclip(bamline_rep) - total_ins - max_varlength) < chrpos_start) {
	   readalign_start = chrpos_start;
	   } */
	if ((readalign_end = Bamline_chrpos_high_noclip(bamline_rep) + total_ins + max_varlength) > chrpos_start + glength_alt) {
	  readalign_end = chrpos_start + glength_alt;
	}
	readalign_length = readalign_end - readalign_start + 1;
#endif

#ifdef NOGAP
	cigar1 = Dynprog_single_nogap(&score1,&finalc1,&md_string_1,dynprog,rsequence,
				      &(genomicseg_alt[readalign_start - chrpos_start]),
				      &(genomicseg_alt[readalign_start - chrpos_start]),
				      &(nindels_alt[readalign_start - chrpos_start]),
				      &(delstrings_alt[readalign_start - chrpos_start]),
				      /*rlength*/readlength,/*glength*/readalign_length,
				      /*extraband_single*/total_ins + max_varlength);
#else
	/* Use gap here to account for sequence errors */
	cigar1 = Dynprog_single_gap(&score1,&finalc1,&md_string_1,dynprog,rsequence,
				    &(genomicseg_alt[readalign_start - chrpos_start]),
				    &(genomicseg_alt[readalign_start - chrpos_start]),
				    &(nindels_alt[readalign_start - chrpos_start]),
				    &(delstrings_alt[readalign_start - chrpos_start]),
				    /*rlength*/readlength,/*glength*/readalign_length,
				    /*extraband_single*/total_ins + max_varlength);
#endif
	/* ? Is this right for chrpos1 */
	chrpos1 = wildtype_chrpos[readalign_start - chrpos_start + finalc1];
	debug2(printf("  alt:%d at %d (%u) %s",score1,finalc1,readalign_start+finalc1,cigar1));
	debug2(printf("\n"));

	if (score1 > score0) {
	  for (linei = linei_start; linei < linei_end; linei++) {
	    bamline = bamlines[linei];

	    orig_cigar = Bamline_cigar_string(bamline);
	    if (strcmp(cigar1,orig_cigar)) {
	      /* Aligns better to genome with indel and changes CIGAR string */
	      orig_nmismatches = read_nmismatches(Bamline_chrpos_low(bamline),Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
						  Bamline_read(bamline),genome,chroffset);

	      new_cigar_types = Samread_parse_cigar(&new_cigar_npositions,&new_readlength,cigar1);
	      new_nmismatches = read_nmismatches(chrpos1,new_cigar_types,new_cigar_npositions,
						 Bamline_read(bamline),genome,chroffset);
	      debug2(printf(" => orig nmismatches:%d vs new_nmismatches:%d\n",orig_nmismatches,new_nmismatches));

	      if (new_nmismatches <= orig_nmismatches) {
		Bamline_print_new_cigar(stdout,bamline,chrpos1,cigar1,md_string_1,/*quality_score_adj*/33);
	      }
	      Uintlist_free(&new_cigar_npositions);
	      Intlist_free(&new_cigar_types);
	    }

	    FREE(orig_cigar);
	  }
	}

	FREE(md_string_1);
	FREE(md_string_0);
	FREE(cigar1);
	FREE(cigar0);
      }
    }

    for (linei = 0; linei < nlines; linei++) {
      bamline = bamlines[linei];
      Bamline_free(&bamline);
    }
    FREE(bamlines);
  }

  Bamread_unlimit_region(realign_bamreader);
  return;
}


static List_T
get_segments (Bamline_T bamline, Genome_T genome, Genomicpos_T chroffset) {
  List_T segments = NULL;
  Genomicpos_T pos;
  char *segment, *p;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;


  pos = Bamline_chrpos_low(bamline) - 1U; /* Bamread reports chrpos as 1-based */

  p = Bamline_read(bamline);

  for (a = Bamline_cigar_types(bamline), b = Bamline_cigar_npositions(bamline);
       a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S' /* && max_softclip == 0 */) {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;

    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */

    } else if (type == 'N') {
      fprintf(stderr,"bam_indelfix can't handle splicing yet\n");
      exit(9);

      /* revise_diffcigar(pos,var_types,var_chrpos,this); */
      pos += mlength;

    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'I') {
      segment = (char *) CALLOC(mlength+1,sizeof(char));
      strncpy(segment,p,mlength);
      segments = List_push(segments,(void *) segment);
      debug9(printf("Saw insertion of %d (%s) at %u\n",mlength,segment,chroffset+pos+1));

      p += mlength;

    } else if (type == 'D') {
      debug(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
      segment = (char *) CALLOC(mlength+1,sizeof(char));
      Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,segment);
      segments = List_push(segments,(void *) segment);
      debug9(printf("Saw deletion of %d (%s) at %u\n",mlength,segment,chroffset+pos+1));

      pos += mlength;

    } else if (type == 'M' /* || (type == 'S' && max_softclip > 0) */) {
      p += mlength;
      pos += mlength;

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }

  return List_reverse(segments);
}


static void
process_neighborhood (int *nsegments, List_T var_neighborhood, Variant_T key_var,
		      Dynprog_T dynprog, Genome_T genome, Genomicpos_T chroffset,
		      Bamreader_T realign_bamreader, char *printchr, bool allow_multiple_p) {
  Intlist_T p;
  Uintlist_T q, r;
  List_T s;
  int nvariants, keyi = -1, i;
  Genomicpos_T chrpos, indelpos_start = -1U, indelpos_end = 0U, chrpos_start, chrpos_end;
  Genomicpos_T *initpos_alt;
  Variant_T var, *var_array;
  List_T ptr;
  int max_readlength, max_varlength, inslength, dellength;
  int c;
  char *genomicseg0, **genomicsegs_alt, *dest, *indel;
  int *nindels0, **nindels_alt, varlength;
  int glength0, glength1, *glengths_alt;
  char ***delstrings_alt;
  Genomicpos_T **wildtype_chrpos, *cptr;
  int besti;


  nvariants = List_length(var_neighborhood);
  for (ptr = var_neighborhood; ptr != NULL; ptr = List_next(ptr)) {
    var = (Variant_T) List_head(ptr);
    if ((chrpos = Uintlist_head(var->chrpos) - Uintlist_head(var->npositions)) < indelpos_start) {
      indelpos_start = chrpos;
    }
    if ((chrpos = Uintlist_last_value(var->chrpos) + Uintlist_last_value(var->npositions)) > indelpos_end) {
      indelpos_end = chrpos;
    }
  }

  max_readlength = find_max_readlength(realign_bamreader,printchr,indelpos_start,indelpos_end);

  chrpos_start = indelpos_start - max_readlength;
  chrpos_end = indelpos_end + max_readlength;
  debug(printf("chrpos_start %u, chrpos_end %u\n",chrpos_start,chrpos_end));
  debug9(printf("In block %u..%u, testing segment at indelpos %u..%u => chrpos %u..%u\n",
		blockstart,blockstart+lasti,indelpos_start,indelpos_end,chrpos_start,chrpos_end));
  *nsegments += 1;

  /* wildtype */
  glength0 = chrpos_end - chrpos_start + 1;
  nindels0 = (int *) CALLOC(glength0,sizeof(int));
  dest = genomicseg0 = (char *) MALLOC((glength0 + 1)*sizeof(char));
  for (chrpos = chrpos_start; chrpos <= chrpos_end; chrpos++) {
    *dest++ = Genome_get_char(genome,chroffset+chrpos-1U);
  }
  *dest = '\0';
  debug(printf("wildtype:     %s\n",genomicseg0));

  var_array = (Variant_T *) List_to_array(var_neighborhood,NULL);
  for (i = 0; i < nvariants; i++) {
    if (var_array[i] == key_var) {
      keyi = i;
    }
  }
  /* qsort(vars_array,nvariants,sizeof(Variant_T),Variant_count_cmp); */

  genomicsegs_alt = (char **) CALLOC(nvariants,sizeof(char *));
  initpos_alt = (Genomicpos_T *) CALLOC(nvariants,sizeof(Genomicpos_T));
  nindels_alt = (int **) CALLOC(nvariants,sizeof(int *));
  glengths_alt = (int *) CALLOC(nvariants,sizeof(int));
  delstrings_alt = (char ***) CALLOC(nvariants,sizeof(char **));
  wildtype_chrpos = (Genomicpos_T **) CALLOC(nvariants,sizeof(Genomicpos_T *));

  max_varlength = 0;
  for (i = 0; i < nvariants; i++) {
    inslength = 0;
    dellength = 0;
    var = var_array[i];

    varlength = 0;
    for (p = var->types, r = var->npositions; p != NULL;
	 p = Intlist_next(p), r = Uintlist_next(r)) {
      if (Intlist_head(p) == 'I') {
	varlength += Uintlist_head(r);
	inslength += Uintlist_head(r);
      } else if (Intlist_head(p) == 'D') {
	varlength -= Uintlist_head(r);
	dellength += Uintlist_head(r);
      } else {
	fprintf(stderr,"Unexpected var type %c\n",Intlist_head(p));
	abort();
      }
      if (inslength > max_varlength) {
	max_varlength = inslength;
      }
      if (dellength > max_varlength) {
	max_varlength = dellength;
      }
    }

    glength1 = glengths_alt[i] = glength0 + varlength;
    initpos_alt[i] = Uintlist_head(var->chrpos);
    nindels_alt[i] = (int *) CALLOC(glength1,sizeof(int));
    delstrings_alt[i] = (char **) CALLOC(glength1,sizeof(char *));
    dest = genomicsegs_alt[i] = (char *) MALLOC((glength1 + 1)*sizeof(char));
    cptr = wildtype_chrpos[i] = (Genomicpos_T *) MALLOC(glength1*sizeof(Genomicpos_T));

    chrpos = chrpos_start;
    c = 0;
    for (p = var->types, q = var->chrpos, r = var->npositions, s = var->segments; p != NULL;
	 p = Intlist_next(p), q = Uintlist_next(q), r = Uintlist_next(r), s = List_next(s)) {

      while (chrpos < Uintlist_head(q)) {
	*dest++ = Genome_get_char(genome,chroffset+chrpos-1U);
	cptr[c++] = chrpos++;
      }

      if (Intlist_head(p) == 'I') {
	indel = (char *) List_head(s);
	assert(strlen(indel) == Uintlist_head(r));
	debug(printf("Insertion segment is %s\n",indel));

	while (*indel != '\0') {
	  *dest++ = *indel++;
	  /* Put mlength instead of 1, so dynprog can use soft-clipping for incomplete insertions */
	  cptr[c] = chrpos;	/* Do not increment chrpos */
	  nindels_alt[i][c++] = +(Uintlist_head(r));
	  debug(printf("Putting %d at nindels_alt[%i][%d]\n",nindels_alt[i][c-1],i,c-1));
	}

      } else if (Intlist_head(p) == 'D') {
	delstrings_alt[i][c-1] = indel = (char *) List_head(s);
	assert(strlen(indel) == Uintlist_head(r));

	nindels_alt[i][c-1] = -(Uintlist_head(r));
	chrpos += strlen(indel);
	debug(printf("Putting %d at nindels_alt[%i][%d]\n",nindels_alt[i][c-1],i,c-1));

      } else {
	abort();
      }
    }

    while (chrpos <= chrpos_end) {
      *dest++ = Genome_get_char(genome,chroffset+chrpos-1U);
      cptr[c++] = chrpos++;
    }

    *dest = '\0';
    assert(c == glength1);
    assert(strlen(genomicsegs_alt[i]) == glength1);
    debug(printf("variant-%d: %s\n",i,genomicsegs_alt[i]));
  }

  if (allow_multiple_p == true) {
    /* Align each read against the alternatives */
    realign_reads_multiple(realign_bamreader,printchr,indelpos_start,indelpos_end,chrpos_start,chrpos_end,
			   genome,chroffset,genomicseg0,genomicsegs_alt,glengths_alt,
			   nindels0,nindels_alt,delstrings_alt,wildtype_chrpos,nvariants,
			   max_varlength,dynprog);

  } else {
    if (nvariants == 1) {
      besti = 0;
    } else {
      besti = tally_reads_multiple(realign_bamreader,printchr,indelpos_start,indelpos_end,chrpos_start,chrpos_end,
				   genomicseg0,genomicsegs_alt,glengths_alt,initpos_alt,
				   nindels0,nindels_alt,delstrings_alt,var_array,nvariants,keyi,
				   max_varlength,dynprog);


    }
    if (besti >= 0) {
      realign_reads_two(realign_bamreader,printchr,indelpos_start,indelpos_end,chrpos_start,chrpos_end,
			genome,chroffset,genomicseg0,genomicsegs_alt[besti],glengths_alt[besti],
			nindels0,nindels_alt[besti],delstrings_alt[besti],wildtype_chrpos[besti],
			/*var_first_chrpos*/var_array[besti]->first_chrpos,
			/*var_last_chrpos*/var_array[besti]->last_chrpos,max_varlength,dynprog);
    }
  }

  FREE(var_array);

  FREE(nindels0);
  FREE(genomicseg0);

  for (i = 0; i < nvariants; i++) {
    FREE(delstrings_alt[i]);
    FREE(nindels_alt[i]);
    FREE(genomicsegs_alt[i]);
    FREE(wildtype_chrpos[i]);
  }
  FREE(delstrings_alt);
  FREE(nindels_alt);
  FREE(initpos_alt);
  FREE(genomicsegs_alt);
  FREE(wildtype_chrpos);
  FREE(glengths_alt);

  return;
}




/* Modeled after Bamtally_run */
void
Indelfix_run (Bamreader_T bamreader, int max_rlength, int max_glength, Genome_T genome, char *printchr,
	      Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      bool ignore_lowend_p, bool ignore_highend_p, Genomicpos_T neighborhood_size,
	      bool allow_multiple_p, char *bamfile) {
  int nsegments = 0;
  Bamline_T bamline;
  Bamreader_T realign_bamreader = NULL;
  Dynprog_T dynprog;

  List_T *variants_by_chr, all_variants, var_neighborhood;
  Variant_T *var_array, var, key_var;
  int nvariants, vari;
  Intlist_T var_types;
  Uintlist_T var_npositions, var_chrpos;
  Genomicpos_T first_chrpos, startpos, endpos, chrpos;
  List_T segments, p;
  char *segment;

  Genomicpos_T last_notice = 0;


  debug0(printf("Indelfix_run on chromosome %s:%u..%u\n",printchr,chrstart,chrend));

  /* 1.  Scan for variants */
  variants_by_chr = (List_T *) CALLOC(chrend+1,sizeof(List_T));

  while ((bamline = Bamread_next_indel_bamline(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
					       need_unique_p,need_primary_p,ignore_duplicates_p,
					       need_concordant_p)) != NULL) {
    if (Bamline_chrpos_low(bamline) > last_notice + 10000000) {
      fprintf(stderr,">%s:%u..%u\n",printchr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline));
      last_notice = Bamline_chrpos_low(bamline);
    }
    debug0(printf(">%s:%u..%u ",printchr,Bamline_chrpos_low(bamline),Bamline_chrpos_high(bamline)));
    debug0(Bamread_print_cigar(stdout,bamline));
    debug0(printf("\n"));

    if (ignore_lowend_p == true && Bamline_lowend_p(bamline) == true) {
      /* Skip */
    } else if (ignore_highend_p == true && Bamline_lowend_p(bamline) == false) {
      /* Skip */
    } else {
      var_types = Bamline_diffcigar(&var_npositions,&var_chrpos,bamline);
      first_chrpos = Uintlist_head(var_chrpos);
      segments = get_segments(bamline,genome,chroffset);
      if (variants_by_chr[first_chrpos] == NULL) {
	variants_by_chr[first_chrpos] = List_push(NULL,(void *) Variant_new(var_types,var_npositions,var_chrpos,segments));
      } else if ((var = Variant_find(variants_by_chr[first_chrpos],var_types,var_npositions,var_chrpos,segments)) == NULL) {
	variants_by_chr[first_chrpos] = List_push(variants_by_chr[first_chrpos],
						  (void *) Variant_new(var_types,var_npositions,var_chrpos,segments));
      } else {
	debug9(printf("Saw this variant already\n"));
	var->count += 1;
	Intlist_free(&var_types);
	Uintlist_free(&var_npositions);
	Uintlist_free(&var_chrpos);
	for (p = segments; p != NULL; p = List_next(p)) {
	  segment = (char *) List_head(p);
	  FREE(segment);
	}
	List_free(&segments);
      }
    }

    Bamline_free(&bamline);
  }

  /* 2.  Obtain array of variants, sorted by count */
  all_variants = (List_T) NULL;
  for (chrpos = 0; chrpos <= chrend; chrpos++) {
    if (variants_by_chr[chrpos] != NULL) {
      for (p = variants_by_chr[chrpos]; p != NULL; p = List_next(p)) {
	all_variants = List_push(all_variants,(void *) List_head(p));
      }
    }
  }

  var_array = (Variant_T *) List_to_array_n(&nvariants,all_variants);
  qsort(var_array,nvariants,sizeof(Variant_T),Variant_count_cmp);
#if 0
  for (vari = 0; vari < nvariants; vari++) {
    Variant_print(var_array[vari]);
  }
#endif
  List_free(&all_variants);
  fprintf(stderr,"Found %d variants in chromosome %s\n",nvariants,printchr);


  /* 3.  Process variants by neighborhood */
  realign_bamreader = Bamread_new(bamfile);
  dynprog = Dynprog_new(max_rlength,max_glength,/*doublep*/false);

  for (vari = 0; vari < nvariants; vari++) {
    key_var = var_array[vari];
    if (key_var->count == 0) {
      /* Processed already */
    } else {
      if (key_var->first_chrpos < neighborhood_size) {
	startpos = 0;
      } else {
	startpos = key_var->first_chrpos - neighborhood_size;
      }
      if (key_var->last_chrpos + neighborhood_size > chrend) {
	endpos = chrend;
      } else {
	endpos = key_var->last_chrpos + neighborhood_size;
      }

      debug8(printf("Neighborhood %u..%u centered by %u..%u with count %d\n",
		    startpos,endpos,key_var->first_chrpos,key_var->last_chrpos,key_var->count));
      var_neighborhood = (List_T) NULL;
      for (chrpos = startpos; chrpos <= endpos; chrpos++) {
	if (variants_by_chr[chrpos] != NULL) {
	  for (p = variants_by_chr[chrpos]; p != NULL; p = List_next(p)) {
	    var = (Variant_T) List_head(p);
	    if (var->count > 0) {
	      debug8(Variant_print(var));
	      var_neighborhood = List_push(var_neighborhood,(void *) var);
	    }
	  }
	}
      }
      debug8(printf("\n"));
      var_neighborhood = List_reverse(var_neighborhood);
      fprintf(stderr,"Neighborhood %u..%u centered by %u..%u with count %d: %d variants\n",
	      startpos,endpos,key_var->first_chrpos,key_var->last_chrpos,key_var->count,List_length(var_neighborhood));

      process_neighborhood(&nsegments,var_neighborhood,key_var,dynprog,
			   genome,chroffset,realign_bamreader,printchr,allow_multiple_p);
      debug0(fprintf(stderr,"Processed segment %d\n",nsegments));

      List_free(&var_neighborhood);
    }
  }

  FREE(var_array);

  for (chrpos = 0; chrpos <= chrend; chrpos++) {
    if (variants_by_chr[chrpos] != NULL) {
      for (p = variants_by_chr[chrpos]; p != NULL; p = List_next(p)) {
	var = (Variant_T) List_head(p);
	Variant_free(&var);
      }
      List_free(&(variants_by_chr[chrpos]));
    }
  }
  FREE(variants_by_chr);


  Dynprog_free(&dynprog);
  Bamread_free(&realign_bamreader);

  return;
}

