static char rcsid[] = "$Id: stage2.c 135441 2014-05-07 22:15:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "intlist.h"
#include "diag.h"
#include "oligoindex_hr.h"
#include "genome_hr.h"
#include "genome_sites.h"
#include "complement.h"
#include "maxent_hr.h"


/* Tests whether genomicseg == query in convert_to_nucleotides, and
   whether oligoindex_hr gives same results as oligoindex */
/* #define EXTRACT_GENOMICSEG 1 */

#define USE_DIAGPOOL 1

/* #define SQUARE 1 */

#define SUFF_PCTCOVERAGE_OLIGOINDEX 0.90

/* #define SUFF_PCTCOVERAGE_STAGE2 0.10 */
#define SUFF_NCOVERED 200
#define SUFF_MAXNCONSECUTIVE 20

#define MAX_NACTIVE 100	/* 100 previously considered too low, but may
			   be okay in conjunction with
			   diagonalization */
#define MAX_GRAND_LOOKBACK 200

/* Penalty for genomic distances */

#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16

#define NINTRON_PENALTY_MISMATCH 8
#define NON_CANONICAL_PENALTY_ENDS 4
#define NON_CANONICAL_PENALTY_MIDDLE 4
#define MISS_BEHIND 16
#define GREEDY_ADVANCE 6
#define FINAL_SCORE_TOLERANCE 8

#define ENOUGH_CONSECUTIVE 32

#define INFINITE 1000000

/* EQUAL_DISTANCE used to be 3 for PMAP and 6 for GMAP, but that
   allowed indels in repetitive regions.  Now have separate
   variables. */
#ifdef PMAP
#define EQUAL_DISTANCE_FOR_CONSECUTIVE 0
#define EQUAL_DISTANCE_NOT_SPLICING 3
#else
#define EQUAL_DISTANCE_FOR_CONSECUTIVE 0
#define EQUAL_DISTANCE_NOT_SPLICING 9
#endif


#define INTRON_DEFN 9		/* Cannot exceed 9 */
#define NEAR_END_LENGTH 20	/* Determines whether to ignore EXON_DEFN at ends */
#define EXON_DEFN 30
#define MAX_SKIPPED 3

#define SCORE_FOR_RESTRICT 10
/* #define SUFFICIENT_ROOTNLINKS 10  */ /* Setting this too low can slow down program considerably */


#ifdef PMAP
#define SAMPLE_INTERVAL 1
#define NT_PER_MATCH 3
#define CONSEC_POINTS_PER_MATCH 3 /* Possible increase to reward consecutiveness */
#define NONCODON_INDEL_PENALTY 15
#else
#define SAMPLE_INTERVAL 2	/* For cases where adjacentp == false.
				   Means that we can find islands of
				   9-mers */
#define NT_PER_MATCH 1
#define NT_PER_CODON 3
#define CONSEC_POINTS_PER_MATCH 1 /* Possible increase to reward consecutiveness */
#define CONSEC_POINTS_PER_CODON 3 /* Possible increase to reward consecutiveness */
#endif

#define SHIFT_EXTRA 15

#define ONE 1.0
#define TEN_THOUSAND 10000.0
#define HUNDRED_THOUSAND 100000.0
#define ONE_MILLION 1000000.0



static bool splicingp;
static bool use_canonical_middle_p;
static bool use_canonical_ends_p;
static int suboptimal_score_end;
static int suboptimal_score_start;
static int min_intronlength;
static Mode_T mode;
static bool snps_p;

void
Stage2_setup (bool splicingp_in, bool cross_species_p,
	      int suboptimal_score_start_in, int suboptimal_score_end_in,
	      int min_intronlength_in, Mode_T mode_in, bool snps_p_in) {
  splicingp = splicingp_in;
  if (splicingp == true) {
    use_canonical_ends_p = true;
  } else {
    use_canonical_ends_p = false;
  }
  if (cross_species_p == true) {
    use_canonical_middle_p = true;
  } else {
    use_canonical_middle_p = false;
  }
  suboptimal_score_start = suboptimal_score_start_in;
  suboptimal_score_end = suboptimal_score_end_in;
  min_intronlength = min_intronlength_in;
  mode = mode_in;
  snps_p = snps_p_in;
  return;
}


/* General */
#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif

/* Final results of stage 2 */
#ifdef DEBUG0
#define debug0(x) x
#else 
#define debug0(x)
#endif

/* Print all links */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* For generating a graph */
#ifdef DEBUG3
#define debug3(x) x
#else 
#define debug3(x)
#endif

/* Converting to nucleotides */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* revise_active */
#ifdef DEBUG6
#define debug6(x) x
#else 
#define debug6(x)
#endif

/* Shifted canonical */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
#endif

/* find_canonical_dinucleotides */
#ifdef DEBUG8
#define debug8(x) x
#else 
#define debug8(x)
#endif

/* Dynamic programming */
/* Can also define debug9(x) as: if (querypos == XX) {x;} */
#ifdef DEBUG9
#define debug9(x) x
#else 
#define debug9(x)
#endif

/* Grand winner */
#ifdef DEBUG10
#define debug10(x) x
#else 
#define debug10(x)
#endif

/* Multiple alignments */
#ifdef DEBUG11
#define debug11(x) x
#else 
#define debug11(x)
#endif

/* Grand winner */
#ifdef DEBUG12
#define debug12(x) x
#else 
#define debug12(x)
#endif


#define T Stage2_T
struct T {
  List_T middle;
  List_T all_starts;
  List_T all_ends;
};


void
Stage2_free (T *old) {
  /* List_free(&(*old)->middle); -- Not necessary because of pairpool */
  List_free(&(*old)->all_starts);
  List_free(&(*old)->all_ends);
  FREE(*old);
  return;
}

static T
Stage2_new (List_T middle, List_T all_starts, List_T all_ends) {
  T new = (T) MALLOC(sizeof(*new));

  new->middle = middle;
  new->all_starts = all_starts;
  new->all_ends = all_ends;

  return new;
}

List_T
Stage2_middle (T this) {
  return this->middle;
}

List_T
Stage2_all_starts (T this) {
  return this->all_starts;
}

List_T
Stage2_all_ends (T this) {
  return this->all_ends;
}


/************************************************************************/

typedef struct Link_T *Link_T;
struct Link_T {
  int fwd_consecutive;
  int fwd_rootposition;
  /*int fwd_rootnlinks;*/		/* Number of links in last branch */
  int fwd_score;

  int fwd_pos;
  int fwd_hit;

#ifdef DEBUG9
  int fwd_tracei;		/* Corresponds to a distinct set of branches */
  int fwd_intronnfwd;
  int fwd_intronnrev;
  int fwd_intronnunk;
#endif

#ifdef SEPARATE_FWD_REV
  /* No longer checking separate fwd/rev directions */
  int rev_consecutive;
  int rev_rootposition;
  /*int rev_rootnlinks;*/		/* Number of links in last branch */
  int rev_score;

  int rev_pos;
  int rev_hit;

#ifdef DEBUG9
  int rev_tracei;		/* Corresponds to a distinct set of branches */
  int rev_intronnfwd;
  int rev_intronnrev;
  int rev_intronnunk;
#endif	/* rev */

#endif
};


/* lengths2 is has length1 entries.  Note that lengths2 may have
   negative entries */
static struct Link_T **
Linkmatrix_1d_new (int length1, int *lengths2, int totallength) {
  struct Link_T **links;
  int i;
  
  /* Outer dimension can be MALLOC, but inner one must be CALLOC */
  links = (struct Link_T **) MALLOC(length1 * sizeof(struct Link_T *));
  links[0] = (struct Link_T *) CALLOC(totallength,sizeof(struct Link_T));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] < 0) {
      links[i] = links[i-1];
    } else {
      links[i] = &(links[i-1][lengths2[i-1]]);
    }
  }
  return links;
}

static void
Linkmatrix_1d_free (struct Link_T ***links) {
  FREE((*links)[0]);
  FREE(*links);
  return;
}


static struct Link_T **
Linkmatrix_2d_new (int length1, int *lengths2) {
  struct Link_T **links;
  int i;
  
  links = (struct Link_T **) CALLOC(length1,sizeof(struct Link_T *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      links[i] = (struct Link_T *) NULL;
    } else {
      links[i] = (struct Link_T *) CALLOC(lengths2[i],sizeof(struct Link_T));
    }
  }
  return links;
}

static void
Linkmatrix_2d_free (struct Link_T ***links, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*links)[i]) {
      FREE((*links)[i]);
    }
  }
  FREE(*links);
  return;
}



#ifdef DEBUG1
#ifdef SEPARATE_FWD_REV
static void
Linkmatrix_print_both (struct Link_T **links, Chrpos_T **mappings, int length1,
		       int *npositions, char *queryseq_ptr, int indexsize) {
  int i, j;
  char *oligo;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  for (i = 0; i <= length1-indexsize; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %d.%u:%d(%d,%d)[%u]-%d(%d,%d)[%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_tracei,
	     links[i][j].rev_score,
	     links[i][j].rev_pos,links[i][j].rev_hit,links[i][j].rev_tracei);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

#else

/* For PMAP, indexsize is in aa */
static void
Linkmatrix_print_fwd (struct Link_T **links, Chrpos_T **mappings, int *npositions,
		      char *queryseq_ptr, int indexsize) {
  int i, j, lastpos;
  char *oligo;
  Intlist_T p, q;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  lastpos = length1 - indexsize;

  for (i = 0; i <= lastpos; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %d.%u:%d(%d,%d)[%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_tracei);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

#endif
#endif

static void
mappings_dump_R (Chrpos_T **mappings, int *npositions, int length1,
		 int **active, int *firstactive, int indexsize, char *varname) {
  int querypos;
  int lastpos, hit;
  bool printp = false;

  lastpos = length1 - indexsize;
  printf("%s <- matrix(c(\n",varname);
  for (querypos = 0; querypos < lastpos; querypos++) {
    if (firstactive) {
      if (mappings[querypos] != NULL) {
	hit = firstactive[querypos];
	while (hit != -1) {
	  /* Last elt is for score */
	  if (printp == false) {
	    printp = true;
	  } else {
	    printf(",\n");
	  }
	  printf("%d,%d,%d,%d",querypos,mappings[querypos][hit],
		 hit,active[querypos][hit]);
	  hit = active[querypos][hit];
	}
      }
    } else {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (printp == false) {
	  printp = true;
	} else {
	  printf(",\n");
	}
	printf("%d,%d,%d",querypos,mappings[querypos][hit],hit);
      }
    }
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}
    

static void
best_path_dump_R (struct Link_T **links, Chrpos_T **mappings,
		  int querypos, int hit, bool fwdp, char *varname) {
  Chrpos_T position;
  int prev_querypos, prevhit, save_querypos, savehit;
  bool printp = false;

  save_querypos = querypos;
  savehit = hit;

  printf("%s <- matrix(c(\n",varname);
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d,%d",querypos,position);

    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifdef SEPARATE_FWD_REV
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
  }
  printf("),ncol=2,byrow=T)\n");

  querypos = save_querypos;
  hit = savehit;

  printp = false;
  printf("%s <- matrix(c(\n","scores");
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    if (fwdp == true) {
      printf("%d,%d",querypos,links[querypos][hit].fwd_score);
#ifdef SEPARATE_FWD_REV
    } else {
      printf("%d,%d",querypos,links[querypos][hit].rev_score);
#endif
    }

    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifdef SEPARATE_FWD_REV
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}

static void
active_bounds_dump_R (Chrpos_T *minactive, Chrpos_T *maxactive,
		      int querylength) {
  int querypos;
  bool printp = false;

  printf("querypos <- 0:%d\n",querylength-1);
  printf("%s <- c(\n","minactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",minactive[querypos]);
  }
  printf(")\n");

  printp = false;
  printf("%s <- c(\n","maxactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",maxactive[querypos]);
  }
  printf(")\n");

  return;
}


#if 0
#define diffdist_penalty_nosplicing(x) (x + 1)
#define diffdist_penalty_splicing(x) (x/TEN_THOUSAND + 1)
#endif


#ifdef PMAP
#define QUERYDIST_PENALTY_FACTOR 2
#else
#define QUERYDIST_PENALTY_FACTOR 8
#endif

#if 0
/* querydistance already has indexsize_nt subtracted */
static int
querydist_penalty (int querydistance) {
#ifdef PMAP
  return querydistance/2;
#else
  return querydistance/8;
#endif
}
#endif




/************************************************************************
 *  Procedures for finding canonical introns quickly
 ************************************************************************/

#ifdef DEBUG8

static void
print_last_dinucl (int *last_dinucl, int genomiclength) {
  int pos;

  for (pos = 0; pos < genomiclength - 3 + SHIFT_EXTRA; pos++) {
    printf("%d %d\n",pos,last_dinucl[pos]);
  }
  printf("\n");

  return;
}

#endif



#if 0
static void
check_canonical_dinucleotides_hr (int *lastGT, int *lastAG, 
#ifndef PMAP
				  int *lastCT, int *lastAC,
#endif
				  Univcoord_T genomicstart, Univcoord_T genomicend, bool plusp,
				  int genomiclength) {

  int pos;
  int lastpos1, lastpos2;

  /* printf("genomicstart %u, genomicend %u, genomiclength %d\n",genomicstart,genomicend,genomiclength); */

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastGT[pos];
    lastpos2 = Genome_prev_donor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("donor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastAG[pos];
    lastpos2 = Genome_prev_acceptor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("acceptor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#ifndef PMAP
#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastAC[pos];
    lastpos2 = Genome_prev_antidonor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("antidonor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastCT[pos];
    lastpos2 = Genome_prev_antiacceptor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("antiacceptor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif
#endif

  return;
}
#endif


/* assert(chrstart < chrend) */
/* For plus, chrinit = chrstart, chrterm = chrend.  For minus, chrinit = (chrhigh - chroffset) - chrend, chrterm = (chrhigh - chroffset) - chrstart. */

#if 0
/* Need this procedure because we are skipping some oligomers */
static bool
find_shifted_canonical (Chrpos_T leftpos, Chrpos_T rightpos, int querydistance, 
			Chrpos_T (*genome_left_position)(Chrpos_T, Chrpos_T, Univcoord_T, Univcoord_T, bool),
			Chrpos_T (*genome_right_position)(Chrpos_T, Chrpos_T, Univcoord_T, Univcoord_T, bool),
			Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp, bool skip_repetitive_p) {
  Chrpos_T leftdi, rightdi;
  Chrpos_T last_leftpos, last_rightpos;
  int shift, leftmiss, rightmiss;
  Chrpos_T left_chrbound, right_chrbound;
  
  /* leftpos = prevposition + querydistance + indexsize_nt - 1; */
  /* rightpos = position; */

  debug7(printf("Looking for shifted canonical at leftpos %u to rightpos %u, chroffset %u, chrhigh %u\n",
		leftpos,rightpos,chroffset,chrhigh));

#if 0
  /* previously checked against genomiclength */
  if (leftpos > genomiclength || rightpos > genomiclength) {
    return false;
  }
#else
  /* Checking just before call to genome_right_position */
#endif

  if (leftpos >= rightpos) {
    debug7(printf("leftpos %u >= rightpos %u, so returning false\n",leftpos,rightpos));
    return false;
  }

  if (leftpos < 103) {
    left_chrbound = 3;	/* Previously 0, but then can find splice site at beginning of segment */
  } else {
    left_chrbound = leftpos - 100;
  }

  if (rightpos < 103) {
    right_chrbound = 3;	/* Previously 0, but then can find splice site at beginning of segment */
  } else {
    right_chrbound = rightpos - 100;
  }

#if 0
  if (skip_repetitive_p == false) {

    last_leftpos = (*genome_left_position)(leftpos,left_chrbound,chroffset,chrhigh,plusp);
    last_rightpos = (*genome_right_position)(rightpos,right_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_leftpos %u, last_rightpos %u\n",last_leftpos,last_rightpos));

    debug7(printf("skip_repetitive_p == false, so returning %u == %u && %u == %u\n",
		  leftpos,last_leftpos,rightpos,last_rightpos));
    return (leftpos == last_leftpos && rightpos == last_rightpos);
  }
#endif

  /* Allow canonical to be to right of match */
  leftpos += SHIFT_EXTRA;
  if (leftpos > chrhigh - 3) {
    leftpos = chrhigh - 3;
  }
  rightpos += SHIFT_EXTRA;
  if (rightpos > chrhigh - 3) {
    rightpos = chrhigh - 3;
  }
  debug7(printf("after shift, leftpos = %u, rightpos = %u\n",leftpos,rightpos));

  shift = 0;
  while (shift <= querydistance + SHIFT_EXTRA + SHIFT_EXTRA) {

#if 0
    if (leftpos < 0) {
      return false;
    } else if (rightpos < 0) {
      /* Shouldn't need to check if leftpos >= 0 and rightpos >= leftpos, in the other two conditions) */
      return false;
    } else if (rightpos >= chrlength) {
      return false;
    }
#endif
    if (leftpos < 3) {
      return false;
    } else if (leftpos > rightpos) {
      return false;
    }

    last_leftpos = (*genome_left_position)(leftpos,left_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_leftpos %u\n",last_leftpos));
    assert(last_leftpos != 0U);
    if ((leftdi = last_leftpos) == -1U) {
      debug7(printf("\n"));
      return false;
    } else {
      leftmiss = (int) (leftpos - leftdi);
    }

    last_rightpos = (*genome_right_position)(rightpos,right_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_rightpos %u\n",last_rightpos));
    assert(last_rightpos != 0U);
    if ((rightdi = last_rightpos) == -1U) {
      debug7(printf("\n"));
      return false;
    } else {
      rightmiss = (int) (rightpos - rightdi);
    }

    debug7(printf("shift %d/left %d (miss %d)/right %d (miss %d)\n",shift,leftpos,leftmiss,rightpos,rightmiss));
    if (leftmiss == rightmiss) {  /* was leftmiss == 0 && rightmiss == 0, which doesn't allow for a shift */
      debug7(printf(" => Success at %u..%u (fwd) or %u..%u (rev)\n\n",
		    leftpos-leftmiss+/*onebasedp*/1U,rightpos-rightmiss+/*onebasedp*/1U,
		    chrhigh-chroffset-(leftpos-leftmiss),chrhigh-chroffset-(rightpos-rightmiss)));
      return true;
    } else if (leftmiss >= rightmiss) {
      shift += leftmiss;
      leftpos -= leftmiss;
      rightpos -= leftmiss;
    } else {
      shift += rightmiss;
      leftpos -= rightmiss;
      rightpos -= rightmiss;
    }
  }

  debug7(printf("\n"));
  return false;
}
#endif


#if 0
/* General case for ranges in score_querypos */
while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
  prevlink = &(links[prev_querypos][prevhit]);
  
  gendistance = position - prevposition - indexsize_nt;
  diffdistance = abs(gendistance - querydistance);

  if (diffdistance < maxintronlen) {
    if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
      debug9(canonicalsgn = 9);
      fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
#ifdef PMAP
      if (diffdistance % 3 != 0) {
	fwd_score -= NONCODON_INDEL_PENALTY;
      }
#endif
    } else if (near_end_p == false && prevlink->fwd_consecutive < EXON_DEFN) {
      debug9(canonicalsgn = 0);
      if (splicingp == true) {
	fwd_score = prevlink->fwd_score - (diffdistance/TEN_THOUSAND + 1) - querydist_penalty - NINTRON_PENALTY_MISMATCH;
      } else {
	fwd_score = prevlink->fwd_score - (diffdistance/ONE + 1) - querydist_penalty - NINTRON_PENALTY_MISMATCH;
      }

    } else if (splicingp == false) {
      debug9(canonicalsgn = 0);
      fwd_score = prevlink->fwd_score - (diffdistance/ONE + 1) - querydist_penalty;

    } else if (use_shifted_canonical_p == true) {
      leftpos = prevposition + querydistance - 1;
      /* printf("leftpos %d, last_leftpos %d, rightpos %d\n",leftpos,last_leftpos,rightpos); */
      if (leftpos == last_leftpos) {
	canonicalp = last_canonicalp;
      } else {
	debug7(printf("Calling find_shift_canonical fwd\n"));
	canonicalp = find_shifted_canonical(leftpos,rightpos,querydistance-indexsize_nt,
					    /* &lastGT,&lastAG, */
					    Genome_prev_donor_position,Genome_prev_acceptor_position,
					    chroffset,chrhigh,plusp,skip_repetitive_p);
	/* And need to check for shift_canonical_rev */

	last_leftpos = leftpos;
	last_canonicalp = canonicalp;
      }
      if (canonicalp == true) {
	debug9(canonicalsgn = +1);
	fwd_score = prevlink->fwd_score - (diffdistance/TEN_THOUSAND + 1) - querydist_penalty;
      } else {
	debug9(canonicalsgn = 0);
	fwd_score = prevlink->fwd_score - (diffdistance/TEN_THOUSAND + 1) - querydist_penalty - NINTRON_PENALTY_MISMATCH;
      }

    } else {
      debug9(canonicalsgn = +1);
      fwd_score = prevlink->fwd_score - (diffdistance/TEN_THOUSAND + 1) - querydist_penalty;
    }

    debug9(printf("\tD. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
		  prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		  prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
		  best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		  gendistance,querydistance,canonicalsgn));
	    
    /* Allow ties, which should favor shorter intron */
    if (fwd_score >= best_fwd_score) {
      if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	best_fwd_consecutive = prevlink->fwd_consecutive + (querydistance + indexsize_nt);
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
      } else {
	best_fwd_consecutive = 0;
	/* best_fwd_rootnlinks = 1; */
      }
      best_fwd_score = fwd_score;
      best_fwd_prevpos = prev_querypos;
      best_fwd_prevhit = prevhit;
#ifdef DEBUG9
      best_fwd_tracei = ++*fwd_tracei;
      best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
      best_fwd_intronnrev = prevlink->fwd_intronnrev;
      best_fwd_intronnunk = prevlink->fwd_intronnunk;
      switch (canonicalsgn) {
      case 1: best_fwd_intronnfwd++; break;
      case 0: best_fwd_intronnunk++; break;
      }
#endif
      debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
    } else {
      debug9(printf(" => Loses to %d\n",best_fwd_score));
    }
  }

  prevhit = active[prev_querypos][prevhit];
 }
#endif


static void
score_querypos_lookback (
#ifdef DEBUG9
			 int *fwd_tracei,
#endif
			 Link_T currlink, int querypos,
			 int querystart, int queryend, unsigned int position,
			 struct Link_T **links, Chrpos_T **mappings,
			 int **active, int *firstactive,
			 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			 int indexsize, Intlist_T processed, int sufflookback, int nsufflookback, int maxintronlen, 
			 bool anchoredp, bool localp, bool splicingp, bool skip_repetitive_p,
			 bool use_canonical_p, int non_canonical_penalty) {
  Link_T prevlink;
  Intlist_T p;

  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
  int best_fwd_rootposition = position;
  /* int best_fwd_rootnlinks = 1; */
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
#ifdef DEBUG9
  int best_fwd_tracei;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
  int canonicalsgn;
#endif
  bool adjacentp = false, donep;
  int prev_querypos, prevhit;
  Chrpos_T prevposition, gendistance;
  Chrpos_T leftpos, last_leftpos, rightpos;
  Univcoord_T donorpos, acceptorpos;
  Univcoord_T prevpos, currpos;
  int querydistance, diffdistance, lookback, nlookback, nseen, indexsize_nt;
  /* int querydist_penalty; */
  int querydist_credit;
  int enough_consecutive;
  /* bool near_end_p; */
  bool canonicalp;
  double donor_score, acceptor_score, antidonor_score, antiacceptor_score;

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  enough_consecutive = 32;

#if 0
  if (querypos < querystart + NEAR_END_LENGTH) {
    near_end_p = true;
  } else if (querypos > queryend - NEAR_END_LENGTH) {
    near_end_p = true;
  } else {
    near_end_p = false;
  }
#endif

  /* A. Evaluate adjacent position (at last one processed) */
  if (processed != NULL) {
    prev_querypos = Intlist_head(processed);
#ifdef PMAP
    querydistance = (querypos - prev_querypos)*3;
#else
    querydistance = querypos - prev_querypos;
#endif
    prevhit = firstactive[prev_querypos];
    while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + querydistance <= position) {
      if (prevposition + querydistance == position) {
	prevlink = &(links[prev_querypos][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	best_fwd_rootposition = prevlink->fwd_rootposition;
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	best_fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH*querydistance;

	best_fwd_prevpos = prev_querypos;
	best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	best_fwd_tracei = prevlink->fwd_tracei;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
#endif
	adjacentp = true;

	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,
		      best_fwd_score,best_fwd_consecutive,best_fwd_tracei,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
	prevhit = -1;		/* Exit loop */
      } else {
	prevhit = active[prev_querypos][prevhit];
      }
    }

  }

  if (best_fwd_consecutive < enough_consecutive) {

    /* D. Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Look not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
    } else {
      /* Look farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
    }

    last_leftpos = 0;		/* if use_shifted_canonical_p is true */
    rightpos = position;	/* if use_shifted_canonical_p is true */

    donep = false;
    nseen = 0;

    p = processed;
    if (anchoredp && querypos - indexsize_nt <= querystart) {
      /* Allow close prevpositions that overlap with anchor */
      /* Can give rise to false positives, and increases amount of dynamic programming work */
    } else if (0 && anchoredp && querypos == queryend) {
      /* Test first position */
    } else {
      while (p != NULL && (prev_querypos = Intlist_head(p)) > querypos - indexsize_nt) {
	debug9(printf("Skipping prev_querypos %d, because too close\n",prev_querypos));
	p = Intlist_next(p);
      }
    }

    for ( ; p != NULL && best_fwd_consecutive < enough_consecutive && donep == false;
	 p = Intlist_next(p), nseen++) {
      prev_querypos = Intlist_head(p);

#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3;
#else
      querydistance = querypos - prev_querypos;
#endif
      if (querydistance <= indexsize_nt) {
	querydist_credit = CONSEC_POINTS_PER_MATCH*querydistance;
      } else {
	querydist_credit = CONSEC_POINTS_PER_MATCH*indexsize_nt - querydistance/QUERYDIST_PENALTY_FACTOR;
      }
      if (nseen > nlookback && querydistance - indexsize_nt > lookback) {
	donep = true;
      }

      if ((prevhit = firstactive[prev_querypos]) != -1) {
	/* querydist_penalty = (querydistance - indexsize_nt)/QUERYDIST_PENALTY_FACTOR; */

	/* Range 1: From Infinity to maxintronlen */
	if (splicingp == true) {
	  /* This is equivalent to diffdistance >= maxintronlen, where
	     diffdistance = abs(gendistance - querydistance) and
	     gendistance = (position - prevposition - indexsize_nt) */
	  while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + maxintronlen + querydistance <= position) {
	    /* Skip */
	    /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	    prevhit = active[prev_querypos][prevhit];
	  }
	}

	/* Range 2: From maxintronlen to (prev_querypos + EQUAL_DISTANCE_NOT_SPLICING) */
	/* This is equivalent to +diffdistance > EQUAL_DISTANCE_NOT_SPLICING */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + EQUAL_DISTANCE_NOT_SPLICING + querydistance < position) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = position - prevposition;
	  assert(gendistance > querydistance); /* True because gendistance > EQUAL_DISTANCE_NOT_SPLICING + querydistance */
	  diffdistance = gendistance - querydistance; /* No need for abs() */

	  fwd_score = prevlink->fwd_score + querydist_credit /*- querydist_penalty*/;
	  if (splicingp == true) {
	    fwd_score -= (diffdistance/TEN_THOUSAND + 1);
	  } else {
	    fwd_score -= (diffdistance/ONE + 1);
	  }

	  if (/*near_end_p == false &&*/ prevlink->fwd_consecutive < EXON_DEFN) {
	    debug9(canonicalsgn = 0);
	    fwd_score -= NINTRON_PENALTY_MISMATCH;

	  } else if (use_canonical_p == true) {

	    /* prevpos is lower genomic coordinate than currpos */
	    /* need to subtract from position and prevposition to compensate for greedy matches */
	    /* need to add to position and prevposition to compensate for missed matches */
	    if (plusp == true) {
	      prevpos = chroffset + prevposition + indexsize_nt;
	      currpos = chroffset + position - querydistance + indexsize_nt;
#if 0
	      donor_score = Maxent_hr_donor_prob(prevpos,chroffset);
	      acceptor_score = Maxent_hr_acceptor_prob(currpos,chroffset);
	      antidonor_score = Maxent_hr_antidonor_prob(currpos,chroffset);
	      antiacceptor_score = Maxent_hr_antiacceptor_prob(prevpos,chroffset);
	      debug9(printf("lookback plus: sense donorpos %u, acceptorpos %u: %f..%f\n",
			    prevpos-chroffset,currpos-chroffset+1,donor_score,acceptor_score));
	      debug9(printf("lookback plus: anti donorpos %u, acceptorpos %u: %f..%f\n",
			    currpos-chroffset+1,prevpos-chroffset,antidonor_score,antiacceptor_score));
#elif 0
	      donorpos = Genome_prev_donor_position(/*right*/prevpos + MISS_BEHIND,
						    /*left*/prevpos - GREEDY_ADVANCE,chroffset);
	      acceptorpos = Genome_prev_acceptor_position(/*right*/currpos + MISS_BEHIND,
							  /*left*/currpos - GREEDY_ADVANCE,chroffset);
	      debug9(printf("lookback plus: sense prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));

	      donorpos = Genome_prev_antidonor_position(/*right*/currpos + MISS_BEHIND,
							/*left*/currpos - GREEDY_ADVANCE,chroffset);
	      acceptorpos = Genome_prev_antiacceptor_position(/*right*/prevpos + MISS_BEHIND,
							      /*left*/prevpos - GREEDY_ADVANCE,chroffset);
	      debug9(printf("lookback plus: anti prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));
#else
	      if (prevpos < GREEDY_ADVANCE || currpos < GREEDY_ADVANCE) {
		canonicalp = false;
	      } else if (Genome_sense_canonicalp(/*donor_rightbound*/prevpos + MISS_BEHIND,
						 /*donor_leftbound*/prevpos - GREEDY_ADVANCE,
						 /*acceptor_rightbound*/currpos + MISS_BEHIND,
						 /*acceptor_leftbound*/currpos - GREEDY_ADVANCE,
						 chroffset) == true) {
		debug9(printf("lookback plus: sense canonical\n"));
		canonicalp = true;
	      } else if (Genome_antisense_canonicalp(/*donor_rightbound*/currpos + MISS_BEHIND,
						     /*donor_leftbound*/currpos - GREEDY_ADVANCE,
						     /*acceptor_rightbound*/prevpos + MISS_BEHIND,
						     /*acceptor_leftbound*/prevpos - GREEDY_ADVANCE,
						     chroffset) == true) {
		debug9(printf("lookback plus: antisense canonical\n"));
		canonicalp = true;
	      } else {
		debug9(printf("lookback plus: not canonical\n"));
		canonicalp = false;
	      }
#endif

	    } else {
	      prevpos = chrhigh + 1 - prevposition - indexsize_nt;
	      currpos = chrhigh + 1 - position + querydistance - indexsize_nt;
#if 0
	      donor_score = Maxent_hr_donor_prob(currpos,chroffset);
	      acceptor_score = Maxent_hr_acceptor_prob(prevpos,chroffset);
	      antidonor_score = Maxent_hr_antidonor_prob(prevpos,chroffset);
	      antiacceptor_score = Maxent_hr_antiacceptor_prob(currpos,chroffset);
	      debug9(printf("lookback minus: sense donorpos %u, acceptorpos %u: %f..%f\n",
			    currpos-chroffset,prevpos-chroffset+1,donor_score,acceptor_score));
	      debug9(printf("lookback minus: anti donorpos %u, acceptorpos %u: %f..%f\n",
			    prevpos-chroffset+1,currpos-chroffset,antidonor_score,antiacceptor_score));
#elif 0
	      donorpos = Genome_prev_donor_position(/*right*/currpos + GREEDY_ADVANCE,
						    /*left*/currpos - MISS_BEHIND,chroffset);
	      acceptorpos = Genome_prev_acceptor_position(/*right*/prevpos + GREEDY_ADVANCE,
							  /*left*/prevpos - MISS_BEHIND,chroffset);
	      debug9(printf("lookback minus: sense prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));

	      donorpos = Genome_prev_antidonor_position(/*right*/prevpos + GREEDY_ADVANCE,
							/*left*/prevpos - MISS_BEHIND,chroffset);
	      acceptorpos = Genome_prev_antiacceptor_position(/*right*/currpos + GREEDY_ADVANCE,
							      /*left*/currpos - MISS_BEHIND,chroffset);
	      debug9(printf("lookback minus: anti prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));
#else
	      if (currpos < MISS_BEHIND || prevpos < MISS_BEHIND) {
		canonicalp = false;
	      } else if (Genome_sense_canonicalp(/*donor_rightbound*/currpos + GREEDY_ADVANCE,
						 /*donor_leftbound*/currpos - MISS_BEHIND,
						 /*acceptor_rightbound*/prevpos + GREEDY_ADVANCE,
						 /*acceptor_leftbound*/prevpos - MISS_BEHIND,
						 chroffset) == true) {
		debug9(printf("lookback minus: sense canonical\n"));
		canonicalp = true;
	      } else if (Genome_antisense_canonicalp(/*donor_rightbound*/prevpos + GREEDY_ADVANCE,
						     /*donor_leftbound*/prevpos - MISS_BEHIND,
						     /*acceptor_rightbound*/currpos + GREEDY_ADVANCE,
						     /*acceptor_leftbound*/currpos - MISS_BEHIND,
						     chroffset) == true) {
		debug9(printf("lookback minus: antisense canonical\n"));
		canonicalp = true;
	      } else {
		debug9(printf("lookback minus: not canonical\n"));
		canonicalp = false;
	      }
#endif
	    }

	    if (canonicalp == false) {
	      fwd_score -= non_canonical_penalty;
	    }

	  }

	  debug9(printf("\tD2. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,canonicalsgn));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    best_fwd_tracei = ++*fwd_tracei;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}

#if 0
	/* Scoring appears to be the same as for range 4, which is rarely called */
	/* Range 3: From (querypos + EQUAL_DISTANCE_NOT_SPLICING) to (querypos - EQUAL_DISTANCE_NOT_SPLICING) */
	/* This is equivalent to -diffdistance > EQUAL_DISTANCE_NOT_SPLICING && prevposition + indexsize_nt <= position */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + querydistance < position + EQUAL_DISTANCE_NOT_SPLICING &&
	       prevposition + indexsize_nt <= position) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = position - prevposition;
	  diffdistance = abs(gendistance - querydistance);

	  /* canonicalsgn = 9; */
	  /* ? Add diffdistance as penalty here */
	  fwd_score = prevlink->fwd_score + querydist_credit - (diffdistance/ONE + 1);
#ifdef PMAP
	  if (diffdistance % 3 != 0) {
	    fwd_score -= NONCODON_INDEL_PENALTY;
	  }
#endif
	  debug9(printf("\tD3. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,/*canonicalsgn*/9));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    best_fwd_tracei = prevlink->fwd_tracei; /* Keep previous trace */
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
#if 0
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}
#endif

	/* Range 4: From (prev_querypos - EQUAL_DISTANCE_NOT_SPLICING) to indexsize_nt */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = position - prevposition;
	  diffdistance = abs(gendistance - querydistance);

	  fwd_score = prevlink->fwd_score + querydist_credit - (diffdistance/ONE + 1) /*- querydist_penalty*/;
#if 0
	  /* Used in range 4 but not in range 3 */
	  if (/*near_end_p == false &&*/ prevlink->fwd_consecutive < EXON_DEFN) {
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }
#endif

	  debug9(printf("\tD4. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,canonicalsgn));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    /* best_fwd_tracei = ++*fwd_tracei; */
	    best_fwd_tracei = prevlink->fwd_tracei; /* Keep previous trace, as in range 3 */
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}
      }

    }
  }

  /* Best_score needs to beat something positive to prevent a
     small local extension from beating a good canonical intron.
     If querypos is too small, don't insert an intron.  */
  /* linksconsecutive already assigned above */
  currlink->fwd_consecutive = best_fwd_consecutive;
  currlink->fwd_rootposition = best_fwd_rootposition;
  /* currlink->fwd_rootnlinks = best_fwd_rootnlinks; */
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  if (currlink->fwd_pos >= 0) {
    debug9(currlink->fwd_tracei = best_fwd_tracei);
    currlink->fwd_score = best_fwd_score;
  } else if (anchoredp == true) {
    debug9(currlink->fwd_tracei = -1);
    currlink->fwd_score = -100000;
  } else if (localp == true) {
    debug9(currlink->fwd_tracei = ++*fwd_tracei);
    currlink->fwd_score = indexsize_nt;
  } else {
    debug9(currlink->fwd_tracei = ++*fwd_tracei);
    currlink->fwd_score = best_fwd_score;
  }

#ifdef DEBUG9
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;
#endif

  debug9(printf("\tChose %d,%d with score %d (fwd) => trace #%d\n",
		currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score,currlink->fwd_tracei));
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));

  return;
}


static void
score_querypos_lookforward (
#ifdef DEBUG9
			    int *fwd_tracei,
#endif
			    Link_T currlink, int querypos,
			    int querystart, int queryend, unsigned int position,
			    struct Link_T **links, Chrpos_T **mappings,
			    int **active, int *firstactive,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			    int indexsize, Intlist_T processed, int sufflookback, int nsufflookback, int maxintronlen, 
			    bool anchoredp, bool localp, bool splicingp, bool skip_repetitive_p,
			    bool use_canonical_p, int non_canonical_penalty) {
  Link_T prevlink;
  Intlist_T p;

  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
  int best_fwd_rootposition = position;
  /* int best_fwd_rootnlinks = 1; */
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
#ifdef DEBUG9
  int best_fwd_tracei;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
  int canonicalsgn;
#endif
  bool adjacentp = false, donep;
  int prev_querypos, prevhit;
  Chrpos_T prevposition, gendistance;
  Chrpos_T leftpos, rightpos, last_rightpos;
  Univcoord_T prevpos, currpos;
  Univcoord_T donorpos, acceptorpos;
  int querydistance, diffdistance, lookback, nlookback, nseen, indexsize_nt;
  /* int querydist_penalty; */
  int querydist_credit;
  int enough_consecutive;
  /* bool near_end_p; */
  bool canonicalp;
  double donor_score, acceptor_score, antidonor_score, antiacceptor_score;

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  enough_consecutive = 32;

#if 0
  if (querypos < querystart + NEAR_END_LENGTH) {
    near_end_p = true;
  } else if (querypos > queryend - NEAR_END_LENGTH) {
    near_end_p = true;
  } else {
    near_end_p = false;
  }
#endif

  /* A. Evaluate adjacent position (at last one processed) */
  if (processed != NULL) {
    prev_querypos = Intlist_head(processed);
#ifdef PMAP
    querydistance = (prev_querypos - querypos)*3;
#else
    querydistance = prev_querypos - querypos;
#endif
    prevhit = firstactive[prev_querypos];
    while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) >= position + querydistance) {
      if (prevposition == position + querydistance) {
	prevlink = &(links[prev_querypos][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	best_fwd_rootposition = prevlink->fwd_rootposition;
	best_fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH*querydistance;

	best_fwd_prevpos = prev_querypos;
	best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	best_fwd_tracei = prevlink->fwd_tracei;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
#endif
	adjacentp = true;

	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,
		      best_fwd_score,best_fwd_consecutive,best_fwd_tracei,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
	prevhit = -1;		/* Exit loop */
      } else {
	prevhit = active[prev_querypos][prevhit];
      }
    }

  }

  if (best_fwd_consecutive < enough_consecutive) {

    /* D. Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Look not so far forward */
      nlookback = 1;
      lookback = sufflookback/2;
    } else {
      /* Look farther forward */
      nlookback = nsufflookback;
      lookback = sufflookback;
    }

    last_rightpos = 0;		/* if use_shifted_canonical_p is true */
    leftpos = position;	      /* if use_shifted_canonical_p is true */

    donep = false;
    nseen = 0; 

    p = processed;
    if (anchoredp && querypos + indexsize_nt >= queryend) {
      /* Allow close prevpositions that overlap with anchor */
      /* Can give rise to false positives, and increases amount of dynamic programming work */
      debug9(printf("No skipping because close to anchor\n"));
    } else if (0 && anchoredp && querypos == querystart) {
      /* Test end position */
    } else {
      while (p != NULL && (prev_querypos = Intlist_head(p)) < querypos + indexsize_nt) {
	debug9(printf("Skipping prev_querypos %d, because too close\n",prev_querypos));
	p = Intlist_next(p);
      }
    }

    for ( ; p != NULL && best_fwd_consecutive < enough_consecutive && donep == false;
	  p = Intlist_next(p), nseen++) {
      prev_querypos = Intlist_head(p);
				  
#ifdef PMAP
      querydistance = (prev_querypos - querypos)*3;
#else
      querydistance = prev_querypos - querypos;
#endif
      /* Allow close prevpositions at beginning of path, to overlap with anchor */
      if (querydistance <= indexsize_nt) {
	querydist_credit = CONSEC_POINTS_PER_MATCH*querydistance;
      } else {
	querydist_credit = CONSEC_POINTS_PER_MATCH*indexsize_nt - querydistance/QUERYDIST_PENALTY_FACTOR;
      }
      if (nseen > nlookback && querydistance - indexsize_nt > lookback) {
	donep = true;
      }

      if ((prevhit = firstactive[prev_querypos]) != -1) {
	/* querydist_penalty = (querydistance - indexsize_nt)/QUERYDIST_PENALTY_FACTOR; */

	/* Range 1: From Infinity to maxintronlen */
	if (splicingp == true) {
	  /* This is equivalent to diffdistance >= maxintronlen, where
	     diffdistance = abs(gendistance - querydistance) and
	     gendistance = (position - prevposition - indexsize_nt) */
	  while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) >= position + maxintronlen + querydistance) {
	    /* Skip */
	    /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	    prevhit = active[prev_querypos][prevhit];
	  }
	}

	/* Range 2: From maxintronlen to (prev_querypos + EQUAL_DISTANCE_NOT_SPLICING) */
	/* This is equivalent to +diffdistance > EQUAL_DISTANCE_NOT_SPLICING */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) > position + EQUAL_DISTANCE_NOT_SPLICING + querydistance) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = prevposition - position;
	  assert(gendistance > querydistance); /* True because gendistance > EQUAL_DISTANCE_NOT_SPLICING + querydistance */
	  diffdistance = gendistance - querydistance; /* No need for abs() */

	  fwd_score = prevlink->fwd_score + querydist_credit /*- querydist_penalty*/;
	  if (splicingp == true) {
	    fwd_score -= (diffdistance/TEN_THOUSAND + 1);
	  } else {
	    fwd_score -= (diffdistance/ONE + 1);
	  }

	  if (/*near_end_p == false &&*/ prevlink->fwd_consecutive < EXON_DEFN) {
	    debug9(canonicalsgn = 0);
	    fwd_score -= NINTRON_PENALTY_MISMATCH;

	  } else if (use_canonical_p == true) {

	    /* prevpos is higher genomic coordinate than currpos */
	    /* need to add to position and prevposition to compensate for greedy matches */
	    /* need to subtract from position and prevposition to compensate for missed matches */
	    if (plusp == true) {
	      prevpos = chroffset + prevposition;
	      currpos = chroffset + position + querydistance;
#if 0
	      donor_score = Maxent_hr_donor_prob(currpos,chroffset);
	      acceptor_score = Maxent_hr_acceptor_prob(prevpos,chroffset);
	      antidonor_score = Maxent_hr_antidonor_prob(prevpos,chroffset);
	      antiacceptor_score = Maxent_hr_antiacceptor_prob(currpos,chroffset);
	      debug9(printf("lookforward plus: sense donorpos %u, acceptorpos %u: %f..%f\n",
			    currpos-chroffset,prevpos-chroffset+1,donor_score,acceptor_score));
	      debug9(printf("lookforward plus: anti donorpos %u, acceptorpos %u: %f..%f\n",
			    prevpos-chroffset+1,currpos-chroffset,antidonor_score,antiacceptor_score));
#elif 0
	      donorpos = Genome_prev_donor_position(/*right*/currpos + GREEDY_ADVANCE,
						    /*left*/currpos - MISS_BEHIND,chroffset);
	      acceptorpos = Genome_prev_acceptor_position(/*right*/prevpos + GREEDY_ADVANCE,
							  /*left*/prevpos - MISS_BEHIND,chroffset);
	      debug9(printf("lookforward plus: sense prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));

	      donorpos = Genome_prev_antidonor_position(/*right*/prevpos + GREEDY_ADVANCE,
							/*left*/prevpos - MISS_BEHIND,chroffset);
	      acceptorpos = Genome_prev_antiacceptor_position(/*right*/currpos + GREEDY_ADVANCE,
							      /*left*/currpos - MISS_BEHIND,chroffset);
	      debug9(printf("lookforward plus: anti prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));
#else
	      if (currpos < MISS_BEHIND || prevpos < MISS_BEHIND) {
		canonicalp = false;
	      } else if (Genome_sense_canonicalp(/*donor_rightbound*/currpos + GREEDY_ADVANCE,
						 /*donor_leftbound*/currpos - MISS_BEHIND,
						 /*acceptor_rightbound*/prevpos + GREEDY_ADVANCE,
						 /*acceptor_leftbound*/prevpos - MISS_BEHIND,
						 chroffset) == true) {
		debug9(printf("lookforward plus: sense canonical\n"));
		canonicalp = true;
	      } else if (Genome_antisense_canonicalp(/*donor_rightbound*/prevpos + GREEDY_ADVANCE,
						     /*donor_leftbound*/prevpos - MISS_BEHIND,
						     /*acceptor_rightbound*/currpos + GREEDY_ADVANCE,
						     /*acceptor_leftbound*/currpos - MISS_BEHIND,
						     chroffset) == true) {
		debug9(printf("lookforward plus: antisense canonical\n"));
		canonicalp = true;
	      } else {
		debug9(printf("lookforward plus: not canonical\n"));
		canonicalp = false;
	      }
#endif
	      
	    } else {
	      prevpos = chrhigh + 1 - prevposition;
	      currpos = chrhigh + 1 - position - querydistance;
#if 0
	      donor_score = Maxent_hr_donor_prob(prevpos,chroffset);
	      acceptor_score = Maxent_hr_acceptor_prob(currpos,chroffset);
	      antidonor_score = Maxent_hr_antidonor_prob(currpos,chroffset);
	      antiacceptor_score = Maxent_hr_antiacceptor_prob(prevpos,chroffset);
	      debug9(printf("lookforward minus: sense donorpos %u, acceptorpos %u: %f..%f\n",
			    prevpos-chroffset,currpos-chroffset+1,donor_score,acceptor_score));
	      debug9(printf("lookforward minus: anti donorpos %u, acceptorpos %u: %f..%f\n",
			    currpos-chroffset+1,prevpos-chroffset,antidonor_score,antiacceptor_score)); 
#elif 0
	      donorpos = Genome_prev_donor_position(/*right*/prevpos + MISS_BEHIND,
						    /*left*/prevpos - GREEDY_ADVANCE,chroffset);
	      acceptorpos = Genome_prev_acceptor_position(/*right*/currpos + MISS_BEHIND,
							  /*left*/currpos - GREEDY_ADVANCE,chroffset);
	      debug9(printf("lookforward minus: sense prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));

	      donorpos = Genome_prev_antidonor_position(/*right*/currpos + MISS_BEHIND,
							/*left*/currpos - GREEDY_ADVANCE,chroffset);
	      acceptorpos = Genome_prev_antiacceptor_position(/*right*/prevpos + MISS_BEHIND,
							      /*left*/prevpos - GREEDY_ADVANCE,chroffset);
	      debug9(printf("lookforward minus: anti prev donorpos %u, prev acceptorpos %u\n",donorpos-chroffset,acceptorpos-chroffset));
#else
	      if (prevpos < GREEDY_ADVANCE || currpos < GREEDY_ADVANCE) {
		canonicalp = false;
	      } else if (Genome_sense_canonicalp(/*donor_rightbound*/prevpos + MISS_BEHIND,
						 /*donor_leftbound*/prevpos - GREEDY_ADVANCE,
						 /*acceptor_rightbound*/currpos + MISS_BEHIND,
						 /*acceptor_leftbound*/currpos - GREEDY_ADVANCE,
						 chroffset) == true) {
		debug9(printf("lookforward minus: sense canonical\n"));
		canonicalp = true;
	      } else if (Genome_antisense_canonicalp(/*donor_rightbound*/currpos + MISS_BEHIND,
						     /*donor_leftbound*/currpos - GREEDY_ADVANCE,
						     /*acceptor_rightbound*/prevpos + MISS_BEHIND,
						     /*acceptor_leftbound*/prevpos - GREEDY_ADVANCE,
						     chroffset) == true) {
		debug9(printf("lookforward minus: antisense canonical\n"));
		canonicalp = true;
	      } else {
		debug9(printf("lookforward minus: not canonical\n"));
		canonicalp = false;
	      }
#endif
	    }

	    if (canonicalp == false) {
	      fwd_score -= non_canonical_penalty;
	    }
	  }

	  debug9(printf("\tD2. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,canonicalsgn));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    best_fwd_tracei = ++*fwd_tracei;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}

#if 0
	/* Scoring appears to be the same as for range 4, which is rarely called */
	/* Range 3: From (querypos + EQUAL_DISTANCE_NOT_SPLICING) to (querypos - EQUAL_DISTANCE_NOT_SPLICING) */
	/* This is equivalent to -diffdistance > EQUAL_DISTANCE_NOT_SPLICING && prevposition + indexsize_nt <= position */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + EQUAL_DISTANCE_NOT_SPLICING > position + querydistance &&
	       prevposition >= position + indexsize_nt) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = prevposition - position;
	  diffdistance = abs(gendistance - querydistance);

	  /* canonicalsgn = 9; */
	  /* ? Add diffdistance as penalty here */
	  fwd_score = prevlink->fwd_score + querydist_credit - (diffdistance/ONE + 1);
#ifdef PMAP
	  if (diffdistance % 3 != 0) {
	    fwd_score -= NONCODON_INDEL_PENALTY;
	  }
#endif
	  debug9(printf("\tD3. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,/*canonicalsgn*/9));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    best_fwd_tracei = prevlink->fwd_tracei; /* Keep previous trace */
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
#if 0
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}
#endif

	/* Range 4: From (prev_querypos - EQUAL_DISTANCE_NOT_SPLICING) to indexsize_nt */
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) >= position + indexsize_nt) {
	  /* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = prevposition - position;
	  diffdistance = abs(gendistance - querydistance);

	  fwd_score = prevlink->fwd_score + querydist_credit - (diffdistance/ONE + 1) /*- querydist_penalty*/;
#if 0
	  if (/*near_end_p == false &&*/ prevlink->fwd_consecutive < EXON_DEFN) {
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }
#endif

	  debug9(printf("\tD4. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from #%d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_tracei,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance-indexsize_nt,querydistance-indexsize_nt,canonicalsgn));
	    
	  /* Disallow ties, which should favor adjacent */
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_rootposition = prevlink->fwd_rootposition;
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
#ifdef DEBUG9
	    /* best_fwd_tracei = ++*fwd_tracei; */
	    best_fwd_tracei = prevlink->fwd_tracei; /* Keep previous trace, as in range 3 */
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
#endif
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }

	  prevhit = active[prev_querypos][prevhit];
	}
      }

    }
  }

  /* Best_score needs to beat something positive to prevent a
     small local extension from beating a good canonical intron.
     If querypos is too small, don't insert an intron.  */
  /* linksconsecutive already assigned above */
  currlink->fwd_consecutive = best_fwd_consecutive;
  currlink->fwd_rootposition = best_fwd_rootposition;
  /* currlink->fwd_rootnlinks = best_fwd_rootnlinks; */
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  if (currlink->fwd_pos >= 0) {
    debug9(currlink->fwd_tracei = best_fwd_tracei);
    currlink->fwd_score = best_fwd_score;
  } else if (anchoredp == true) {
    debug9(currlink->fwd_tracei = -1);
    currlink->fwd_score = -100000;
  } else if (localp == true) {
    debug9(currlink->fwd_tracei = ++*fwd_tracei);
    currlink->fwd_score = indexsize_nt;
  } else {
    debug9(currlink->fwd_tracei = ++*fwd_tracei);
    currlink->fwd_score = best_fwd_score;
  }

#ifdef DEBUG9
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;
#endif

  debug9(printf("\tChose %d,%d with score %d (fwd) => trace #%d\n",
		currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score,currlink->fwd_tracei));
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));

  return;
}


static void
revise_active_lookback (int **active, int *firstactive, int *nactive, 
			int low_hit, int high_hit, struct Link_T **links, int querypos,
			Chrpos_T **mappings) {
  int best_score, threshold, score;
  int hit, *ptr;

  debug6(printf("Revising querypos %d from low_hit %d to high_hit %d.  Scores:\n",querypos,low_hit,high_hit));
  if ((hit = low_hit) >= high_hit) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  } else {
    debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
    best_score = links[querypos][hit].fwd_score;
#ifdef SEPARATE_FWD_REV
    debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
    if ((score = links[querypos][hit].rev_score) > best_score) {
      best_score = score;
    }
#endif
    debug6(printf("\n"));

    for (hit++; hit < high_hit; hit++) {
      debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
      }
#ifdef SEPARATE_FWD_REV
      debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
      }
#endif
      debug6(printf("\n"));
    }

    threshold = best_score - SCORE_FOR_RESTRICT;
    if (threshold < 0) {
      threshold = 0;
    }

    nactive[querypos] = 0;
    ptr = &(firstactive[querypos]);
    hit = low_hit;
    while (hit < high_hit) {
      while (hit < high_hit && links[querypos][hit].fwd_score <= threshold
#ifdef SEPARATE_FWD_REV
	     && links[querypos][hit].rev_score <= threshold
#endif
	     ) {
	hit++;
      }
      *ptr = hit;
      if (hit < high_hit) {
	nactive[querypos] += 1;
	ptr = &(active[querypos][hit]);
	hit++;
      }
    }
    *ptr = -1;
  }

  debug6(
	 printf("Valid hits (%d) at querypos %d (firstactive %d):",nactive[querypos],querypos,firstactive[querypos]);
	 hit = firstactive[querypos];
	 while (hit != -1) {
	   printf(" %d",hit);
	   hit = active[querypos][hit];
	 }
	 printf("\n");
	 );

  return;
}


static void
revise_active_lookforward (int **active, int *firstactive, int *nactive, 
			   int low_hit, int high_hit, struct Link_T **links, int querypos,
			   Chrpos_T **mappings) {
  int best_score, threshold, score;
  int hit, *ptr;

  debug6(printf("Revising querypos %d from high_hit %d to low_hit %d.  Scores:\n",querypos,high_hit,low_hit));
  if ((hit = high_hit - 1) < low_hit) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  } else {
    debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
    best_score = links[querypos][hit].fwd_score;
#ifdef SEPARATE_FWD_REV
    debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
    if ((score = links[querypos][hit].rev_score) > best_score) {
      best_score = score;
    }
#endif
    debug6(printf("\n"));

    for (--hit; hit >= low_hit; --hit) {
      debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
      }
#ifdef SEPARATE_FWD_REV
      debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
      }
#endif
      debug6(printf("\n"));
    }

    threshold = best_score - SCORE_FOR_RESTRICT;
    if (threshold < 0) {
      threshold = 0;
    }

    nactive[querypos] = 0;
    ptr = &(firstactive[querypos]);
    hit = high_hit - 1;
    while (hit >= low_hit) {
      while (hit >= low_hit && links[querypos][hit].fwd_score <= threshold
#ifdef SEPARATE_FWD_REV
	     && links[querypos][hit].rev_score <= threshold
#endif
	     ) {
	--hit;
      }
      *ptr = hit;
      if (hit >= low_hit) {
	nactive[querypos] += 1;
	ptr = &(active[querypos][hit]);
	--hit;
      }
    }
    *ptr = -1;
  }

  debug6(
	 printf("Valid hits (%d) at querypos %d (firstactive %d):",nactive[querypos],querypos,firstactive[querypos]);
	 hit = firstactive[querypos];
	 while (hit != -1) {
	   printf(" %d",hit);
	   hit = active[querypos][hit];
	 }
	 printf("\n");
	 );

  return;
}



static int **
intmatrix_1d_new (int length1, int *lengths2, int totallength) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  matrix[0] = (int *) CALLOC(totallength,sizeof(int));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] <= 0) {
      matrix[i] = matrix[i-1];
    } else {
      matrix[i] = &(matrix[i-1][lengths2[i-1]]);
    }
  }
  return matrix;
}

static void
intmatrix_1d_free (int ***matrix) {
  FREE((*matrix)[0]);
  FREE(*matrix);
  return;
}


static int **
intmatrix_2d_new (int length1, int *lengths2) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      matrix[i] = (int *) NULL;
    } else {
      matrix[i] = (int *) CALLOC(lengths2[i],sizeof(int));
    }
  }
  return matrix;
}

static void
intmatrix_2d_free (int ***matrix, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*matrix)[i]) {
      FREE((*matrix)[i]);
    }
  }
  FREE(*matrix);
  return;
}


/************************************************************************
 *   Cells used for ranking hits
 ************************************************************************/

#if 0
typedef struct Cell_T *Cell_T;
struct Cell_T {
  int rootposition;
  int querypos;
  int hit;
  bool fwdp;
  int score;
};

/* Replaced by Cellpool_T routines */
static void
Cell_free (Cell_T *old) {
  FREE(*old);
  return;
}


static Cell_T
Cell_new (int rootposition, int querypos, int hit, bool fwdp, int score) {
  Cell_T new = (Cell_T) MALLOC(sizeof(*new));

  new->rootposition = rootposition;
  new->querypos = querypos;
  new->hit = hit;
  new->fwdp = fwdp;
  new->score = score;
  return new;
}
#endif


static int
Cell_rootposition_left_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->rootposition < y->rootposition) {
    return -1;
  } else if (y->rootposition < x->rootposition) {
    return +1;
#if 0
  } else if (x->tracei < y->tracei) {
    return -1;
  } else if (y->tracei < x->tracei) {
    return +1;
#endif
  } else if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else if (x->querypos > y->querypos) {
    return -1;
  } else if (y->querypos > x->querypos) {
    return +1;
  } else if (x->hit < y->hit) {
    return -1;
  } else if (y->hit < x->hit) {
    return +1;
  } else if (x->fwdp == true && y->fwdp == false) {
    return -1;
  } else if (y->fwdp == true && x->fwdp == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_rootposition_right_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->rootposition < y->rootposition) {
    return -1;
  } else if (y->rootposition < x->rootposition) {
    return +1;
#if 0
  } else if (x->tracei < y->tracei) {
    return -1;
  } else if (y->tracei < x->tracei) {
    return +1;
#endif
  } else if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else if (x->querypos > y->querypos) {
    return -1;
  } else if (y->querypos > x->querypos) {
    return +1;
  } else if (x->hit > y->hit) {
    return -1;
  } else if (y->hit > x->hit) {
    return +1;
  } else if (x->fwdp == true && y->fwdp == false) {
    return -1;
  } else if (y->fwdp == true && x->fwdp == false) {
    return +1;
  } else {
    return 0;
  }
}


static int
Cell_score_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else {
    return 0;
  }
}


#ifdef USE_THRESHOLD_SCORE
/* Doesn't work well for short dynamic programming at the ends of a read */
static Cell_T *
Linkmatrix_get_cells_fwd (int *nunique, struct Link_T **links, int querystart, int queryend, int *npositions,
			  int indexsize, int bestscore, bool favor_right_p, Cellpool_T cellpool) {
  Cell_T *sorted, *cells;
  List_T celllist = NULL;
  int querypos, hit;
  int rootposition, last_rootposition;
  int threshold_score;
  int ngood, ncells, i, k;

  if (bestscore > 2*suboptimal_score_end) {
    threshold_score = bestscore - suboptimal_score_end;
  } else {
    threshold_score = bestscore/2;
  }
  if (threshold_score <= indexsize) {
    threshold_score = indexsize + 1;
  }

  ncells = 0;
  for (querypos = querystart; querypos <= queryend; querypos++) {
    ngood = 0;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if (links[querypos][hit].fwd_score >= threshold_score) {
	ngood++;
      }
    }
    if (ngood > 0 && ngood <= 10) {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	debug11(printf("  At %d,%d, comparing score %d with threshold_score %d\n",
		       querypos,hit,links[querypos][hit].fwd_score,threshold_score));
	if (links[querypos][hit].fwd_score >= threshold_score) {
	  rootposition = links[querypos][hit].fwd_rootposition;
	  /* tracei = links[querypos][hit].fwd_tracei; */
	  celllist = Cellpool_push(celllist,cellpool,rootposition,querypos,hit,/*fwdp*/true,
				   links[querypos][hit].fwd_score);
	  ncells++;
	}
      }
    }
  }

  if (ncells == 0) {
    *nunique = 0;
    return (Cell_T *) NULL;

  } else {
    /* Take best result for each tracei */
    cells = (Cell_T *) List_to_array(celllist,NULL);
    /* List_free(&celllist); -- No need with cellpool */

    if (favor_right_p == true) {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_right_cmp);
    } else {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_left_cmp);
    }

    sorted = (Cell_T *) CALLOC(ncells,sizeof(Cell_T));
    k = 0;

    last_rootposition = -1;
    for (i = 0; i < ncells; i++) {
      if (cells[i]->rootposition == last_rootposition) {
	/* Cell_free(&(cells[i])); -- no need with cellpool */
      } else {
	debug11(printf("Pushing rootposition %d, trace #%d, score %d, pos %d, hit %d\n",
		       cells[i]->rootposition,cells[i]->tracei,cells[i]->score,cells[i]->querypos,cells[i]->hit));
	sorted[k++] = cells[i];
	last_rootposition = cells[i]->rootposition;
      }
    }
    debug11(printf("\n"));
    FREE(cells);
  
    *nunique = k;
    qsort(sorted,*nunique,sizeof(Cell_T),Cell_score_cmp);

    return sorted;
  }
}

#else

static Cell_T *
Linkmatrix_get_cells_fwd (int *nunique, struct Link_T **links, int querystart, int queryend, int *npositions,
			  int indexsize, int bestscore, bool favor_right_p, Cellpool_T cellpool) {
  Cell_T *sorted, *cells;
  List_T celllist = NULL;
  int querypos, hit;
  int rootposition, last_rootposition;
  int ngood, ncells, i, k;

  ncells = 0;
  for (querypos = querystart; querypos <= queryend; querypos++) {
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if (links[querypos][hit].fwd_score > 0) {
	rootposition = links[querypos][hit].fwd_rootposition;
	/* tracei = links[querypos][hit].fwd_tracei; */
	celllist = Cellpool_push(celllist,cellpool,rootposition,querypos,hit,/*fwdp*/true,
				 links[querypos][hit].fwd_score);
	ncells++;
      }
    }
  }

  if (ncells == 0) {
    *nunique = 0;
    return (Cell_T *) NULL;

  } else {
    /* Take best result for each tracei */
    cells = (Cell_T *) List_to_array(celllist,NULL);
    /* List_free(&celllist); -- No need with cellpool */

    if (favor_right_p == true) {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_right_cmp);
    } else {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_left_cmp);
    }

    sorted = (Cell_T *) CALLOC(ncells,sizeof(Cell_T));
    k = 0;

    last_rootposition = -1;
    for (i = 0; i < ncells; i++) {
      if (cells[i]->rootposition == last_rootposition) {
	/* Cell_free(&(cells[i])); -- no need with cellpool */
      } else {
	debug11(printf("Pushing rootposition %d, score %d, pos %d, hit %d\n",
		       cells[i]->rootposition,cells[i]->score,cells[i]->querypos,cells[i]->hit));
	sorted[k++] = cells[i];
	last_rootposition = cells[i]->rootposition;
      }
    }
    debug11(printf("\n"));
    FREE(cells);
  
    *nunique = k;
    qsort(sorted,*nunique,sizeof(Cell_T),Cell_score_cmp);

    return sorted;
  }
}

#endif

static Cell_T *
Linkmatrix_get_cells_both (int *nunique, struct Link_T **links, int querystart, int queryend, int *npositions,
			   int indexsize, int bestscore, bool favor_right_p, Cellpool_T cellpool) {
  Cell_T *sorted, *cells;
  List_T celllist = NULL;
  int querypos, hit;
  int rootposition, last_rootposition;
  int threshold_score;
  int ngood, ncells, i, k;

  if (bestscore > 2*suboptimal_score_end) {
    threshold_score = bestscore - suboptimal_score_end;
  } else {
    threshold_score = bestscore/2;
  }
  if (threshold_score <= indexsize) {
    threshold_score = indexsize + 1;
  }

  debug11(printf("Entered Linkmatrix_get_cells_both with querystart %d, queryend %d, threshold score %d\n",
		 querystart,queryend,threshold_score));

  ncells = 0;
  for (querypos = querystart; querypos <= queryend; querypos++) {
    ngood = 0;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if (links[querypos][hit].fwd_score >= threshold_score) {
	ngood++;
      }
#ifdef SEPARATE_FWD_REV
      if (links[querypos][hit].rev_score >= threshold_score) {
	ngood++;
      }
#endif
    }
    if (ngood > 0 && ngood <= 10) {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (links[querypos][hit].fwd_score >= threshold_score) {
	  rootposition = links[querypos][hit].fwd_rootposition;
	  /* tracei = links[querypos][hit].fwd_tracei; */
	  celllist = Cellpool_push(celllist,cellpool,rootposition,querypos,hit,/*fwdp*/true,
				   links[querypos][hit].fwd_score);
	  ncells++;
	}
#ifdef SEPARATE_FWD_REV
	if (links[querypos][hit].rev_score >= threshold_score) {
	  rootposition = links[querypos][hit].rev_rootposition;
	  /* tracei = links[querypos][hit].rev_tracei; */
	  celllist = Cellpool_push(celllist,cellpool,rootposition,querypos,hit,/*fwdp*/false,
				   links[querypos][hit].rev_score);
	  ncells++;
	}
#endif
      }
    }
  }

  if (ncells == 0) {
    *nunique = 0;
    return (Cell_T *) NULL;

  } else {
    /* Take best result for each tracei */
    cells = (Cell_T *) List_to_array(celllist,NULL);
    /* List_free(&celllist); -- no need with cellpool */

    if (favor_right_p == true) {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_right_cmp);
    } else {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_left_cmp);
    }

    sorted = (Cell_T *) CALLOC(ncells,sizeof(Cell_T));
    k = 0;

    last_rootposition = -1;
    for (i = 0; i < ncells; i++) {
      if (cells[i]->rootposition == last_rootposition) {
	/* Cell_free(&(cells[i])); -- no need with cellpool */
      } else {
	debug11(printf("rootposition %d, score %d, pos %d, hit %d\n",
		       cells[i]->rootposition,cells[i]->score,cells[i]->querypos,cells[i]->hit));
	sorted[k++] = cells[i];
	last_rootposition = cells[i]->rootposition;
      }
    }
    debug11(printf("\n"));
    FREE(cells);

    *nunique = k;
    qsort(sorted,*nunique,sizeof(Cell_T),Cell_score_cmp);

    return sorted;
  }
}



static int
binary_search (int lowi, int highi, Chrpos_T *mappings, Chrpos_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));
  if (mappings == NULL) {
    return -1;
  } else {
    while (lowi < highi) {
      middlei = (lowi+highi)/2;
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,mappings[lowi],middlei,mappings[middlei],
		     highi,mappings[highi],goal));
      if (goal < mappings[middlei]) {
	highi = middlei;
      } else if (goal > mappings[middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d\n",middlei));
	return middlei;
      }
    }

    debug10(printf("binary search returns %d\n",highi));
    return highi;
  }
}


/* Returns celllist */
/* For PMAP, indexsize is in aa. */
static Cell_T *
align_compute_scores_lookback (int *ncells, struct Link_T **links, Chrpos_T **mappings, int *npositions, int totalpositions,
			       bool oned_matrix_p, Chrpos_T *minactive, Chrpos_T *maxactive, Cellpool_T cellpool,
			       int querystart, int queryend, int querylength,

			       Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,

			       int indexsize, int sufflookback, int nsufflookback, int maxintronlen,
#ifdef DEBUG9
			       char *queryseq_ptr,
#endif
			       bool anchoredp, int anchor_querypos, Chrpos_T anchor_position,
			       bool localp, bool skip_repetitive_p, 
			       bool use_canonical_p, int non_canonical_penalty, bool debug_graphic_p, bool favor_right_p) {
  Cell_T *cells;
  Link_T currlink, prevlink;
  int querypos, indexsize_nt, hit, low_hit, high_hit;
  int nskipped, min_hits, specific_querypos, specific_low_hit, specific_high_hit, next_querypos;
  debug9(int fwd_tracei = 0);
  Intlist_T processed = NULL;
  int best_overall_score = 0;
  int grand_fwd_score, grand_fwd_querypos, grand_fwd_hit, best_fwd_hit, best_fwd_score;
#ifdef SEPARATE_FWD_REV
  int grand_rev_score, grand_rev_querypos, grand_rev_hit, best_rev_hit, best_rev_score;
  debug9(int rev_tracei = 0);
#endif
  int **active, *firstactive, *nactive;
  Chrpos_T position, prevposition;
#if 0
  int *lastGT, *lastAG;
#ifndef PMAP
  int *lastCT, *lastAC;
#endif
#endif
#ifdef DEBUG9
  char *oligo;
#endif
#ifdef DEBUG10
  Link_T termlink = NULL;
#endif

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

#ifdef DEBUG9
  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
#endif
  debug0(printf("Lookback: querystart = %d, queryend = %d, indexsize = %d\n",querystart,queryend,indexsize));

  if (oned_matrix_p == true) {
    active = intmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    active = intmatrix_2d_new(querylength,npositions);
  }

  firstactive = (int *) MALLOC(querylength * sizeof(int));
  nactive = (int *) MALLOC(querylength * sizeof(int));

  /* Initialize */
  for (querypos = 0; querypos < querystart; querypos++) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  }
  while (querypos <= queryend && npositions[querypos] <= 0) {
    debug9(printf("Skipping querypos %d which has no positions\n",querypos));
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
    querypos++;
  }

  if (anchoredp == true) {
    /* Guaranteed to find a hit */
    hit = binary_search(0,npositions[anchor_querypos],mappings[anchor_querypos],/*goal*/anchor_position);
    if (mappings[anchor_querypos] == NULL) {
      printf("mappings at anchor_querypos %d is NULL.  mappings = %p\n",anchor_querypos,mappings);
      abort();
    }

    currlink = &(links[anchor_querypos][hit]);	
#ifndef SEPARATE_FWD_REV
    currlink->fwd_pos = currlink->fwd_hit = -1;
    currlink->fwd_score = indexsize_nt;
    currlink->fwd_consecutive = EXON_DEFN;
    debug9(currlink->fwd_tracei = 0);
#else
    fprintf(stderr,"Not implemented yet\n");
    abort();
#endif

    debug6(printf("Setting firstactive for anchorpos %d to be %d\n",anchor_querypos,hit));
    firstactive[anchor_querypos] = hit;
    nactive[anchor_querypos] = 1;
    active[anchor_querypos][hit] = -1;

    debug6(printf("Pushing anchorpos %d as processed\n",anchor_querypos));
    processed = Intlist_push(processed,anchor_querypos);

  } else if (querypos <= queryend) {
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
#ifndef SEPARATE_FWD_REV
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      debug9(currlink->fwd_tracei = -1);
      /* currlink->fwd_rootnlinks = 1; */
#else
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      debug9(currlink->fwd_tracei = -1);
      /* currlink->fwd_rootnlinks = 1; */
      if (splicingp == true) {
	currlink->rev_pos = currlink->rev_hit = -1;
	currlink->rev_score = indexsize_nt;
	currlink->rev_consecutive = indexsize_nt;
	debug9(currlink->rev_tracei = -1);
	/* currlink->rev_rootnlinks = 1; */
      }
#endif
    }
    revise_active_lookback(active,firstactive,nactive,0,npositions[querypos],links,querypos,mappings);
  }

  grand_fwd_score = 0;
  grand_fwd_querypos = -1;
  grand_fwd_hit = -1;
#ifdef SEPARATE_FWD_REV
  if (splicingp == true) {
    grand_rev_score = 0;
    grand_rev_querypos = -1;
    grand_rev_hit = -1;
  }
#endif
  
  nskipped = 0;
  min_hits = 1000000;
  specific_querypos = -1;

  /* querypos += 1; -- this causes querypos at querystart to be ignored */
  while (querypos <= queryend) {
    best_fwd_score = 0;
    best_fwd_hit = -1;
#ifdef SEPARATE_FWD_REV
    best_rev_score = 0;
    best_rev_hit = -1;
#endif
    
#ifdef DEBUG9
    printf("Positions (forward order):");
    for (hit = 0; hit < npositions[querypos]; hit++) {
      printf(" %u",mappings[querypos][hit]);
    }
    printf("\n");
#endif

    hit = 0;
    while (hit < npositions[querypos] && mappings[querypos][hit] < minactive[querypos]) {
      hit++;
    }
    low_hit = hit;
    while (hit < npositions[querypos] && mappings[querypos][hit] <= maxactive[querypos]) {
      hit++;
    }
    high_hit = hit;
    debug9(printf("Querypos %d has hit %d..%d out of %d (minactive = %u, maxactive = %u)\n",
		  querypos,low_hit,high_hit-1,npositions[querypos],minactive[querypos],maxactive[querypos]));

    /* Can't use nactive yet, so use high_hit - low_hit */
    if (skip_repetitive_p && high_hit - low_hit >= MAX_NACTIVE && nskipped <= MAX_SKIPPED) { /* Previously turned off */
      debug6(printf("Too many active (%d - %d) at querypos %d.  Setting firstactive to be -1\n",high_hit,low_hit,querypos));
      firstactive[querypos] = -1;
      nactive[querypos] = 0;
      nskipped++;
      debug9(printf("  %d skipped because of %d hits\n",nskipped,high_hit - low_hit + 1));

      /* Store most specific querypos in section of skipped */
      if (high_hit - low_hit < min_hits) {
	min_hits = high_hit - low_hit;
	specific_querypos = querypos;
	specific_low_hit = low_hit;
	specific_high_hit = high_hit;
      }
      querypos++;

    } else {
      if (nskipped > MAX_SKIPPED) {
	debug9(printf("Too many skipped.  Going back to specific querypos %d\n",specific_querypos));
	next_querypos = querypos;
	querypos = specific_querypos;
	low_hit = specific_low_hit;
	high_hit = specific_high_hit;
      } else {
	next_querypos = querypos + 1;
      }

      for (hit = low_hit; hit < high_hit; hit++) {
	currlink = &(links[querypos][hit]);
	position = mappings[querypos][hit];
	  
	debug9(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	debug9(printf("Finding link at querypos %d,%d at %ux%d (%s).  prev_querypos was %d\n",
		      querypos,hit,position,active[querypos][hit],oligo,processed ? Intlist_head(processed) : -1));

	score_querypos_lookback(
#ifdef DEBUG9
				&fwd_tracei,
#endif
				currlink,querypos,querystart,queryend,position,
				links,mappings,active,firstactive,
				chroffset,chrhigh,plusp,
				indexsize,processed,sufflookback,nsufflookback,maxintronlen,
				anchoredp,localp,splicingp,skip_repetitive_p,use_canonical_p,
				non_canonical_penalty);

	if (currlink->fwd_score > best_fwd_score) {
	  best_fwd_score = currlink->fwd_score;
	  best_fwd_hit = hit;
	}
      }

      if (best_fwd_score > best_overall_score) {
	best_overall_score = best_fwd_score;
      }

      nskipped = 0;
      min_hits = 1000000;
      specific_querypos = -1;
      
#ifndef SEPARATE_FWD_REV
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d\n",
		    querypos,best_fwd_hit));
#else
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d and best_rev_hit %d\n",
		    querypos,best_fwd_hit,best_rev_hit));
#endif

      if (splicingp == true && best_fwd_hit >= 0 && links[querypos][best_fwd_hit].fwd_hit < 0 && 
	  grand_fwd_querypos >= 0 && querypos >= grand_fwd_querypos + indexsize_nt) {
	prevlink = &(links[grand_fwd_querypos][grand_fwd_hit]);
	if ((best_fwd_score = prevlink->fwd_score - (querypos - grand_fwd_querypos)) > 0) {
	  prevposition = mappings[grand_fwd_querypos][grand_fwd_hit];
	  debug12(printf("Considering prevposition %u to position %u as a grand fwd lookback\n",prevposition,position));
	  if (position > prevposition + maxintronlen) {
	    debug12(printf("  => Too long\n"));
	  } else {
	    for (hit = low_hit; hit < high_hit; hit++) {
	      currlink = &(links[querypos][hit]);
	      if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
		currlink->fwd_consecutive = indexsize_nt;
		/* currlink->fwd_rootnlinks = 1; */
		currlink->fwd_pos = grand_fwd_querypos;
		currlink->fwd_hit = grand_fwd_hit;
		currlink->fwd_score = best_fwd_score;
#ifdef DEBUG9
		currlink->fwd_tracei = ++fwd_tracei;
		currlink->fwd_intronnfwd = prevlink->fwd_intronnfwd;
		currlink->fwd_intronnrev = prevlink->fwd_intronnrev;
		currlink->fwd_intronnunk = prevlink->fwd_intronnunk + 1;
#endif
	      }
	    }
	    debug12(printf("At querypos %d, setting all fwd hits to point back to grand_fwd %d,%d with a score of %d\n",
			 querypos,grand_fwd_querypos,grand_fwd_hit,prevlink->fwd_score));
	  }
	}
      }

      /* Use >= to favor longer path in case of ties */
      if (best_fwd_hit >= 0 && best_fwd_score >= grand_fwd_score && 
	  links[querypos][best_fwd_hit].fwd_consecutive > EXON_DEFN) {
	grand_fwd_score = best_fwd_score;
	grand_fwd_querypos = querypos;
	grand_fwd_hit = best_fwd_hit;
	debug12(termlink = &(links[querypos][best_fwd_hit]));
	debug12(printf("At querypos %d, revising grand fwd to be hit %d with score of %d (pointing back to %d,%d)\n",
		       querypos,best_fwd_hit,best_fwd_score,termlink->fwd_pos,termlink->fwd_hit));
      }

#ifdef SEPARATE_FWD_REV
      if (best_rev_score > best_overall_score) {
	best_overall_score = best_rev_score;
      }

      if (splicingp == false || use_canonical_p == false) {
	/* rev scores should be the same as the fwd scores */
      } else {
	if (best_rev_hit >= 0 && links[querypos][best_rev_hit].rev_hit < 0 && 
	    grand_rev_querypos >= 0 && querypos >= grand_rev_querypos + indexsize_nt) {
	  prevlink = &(links[grand_rev_querypos][grand_rev_hit]);
	  if ((best_rev_score = prevlink->rev_score - (querypos - grand_rev_querypos)) > 0) {
	    prevposition = mappings[grand_rev_querypos][grand_rev_hit];
	    debug12(printf("Considering prevposition %u to position %u as a grand rev lookback\n",prevposition,position));
	    if (position > prevposition + maxintronlen) {
	      debug12(printf("  => Too long\n"));
	    } else {
	      for (hit = low_hit; hit < high_hit; hit++) {
		currlink = &(links[querypos][hit]);
		if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
		  currlink->rev_consecutive = indexsize_nt;
		  /* currlink->rev_rootnlinks = 1; */
		  currlink->rev_pos = grand_rev_querypos;
		  currlink->rev_hit = grand_rev_hit;
		  currlink->rev_score = best_rev_score;
#ifdef DEBUG9
		  currlink->rev_tracei = ++rev_tracei;
		  currlink->rev_intronnrev = prevlink->rev_intronnfwd;
		  currlink->rev_intronnrev = prevlink->rev_intronnrev;
		  currlink->rev_intronnunk = prevlink->rev_intronnunk + 1;
#endif
		}
	      }
	      debug12(printf("At querypos %d, setting all rev hits to point back to grand_rev %d,%d with a score of %d\n",
			     querypos,grand_rev_querypos,grand_rev_hit,prevlink->rev_score));
	    }
	  }
	}

	/* Use >= to favor longer path in case of ties */
	if (best_rev_hit >= 0 && best_rev_score >= grand_rev_score &&
	    links[querypos][best_rev_hit].rev_consecutive > EXON_DEFN) {
	  grand_rev_score = best_rev_score;
	  grand_rev_querypos = querypos;
	  grand_rev_hit = best_rev_hit;
	}
      }
#endif

      revise_active_lookback(active,firstactive,nactive,low_hit,high_hit,links,querypos,mappings);
      /* Need to push querypos, even if firstactive[querypos] == -1 */
      debug6(printf("Pushing querypos %d onto processed\n",querypos));
      processed = Intlist_push(processed,querypos);
      querypos = next_querypos;
    }
  }

  Intlist_free(&processed);

  /* These are the final active oligomers, after pruning by score */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,active,firstactive,indexsize,"active.mers");
  }

  FREE(nactive);
  FREE(firstactive);

  if (oned_matrix_p == true) {
    intmatrix_1d_free(&active);
  } else {
    intmatrix_2d_free(&active,querylength);
  }


  /* Grand winners */
  debug12(printf("Finding grand winners, using root position method\n"));
#ifdef SEPARATE_FWD_REV
  if (splicingp == false || use_canonical_p == false) {
    cells = Linkmatrix_get_cells_fwd(&(*ncells),links,querystart,queryend,npositions,
				     indexsize,best_overall_score,favor_right_p,cellpool);
  } else {
    cells = Linkmatrix_get_cells_both(&(*ncells),links,querystart,queryend,npositions,
				      indexsize,best_overall_score,favor_right_p,cellpool);
  }
#else
  cells = Linkmatrix_get_cells_fwd(&(*ncells),links,querystart,queryend,npositions,
				   indexsize,best_overall_score,favor_right_p,cellpool);
#endif

  debug9(FREE(oligo));

  return cells;
}


static char complCode[128] = COMPLEMENT_LC;

/* genomicstart == chroffset + chrpos */
/* arguments were genomicpos, genomicstart, genomiclength */

static char
get_genomic_nt (char *g_alt, Chrpos_T chrpos, Univcoord_T chroffset,
		Univcoord_T chrhigh, bool watsonp) {
  char c2, c2_alt;

  if (watsonp) {
    return Genome_get_char_blocks(&(*g_alt),chroffset + chrpos);

  } else {
    c2 = Genome_get_char_blocks(&c2_alt,chrhigh - chrpos);
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}


static List_T
traceback_one (int querypos, int hit, struct Link_T **links, Chrpos_T **mappings,
	       char *queryseq_ptr, char *queryuc_ptr, 
#ifdef PMAP
	       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
#ifdef DEBUG0
	       int indexsize,
#endif	       
	       Pairpool_T pairpool, bool fwdp) {
  List_T path = NULL;
  Chrpos_T position;
  int prev_querypos, prevhit;
  char c2;
#ifdef PMAP
  char c2_alt;
#endif

#ifdef DEBUG0
  char *oligo;
#endif


  while (querypos >= 0) {
    position = mappings[querypos][hit];

#ifdef PMAP
    /* Change querypos positions from protein to nucleotide */
    c2 = get_genomic_nt(&c2_alt,position+2,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+2,position+2,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position+1,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+1,position+1,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3,position,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#else
    /* genomic nucleotide same as queryseq */
    c2 = queryuc_ptr[querypos];
    path = Pairpool_push(path,pairpool,querypos,position,queryseq_ptr[querypos],MATCH_COMP,
			 c2,/*genomealt*/c2,/*dynprogindex*/0);
#endif


#ifdef DEBUG0
    debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
    debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
    if (fwdp == true) {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d",
		    querypos,hit,oligo,position,
		    links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive));
      debug9(printf(" (from #%d), intr = %d(+)/%d(-)/%d(?)",
		    links[querypos][hit].fwd_tracei,links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
      debug0(printf("\n"));

#ifdef SEPARATE_FWD_REV
    } else {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d",
		    querypos,hit,oligo,position,
		    links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive));
      debug9(printf(" (from #%d), intr = %d(+)/%d(-)/%d(?)",
		    links[querypos][hit].rev_tracei,links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		    links[querypos][hit].rev_intronnunk));
      debug0(printf("\n"));

#endif
    }
#endif
    debug0(FREE(oligo));

    /* prevposition = position; */
    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp == true) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifdef SEPARATE_FWD_REV
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
    debug3(printf("%d %d  %d %d  3\n",prev_querypos,prevhit,querypos,hit));
  }

  return path;
}


static List_T
traceback_one_snps (int querypos, int hit, struct Link_T **links, Chrpos_T **mappings,
		    char *queryseq_ptr, char *queryuc_ptr, 

		    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#ifdef DEBUG0
		    int indexsize,
#endif
		    Pairpool_T pairpool, bool fwdp) {
  List_T path = NULL;
  Chrpos_T position;
  int prev_querypos, prevhit;
  char c2, c2_alt;

#ifdef DEBUG0
  char *oligo;
#endif


  while (querypos >= 0) {
    position = mappings[querypos][hit];

#ifdef PMAP
    /* Change querypos positions from protein to nucleotide */
    c2 = get_genomic_nt(&c2_alt,position+2,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+2,position+2,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position+1,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+1,position+1,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3,position,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#else
    /* genomic nucleotide same as queryseq */
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos,position,queryseq_ptr[querypos],MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#endif


#ifdef DEBUG0
    debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
    debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
    if (fwdp == true) {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d",
		    querypos,hit,oligo,position,
		    links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive));
      debug9(printf(" (from #%d), intr = %d(+)/%d(-)/%d(?)",
		    links[querypos][hit].fwd_tracei,links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
      debug0(printf("\n"));
#ifdef SEPARATE_FWD_REV
    } else {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d",
		    querypos,hit,oligo,position,
		    links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive));
      debug9(printf(" (from #%d), intr = %d(+)/%d(-)/%d(?)",
		    links[querypos][hit].rev_tracei,links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		    links[querypos][hit].rev_intronnunk));
      debug0(printf("\n"));
#endif
    }
#endif
    debug0(FREE(oligo));

    /* prevposition = position; */
    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp == true) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifdef SEPARATE_FWD_REV
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
    debug3(printf("%d %d  %d %d  3\n",prev_querypos,prevhit,querypos,hit));
  }

  return path;
}


/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute_lookback (Chrpos_T **mappings, int *npositions, int totalpositions,
			bool oned_matrix_p, Chrpos_T *minactive, Chrpos_T *maxactive, Cellpool_T cellpool,
			char *queryseq_ptr, char *queryuc_ptr, int querylength, int querystart, int queryend,
			Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			int indexsize, int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool,
			bool anchoredp, int anchor_querypos, Chrpos_T anchor_position,
			bool localp, bool skip_repetitive_p, bool use_canonical_p, int non_canonical_penalty,
			bool favor_right_p, int max_nalignments, bool debug_graphic_p) {
  List_T all_paths = NULL;
  struct Link_T **links;

  Cell_T *cells, cell;
  int ncells, i;

  bool fwdp;
  int querypos, hit;
  int bestscore;


  if (oned_matrix_p == true) {
    links = Linkmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    links = Linkmatrix_2d_new(querylength,npositions);
  }

  /* These are all oligomers */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,/*active*/NULL,/*firstactive*/NULL,indexsize,"all.mers");
  }
  
  cells = align_compute_scores_lookback(&ncells,links,mappings,npositions,totalpositions,
					oned_matrix_p,minactive,maxactive,cellpool,
					querystart,queryend,querylength,
			       
					chroffset,chrhigh,plusp,

					indexsize,sufflookback,nsufflookback,maxintronlen,
#ifdef DEBUG9
					queryseq_ptr,
#endif
					anchoredp,anchor_querypos,anchor_position,localp,
					skip_repetitive_p,use_canonical_p,non_canonical_penalty,
					debug_graphic_p,favor_right_p);

#ifdef PMAP
  debug1(Linkmatrix_print_fwd(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#else
  debug1(Linkmatrix_print_both(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#endif

  if (ncells == 0) {
    all_paths = (List_T) NULL;

  } else {
    bestscore = cells[0]->score;

    debug11(printf("Looping on %d cells, allowing up to %d alignments, plus any with best score %d\n",
		   ncells,max_nalignments,bestscore));

    for (i = 0; i < ncells && (i < max_nalignments || cells[i]->score == bestscore)
	   && cells[i]->score > bestscore - FINAL_SCORE_TOLERANCE; i++) {
      cell = cells[i];
      querypos = cell->querypos;
      hit = cell->hit;
      fwdp = cell->fwdp;
      debug11(printf("Starting subpath %d for rootposition %d with score %d, querypos %d, hit %d, position %u\n",
		     i,cell->rootposition,cell->score,querypos,hit,mappings[querypos][hit]));


      if (debug_graphic_p == true) {
	best_path_dump_R(links,mappings,querypos,hit,fwdp,"best.path");
	printf("plot(all.mers,col=\"black\",pch=\".\",xlab=\"Query\",ylab=\"Genomic\")\n");
	printf("points(active.mers,col=\"red\",pch=\".\")\n");
	printf("points(best.path,col=\"green\",pch=\".\")\n");
	printf("lines(querypos,minactive,col=\"blue\")\n");
	printf("lines(querypos,maxactive,col=\"blue\")\n");
      }

      if (snps_p == true) {
	all_paths = List_push(all_paths,(void *) traceback_one_snps(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
								    chroffset,chrhigh,/*watsonp*/plusp,
#ifdef DEBUG0
								    indexsize,
#endif
								    pairpool,fwdp));
      } else {
	all_paths = List_push(all_paths,(void *) traceback_one(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
#ifdef PMAP
							       chroffset,chrhigh,/*watsonp*/plusp,
#endif
#ifdef DEBUG0
							       indexsize,
#endif
							       pairpool,fwdp));
      }

    }
    debug11(printf("\n"));

#if 0
    /* No need with cellpool */
    for (i = 0; i < ncells; i++) {
      cell = cells[i];
      Cell_free(&cell);
    }
#endif
    FREE(cells);
  }


  if (oned_matrix_p == true) {
    Linkmatrix_1d_free(&links);
  } else {
    Linkmatrix_2d_free(&links,querylength);
  }

#if 0
  for (p = all_paths; p != NULL; p = List_next(p)) {
    Pair_dump_list(List_head(p),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  return all_paths;
}



/* Returns celllist */
/* For PMAP, indexsize is in aa. */
static Cell_T *
align_compute_scores_lookforward (int *ncells, struct Link_T **links, Chrpos_T **mappings, int *npositions, int totalpositions,
				  bool oned_matrix_p, Chrpos_T *minactive, Chrpos_T *maxactive, Cellpool_T cellpool,
				  int querystart, int queryend, int querylength,
				  Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
				  int indexsize, int sufflookback, int nsufflookback, int maxintronlen,
#ifdef DEBUG9
				  char *queryseq_ptr,
#endif
				  bool anchoredp, int anchor_querypos, Chrpos_T anchor_position,
				  bool localp, bool skip_repetitive_p, 
				  bool use_canonical_p, int non_canonical_penalty,
				  bool debug_graphic_p, bool favor_right_p) {
  Cell_T *cells;
  Link_T currlink, prevlink;
  int querypos, indexsize_nt, hit, low_hit, high_hit;
  int nskipped, min_hits, specific_querypos, specific_low_hit, specific_high_hit, next_querypos;
  debug9(int fwd_tracei = 0);
  Intlist_T processed = NULL;
  int best_overall_score = 0;
  int grand_fwd_score, grand_fwd_querypos, grand_fwd_hit, best_fwd_hit, best_fwd_score;
#ifdef SEPARATE_FWD_REV
  int grand_rev_score, grand_rev_querypos, grand_rev_hit, best_rev_hit, best_rev_score;
  debug9(int rev_tracei = 0);
#endif
  int **active, *firstactive, *nactive;
  Chrpos_T position, prevposition;
#if 0
  int *lastGT, *lastAG;
#ifndef PMAP
  int *lastCT, *lastAC;
#endif
#endif
#ifdef DEBUG9
  char *oligo;
#endif
#ifdef DEBUG10
  Link_T termlink = NULL;
#endif

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

#ifdef DEBUG9
  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
#endif
  debug0(printf("Lookforward: querystart = %d, queryend = %d, indexsize = %d\n",querystart,queryend,indexsize));

  if (oned_matrix_p == true) {
    active = intmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    active = intmatrix_2d_new(querylength,npositions);
  }

  firstactive = (int *) MALLOC(querylength * sizeof(int));
  nactive = (int *) MALLOC(querylength * sizeof(int));

  /* Initialize */
  for (querypos = querylength - 1; querypos > queryend; querypos--) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  }
  while (querypos >= querystart && npositions[querypos] <= 0) {
    debug9(printf("Skipping querypos %d which has no positions\n",querypos));
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
    querypos--;
  }
  if (anchoredp == true) {
    /* Guaranteed to find a hit */
    hit = binary_search(0,npositions[anchor_querypos],mappings[anchor_querypos],/*goal*/anchor_position);
    if (mappings[anchor_querypos] == NULL) {
      printf("mappings at anchor_querypos %d is NULL.  mappings = %p\n",anchor_querypos,mappings);
      abort();
    }

    currlink = &(links[anchor_querypos][hit]);	
#ifndef SEPARATE_FWD_REV
    currlink->fwd_pos = currlink->fwd_hit = -1;
    currlink->fwd_score = indexsize_nt;
    currlink->fwd_consecutive = EXON_DEFN;
    debug9(currlink->fwd_tracei = 0);
#else
    fprintf(stderr,"Not implemented yet\n");
    abort();
#endif

    debug6(printf("Setting firstactive for anchorpos %d to be %d\n",anchor_querypos,hit));
    firstactive[anchor_querypos] = hit;
    nactive[anchor_querypos] = 1;
    active[anchor_querypos][hit] = -1;

    debug6(printf("Pushing anchorpos %d as processed\n",anchor_querypos));
    processed = Intlist_push(processed,anchor_querypos);

  } else if (querypos >= querystart) {
    for (hit = npositions[querypos] - 1; hit >= 0; --hit) {
      currlink = &(links[querypos][hit]);
#ifndef SEPARATE_FWD_REV
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      debug9(currlink->fwd_tracei = -1);
      /* currlink->fwd_rootnlinks = 1; */
#else
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      debug9(currlink->fwd_tracei = -1);
      /* currlink->fwd_rootnlinks = 1; */
      if (splicingp == true) {
	currlink->rev_pos = currlink->rev_hit = -1;
	currlink->rev_score = indexsize_nt;
	currlink->rev_consecutive = indexsize_nt;
	currlink->rev_tracei = -1;
	/* currlink->rev_rootnlinks = 1; */
      }
#endif
    }
    revise_active_lookforward(active,firstactive,nactive,0,npositions[querypos],links,querypos,mappings);
  }


  grand_fwd_score = 0;
  grand_fwd_querypos = -1;
  grand_fwd_hit = -1;
#ifdef SEPARATE_FWD_REV
  if (splicingp == true) {
    grand_rev_score = 0;
    grand_rev_querypos = -1;
    grand_rev_hit = -1;
  }
#endif

  nskipped = 0;
  min_hits = 1000000;
  specific_querypos = -1;

  /* querypos -= 1; -- this causes querypos at queryend to be ignored */
  while (querypos >= querystart) {
    best_fwd_score = 0;
    best_fwd_hit = -1;
#ifdef SEPARATE_FWD_REV
    best_rev_score = 0;
    best_rev_hit = -1;
#endif
    
#ifdef DEBUG9
    printf("Positions (reverse order):");
    for (hit = npositions[querypos] - 1; hit >= 0; --hit) {
      printf(" %u",mappings[querypos][hit]);
    }
    printf("\n");
#endif

    hit = npositions[querypos] - 1;
    while (hit >= 0 && mappings[querypos][hit] > maxactive[querypos]) {
      --hit;
    }
    high_hit = hit + 1;
    while (hit >= 0 && mappings[querypos][hit] >= minactive[querypos]) {
      --hit;
    }
    low_hit = hit + 1;
    debug9(printf("Querypos %d has hit %d..%d out of %d (minactive = %u, maxactive = %u)\n",
		  querypos,high_hit-1,low_hit,npositions[querypos],minactive[querypos],maxactive[querypos]));

    /* Can't use nactive yet, so use high_hit - low_hit */
    if (skip_repetitive_p && high_hit - low_hit >= MAX_NACTIVE && nskipped <= MAX_SKIPPED) { /* Previously turned off */
      debug6(printf("Too many active (%d - %d) at querypos %d.  Setting firstactive to be -1\n",high_hit,low_hit,querypos));
      firstactive[querypos] = -1;
      nactive[querypos] = 0;
      nskipped++;
      debug9(printf("  %d skipped because of %d hits\n",nskipped,high_hit - low_hit + 1));

      /* Store most specific querypos in section of skipped */
      if (high_hit - low_hit < min_hits) {
	min_hits = high_hit - low_hit;
	specific_querypos = querypos;
	specific_low_hit = low_hit;
	specific_high_hit = high_hit;
      }
      querypos--;

    } else {
      if (nskipped > MAX_SKIPPED) {
	debug9(printf("Too many skipped.  Going back to specific querypos %d\n",specific_querypos));
	next_querypos = querypos;
	querypos = specific_querypos;
	low_hit = specific_low_hit;
	high_hit = specific_high_hit;
      } else {
	next_querypos = querypos - 1;
      }

      for (hit = high_hit - 1; hit >= low_hit; --hit) {
	currlink = &(links[querypos][hit]);
	position = mappings[querypos][hit];
	
	debug9(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	debug9(printf("Finding link at querypos %d,%d at %ux%d (%s).  prev_querypos was %d\n",
		      querypos,hit,position,active[querypos][hit],oligo,processed ? Intlist_head(processed) : -1));
	
	score_querypos_lookforward(
#ifdef DEBUG9
				   &fwd_tracei,
#endif
				   currlink,querypos,querystart,queryend,position,
				   links,mappings,active,firstactive,
				   chroffset,chrhigh,plusp,
				   indexsize,processed,sufflookback,nsufflookback,maxintronlen,
				   anchoredp,localp,splicingp,skip_repetitive_p,use_canonical_p,
				   non_canonical_penalty);

	if (currlink->fwd_score > best_fwd_score) {
	  best_fwd_score = currlink->fwd_score;
	  best_fwd_hit = hit;
	}
      }

      if (best_fwd_score > best_overall_score) {
	best_overall_score = best_fwd_score;
      }

      nskipped = 0;
      min_hits = 1000000;
      specific_querypos = -1;
      
#ifndef SEPARATE_FWD_REV
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d\n",
		   querypos,best_fwd_hit));
#else
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d and best_rev_hit %d\n",
		   querypos,best_fwd_hit,best_rev_hit));
#endif

      if (splicingp == true && best_fwd_hit >= 0 && links[querypos][best_fwd_hit].fwd_hit < 0 && 
	  grand_fwd_querypos <= querylength - indexsize_nt && querypos + indexsize_nt <= grand_fwd_querypos) {
	prevlink = &(links[grand_fwd_querypos][grand_fwd_hit]);
	if ((best_fwd_score = prevlink->fwd_score - (grand_fwd_querypos - querypos)) > 0) {
	  prevposition = mappings[grand_fwd_querypos][grand_fwd_hit];
	  debug12(printf("Considering prevposition %u to position %u as a grand fwd lookback\n",prevposition,position));
	  if (position > prevposition + maxintronlen) {
	    debug12(printf("  => Too long\n"));
	  } else {
	    for (hit = high_hit - 1; hit >= low_hit; --hit) {
	      currlink = &(links[querypos][hit]);
	      if ((position = mappings[querypos][hit]) + indexsize_nt <= prevposition) {
		currlink->fwd_consecutive = indexsize_nt;
		/* currlink->fwd_rootnlinks = 1; */
		currlink->fwd_pos = grand_fwd_querypos;
		currlink->fwd_hit = grand_fwd_hit;
		currlink->fwd_score = best_fwd_score;
#ifdef DEBUG9
		currlink->fwd_tracei = ++fwd_tracei;
		currlink->fwd_intronnfwd = prevlink->fwd_intronnfwd;
		currlink->fwd_intronnrev = prevlink->fwd_intronnrev;
		currlink->fwd_intronnunk = prevlink->fwd_intronnunk + 1;
#endif
	      }
	    }
	    debug12(printf("At querypos %d, setting all fwd hits to point back to grand_fwd %d,%d with a score of %d\n",
			   querypos,grand_fwd_querypos,grand_fwd_hit,prevlink->fwd_score));
	  }
	}
      }

      /* Use >= to favor longer path in case of ties */
      if (best_fwd_hit >= 0 && best_fwd_score >= grand_fwd_score && 
	  links[querypos][best_fwd_hit].fwd_consecutive > EXON_DEFN) {
	grand_fwd_score = best_fwd_score;
	grand_fwd_querypos = querypos;
	grand_fwd_hit = best_fwd_hit;
	debug12(termlink = &(links[querypos][best_fwd_hit]));
	debug12(printf("At querypos %d, revising grand fwd to be hit %d with score of %d (pointing back to %d,%d)\n",
		     querypos,best_fwd_hit,best_fwd_score,termlink->fwd_pos,termlink->fwd_hit));
      }

#ifdef SEPARATE_FWD_REV
      if (best_rev_score > best_overall_score) {
	best_overall_score = best_rev_score;
      }

      if (splicingp == false || use_canonical_p == false) {
	/* rev scores should be the same as the fwd scores */
      } else {
	if (best_rev_hit >= 0 && links[querypos][best_rev_hit].rev_hit < 0 && 
	    grand_rev_querypos <= querylength - indexsize_nt && querypos + indexsize_nt <= grand_rev_querypos) {
	  prevlink = &(links[grand_rev_querypos][grand_rev_hit]);
	  if ((best_rev_score = prevlink->rev_score - (grand_rev_querypos - querypos)) > 0) {
	    prevposition = mappings[grand_rev_querypos][grand_rev_hit];
	    debug12(printf("Considering prevposition %u to position %u as a grand rev lookback\n",prevposition,position));
	    if (position > prevposition + maxintronlen) {
	      debug12(printf("  => Too long\n"));
	    } else {
	      for (hit = high_hit - 1; hit >= low_hit; --hit) {
		currlink = &(links[querypos][hit]);
		if ((position = mappings[querypos][hit]) + indexsize_nt <= prevposition) {
		  currlink->rev_consecutive = indexsize_nt;
		  /* currlink->rev_rootnlinks = 1; */
		  currlink->rev_pos = grand_rev_querypos;
		  currlink->rev_hit = grand_rev_hit;
		  currlink->rev_score = best_rev_score;
#ifdef DEBUG9
		  currlink->rev_tracei = ++rev_tracei;
		  currlink->rev_intronnrev = prevlink->rev_intronnfwd;
		  currlink->rev_intronnrev = prevlink->rev_intronnrev;
		  currlink->rev_intronnunk = prevlink->rev_intronnunk + 1;
#endif
		}
	      }
	      debug12(printf("At querypos %d, setting all rev hits to point back to grand_rev %d,%d with a score of %d\n",
			     querypos,grand_rev_querypos,grand_rev_hit,prevlink->rev_score));
	    }
	  }
	}

	/* Use >= to favor longer path in case of ties */
	if (best_rev_hit >= 0 && best_rev_score >= grand_rev_score &&
	    links[querypos][best_rev_hit].rev_consecutive > EXON_DEFN) {
	  grand_rev_score = best_rev_score;
	  grand_rev_querypos = querypos;
	  grand_rev_hit = best_rev_hit;
	}
      }
#endif

      revise_active_lookforward(active,firstactive,nactive,low_hit,high_hit,links,querypos,mappings);
      /* Need to push querypos, even if firstactive[querypos] == -1 */
      debug6(printf("Pushing querypos %d onto processed\n",querypos));
      processed = Intlist_push(processed,querypos);
      querypos = next_querypos;
    }
  }

  Intlist_free(&processed);

  /* These are the final active oligomers, after pruning by score */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,active,firstactive,indexsize,"active.mers");
  }

  FREE(nactive);
  FREE(firstactive);

  if (oned_matrix_p == true) {
    intmatrix_1d_free(&active);
  } else {
    intmatrix_2d_free(&active,querylength);
  }


  /* Grand winners */
  debug12(printf("Finding grand winners, using root position method\n"));
#ifdef SEPARATE_FWD_REV
  if (splicingp == false || use_canonical_p == false) {
    cells = Linkmatrix_get_cells_fwd(&(*ncells),links,querystart,queryend,npositions,
				     indexsize,best_overall_score,favor_right_p,cellpool);
  } else {
    cells = Linkmatrix_get_cells_both(&(*ncells),links,querystart,queryend,npositions,
				      indexsize,best_overall_score,favor_right_p,cellpool);
  }
#else
  cells = Linkmatrix_get_cells_fwd(&(*ncells),links,querystart,queryend,npositions,
				   indexsize,best_overall_score,favor_right_p,cellpool);
#endif

  debug9(FREE(oligo));

  return cells;
}


/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute_lookforward (Chrpos_T **mappings, int *npositions, int totalpositions,
			   bool oned_matrix_p, Chrpos_T *minactive, Chrpos_T *maxactive, Cellpool_T cellpool,
			   char *queryseq_ptr, char *queryuc_ptr, int querylength, int querystart, int queryend,

			   Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			   int indexsize, int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool,
			   bool anchoredp, int anchor_querypos, Chrpos_T anchor_position,
			   bool localp, bool skip_repetitive_p, bool use_canonical_p, int non_canonical_penalty,
			   bool favor_right_p, int max_nalignments, bool debug_graphic_p) {
  List_T all_paths = NULL;
  struct Link_T **links;

  Cell_T *cells, cell;
  int ncells, i;

  bool fwdp;
  int querypos, hit;
  int bestscore;

  if (oned_matrix_p == true) {
    links = Linkmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    links = Linkmatrix_2d_new(querylength,npositions);
  }

  /* These are all oligomers */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,/*active*/NULL,/*firstactive*/NULL,indexsize,"all.mers");
  }
  
  cells = align_compute_scores_lookforward(&ncells,links,mappings,npositions,totalpositions,
					   oned_matrix_p,minactive,maxactive,cellpool,
					   querystart,queryend,querylength,
					   
					   chroffset,chrhigh,plusp,

					   indexsize,sufflookback,nsufflookback,maxintronlen,
#ifdef DEBUG9
					   queryseq_ptr,
#endif
					   anchoredp,anchor_querypos,anchor_position,localp,
					   skip_repetitive_p,use_canonical_p,non_canonical_penalty,
					   debug_graphic_p,favor_right_p);

#ifdef PMAP
  debug1(Linkmatrix_print_fwd(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#else
  debug1(Linkmatrix_print_both(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#endif

  if (ncells == 0) {
    all_paths = (List_T) NULL;

  } else {
    bestscore = cells[0]->score;

    debug11(printf("Looping on %d cells, allowing up to %d alignments, plus any with best score %d\n",
		   ncells,max_nalignments,bestscore));

    for (i = 0; i < ncells && (i < max_nalignments || cells[i]->score == bestscore)
	   && cells[i]->score > bestscore - FINAL_SCORE_TOLERANCE; i++) {
      cell = cells[i];
      querypos = cell->querypos;
      hit = cell->hit;
      fwdp = cell->fwdp;
      debug11(printf("Starting subpath %d at rootposition %d with score %d, querypos %d, hit %d, position %u\n",
		     i,cell->rootposition,cell->score,querypos,hit,mappings[querypos][hit]));


      if (debug_graphic_p == true) {
	best_path_dump_R(links,mappings,querypos,hit,fwdp,"best.path");
	printf("plot(all.mers,col=\"black\",pch=\".\",xlab=\"Query\",ylab=\"Genomic\")\n");
	printf("points(active.mers,col=\"red\",pch=\".\")\n");
	printf("points(best.path,col=\"green\",pch=\".\")\n");
	printf("lines(querypos,minactive,col=\"blue\")\n");
	printf("lines(querypos,maxactive,col=\"blue\")\n");
      }

      if (snps_p == true) {
	all_paths = List_push(all_paths,(void *) traceback_one_snps(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
								    chroffset,chrhigh,/*watsonp*/plusp,
#ifdef DEBUG0
								    indexsize,
#endif
								    pairpool,fwdp));
      } else {
	all_paths = List_push(all_paths,(void *) traceback_one(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
#ifdef PMAP
							       chroffset,chrhigh,/*watsonp*/plusp,
#endif
#ifdef DEBUG0
							       indexsize,
#endif
							       pairpool,fwdp));
      }

    }
    debug11(printf("\n"));

#if 0
    /* No need with cellpool */
    for (i = 0; i < ncells; i++) {
      cell = cells[i];
      Cell_free(&cell);
    }
#endif
    FREE(cells);
  }


  if (oned_matrix_p == true) {
    Linkmatrix_1d_free(&links);
  } else {
    Linkmatrix_2d_free(&links,querylength);
  }

#if 0
  for (p = all_paths; p != NULL; p = List_next(p)) {
    Pair_dump_list(List_head(p),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  return all_paths;
}



/* Modified from stage3.c */
static void
get_splicesite_probs (double *left_prob, double *right_prob,
		      Chrpos_T left_genomepos, Chrpos_T right_genomepos,
		      Univcoord_T chroffset, Univcoord_T chrhigh, int genestrand, bool watsonp) {
  Univcoord_T splicesitepos;
  
  if (watsonp == true) {
    splicesitepos = chroffset + left_genomepos;
    if (genestrand > 0) {
      *left_prob = Maxent_hr_donor_prob(splicesitepos /*?*/+ 1,chroffset);
      debug5(printf("1. donor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*left_prob));

    } else {
      *left_prob = Maxent_hr_antiacceptor_prob(splicesitepos /**/+ 1,chroffset);
      debug5(printf("2. antiacceptor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*left_prob));

    }
  } else {
    splicesitepos = chrhigh - left_genomepos + 1;
    if (genestrand > 0) {
      *left_prob = Maxent_hr_acceptor_prob(splicesitepos /*?*/- 1,chroffset);
      debug5(printf("4. acceptor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*left_prob));
    } else {
      *left_prob = Maxent_hr_antidonor_prob(splicesitepos /**/- 1,chroffset);
      debug5(printf("3. antidonor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*left_prob));
    }
  }

  if (watsonp == true) {
    splicesitepos = chroffset + right_genomepos + 1;
    if (genestrand > 0) {
      *right_prob = Maxent_hr_acceptor_prob(splicesitepos /*?*/- 1,chroffset);
      debug5(printf("5. acceptor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*right_prob));
    } else {
      *right_prob = Maxent_hr_antidonor_prob(splicesitepos /**/- 1,chroffset);
      debug5(printf("6. antidonor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*right_prob));

    }
  } else {
    splicesitepos = chrhigh - right_genomepos;
    if (genestrand > 0) {
      *right_prob = Maxent_hr_donor_prob(splicesitepos /*?*/+ 1,chroffset);
      debug5(printf("8. donor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*right_prob));
    } else {
      *right_prob = Maxent_hr_antiacceptor_prob(splicesitepos /**/+ 1,chroffset);
      debug5(printf("7. antiacceptor splicesitepos is %u (%u), prob %f\n",
		    splicesitepos,splicesitepos-chroffset,*right_prob));
    }
  }

  return;
}



/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides (List_T path,
#ifndef PMAP
			char *queryseq_ptr, char *queryuc_ptr, 
#endif
			Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			int query_offset, Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL, pairptr;
  Pair_T pair;
  int querypos, lastquerypos, queryjump, genomejump, fill, default_fill;
  int genomepos, lastgenomepos;
  char c, c_alt;
  double left_prob, right_prob;

  debug5(printf("Beginning convert_to_nucleotides with %d pairs\n",List_length(path)));

  /* pairptr = path; */
  /* path = Pairpool_pop(path,&pair); */
  pair = (Pair_T) path->first;
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  default_fill = indexsize_nt - 3;
#else
  default_fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + default_fill;
  lastgenomepos = genomepos + default_fill;
  while (lastquerypos > querypos) {
    debug5(printf("lastquerypos %d, lastgenomepos %d\n",
		  lastquerypos,lastgenomepos));

#ifdef PMAP
    c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			  /*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
    if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      abort();
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos));
    }
#else
    if (mode == STANDARD) {
      c = queryuc_ptr[lastquerypos];
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MATCH_COMP,c,/*genomealt*/c,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      if (queryuc_ptr[lastquerypos] == c) {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      } else {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      }
    }
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of observed first pair in oligomer */
  if (mode == STANDARD) {
    pair->querypos += query_offset; /* Revise coordinates */
    /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_transfer_one(pairs,&path);
#endif
  } else {
    c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
    if (pair->cdna == c) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    } else {
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			    pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
      debug5(printf("Pushing %c : %c at %d,%d (first pair)\n",pair->cdna,c,
		    pair->querypos+query_offset,pair->genomepos));
    }
  }

  lastquerypos = querypos;
  lastgenomepos = genomepos;

  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - 1 - querypos;
    genomejump = lastgenomepos - 1 - genomepos;

    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else {
      debug5(printf("At querypos %d, saw queryjump of %d and genomejump of %d\n",querypos,queryjump,genomejump));

      if (querypos + default_fill >= lastquerypos || genomepos + default_fill >= lastgenomepos) {
	if (lastquerypos - querypos < (int) (lastgenomepos - genomepos)) {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
#endif
	  fill = lastquerypos - querypos - 1;
	} else {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
#endif
	  fill = lastgenomepos - genomepos - 1;
	}
      } else {
	fill = default_fill;
      }

      /* Recompute queryjump and genomejump */
      queryjump -= fill;
      genomejump -= fill;
      debug5(printf("  Revised queryjump of %d and genomejump of %d\n",queryjump,genomejump));
      if (genomejump > 0 || queryjump > 0) {
	debug5(printf("  Pushing gapholder\n"));
	pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
					/*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#if 0
	/* Need to run on both genestrands, and save both results in pair */
	if (queryjump == 0) {
	  get_splicesite_probs(&left_prob,&right_prob,genomepos+fill,lastgenomepos,
			       chroffset,chrhigh,/*genestrand*/+1,watsonp);
	}
#endif
      }

      /* Fill rest of oligomer */
      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      debug5(printf("  Fill from querypos %d down to %d\n",lastquerypos,querypos));
      while (lastquerypos > querypos) {
#ifdef PMAP
	c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			      /*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos+query_offset,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
	if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,
				queryuc_ptr[lastquerypos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  abort();
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MISMATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	}
#else
	if (mode == STANDARD) {
	  /* assert(queryuc_ptr[lastquerypos] == get_genomic_nt(&c_alt,lastgenomepos,genomicstart,genomiclength,watsonp)); */
	  c = queryuc_ptr[lastquerypos];
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,c,/*genomealt*/c,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  if (queryuc_ptr[lastquerypos] == c) {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  } else {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  }
	}
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    if (mode == STANDARD) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    } else {
      c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
      if (pair->cdna == c) {
	pair->querypos += query_offset; /* Revise coordinates */
	/*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif
      } else {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			      pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d (observed)\n",pair->cdna,c,
		      pair->querypos+query_offset,pair->genomepos));
      }
    }

    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }

  debug5(Pair_dump_list(pairs,true));
  /* pairs is in ascending querypos order */
  return pairs;		      /* Used to return List_reverse(pairs) */
}


/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides_snps (List_T path,
#ifndef PMAP
			     char *queryseq_ptr, char *queryuc_ptr, 
#endif
			     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			     int query_offset, Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL, pairptr;
  Pair_T pair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, fill, default_fill;
  char c, c_alt;

  debug5(printf("Beginning convert_to_nucleotides with %d pairs\n",List_length(path)));

  /* pairptr = path; */
  /* path = Pairpool_pop(path,&pair); */
  pair = (Pair_T) path->first;
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  default_fill = indexsize_nt - 3;
#else
  default_fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + default_fill;
  lastgenomepos = genomepos + default_fill;
  while (lastquerypos > querypos) {
    debug5(printf("lastquerypos %d, lastgenomepos %d\n",
		  lastquerypos,lastgenomepos));

#ifdef PMAP
    c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			  /*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
    if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      abort();
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos));
    }
#else
    if (mode == STANDARD) {
      /* assert(queryuc_ptr[lastquerypos] == get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp)); */
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      if (queryuc_ptr[lastquerypos] == c) {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      } else {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      }
    }
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of observed first pair in oligomer */
  if (mode == STANDARD) {
    pair->querypos += query_offset; /* Revise coordinates */
    /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_transfer_one(pairs,&path);
#endif
  } else {
    c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
    if (pair->cdna == c) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    } else {
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			    pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
      debug5(printf("Pushing %c : %c at %d,%d (first pair)\n",pair->cdna,c,
		    pair->querypos+query_offset,pair->genomepos));
    }
  }

  lastquerypos = querypos;
  lastgenomepos = genomepos;

  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - 1 - querypos;
    genomejump = lastgenomepos - 1 - genomepos;

    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else {
      debug5(printf("At querypos %d, saw queryjump of %d and genomejump of %d\n",querypos,queryjump,genomejump));

      if (querypos + default_fill >= lastquerypos || genomepos + default_fill >= lastgenomepos) {
	if (lastquerypos - querypos < lastgenomepos - genomepos) {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
#endif
	  fill = lastquerypos - querypos - 1;
	} else {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
#endif
	  fill = lastgenomepos - genomepos - 1;
	}
      } else {
	fill = default_fill;
      }

      /* Recompute queryjump and genomejump */
      queryjump -= fill;
      genomejump -= fill;
      debug5(printf("  Revised queryjump of %d and genomejump of %d\n",queryjump,genomejump));
      if (genomejump > 0 || queryjump > 0) {
	debug5(printf("  Pushing gapholder\n"));
	pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
					/*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      }

      /* Fill rest of oligomer */
      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      debug5(printf("  Fill from querypos %d down to %d\n",lastquerypos,querypos));
      while (lastquerypos > querypos) {
#ifdef PMAP
	c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			      /*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos+query_offset,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
	if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,
				queryuc_ptr[lastquerypos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  abort();
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MISMATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	}
#else
	if (mode == STANDARD) {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  if (queryuc_ptr[lastquerypos] == c) {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  } else {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  }
	}
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    if (mode == STANDARD) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    } else {
      c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
      if (pair->cdna == c) {
	pair->querypos += query_offset; /* Revise coordinates */
	/*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif
      } else {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			      pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d (observed)\n",pair->cdna,c,
		      pair->querypos+query_offset,pair->genomepos));
      }
    }

    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }

  debug5(Pair_dump_list(pairs,true));
  /* pairs is in ascending querypos order */
  return pairs;		      /* Used to return List_reverse(pairs) */
}



/* Returns ncovered */
int
Stage2_scan (int *stage2_source, char *queryuc_ptr, int querylength,
#ifdef PMAP
	     char *genomicuc_ptr,
#endif
	     Chrpos_T chrstart, Chrpos_T chrend,
	     Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
	     int genestrand, Oligoindex_array_T oligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp) {
  int ncovered;
  int source;
  int indexsize;
  Oligoindex_T oligoindex;
  Chrpos_T **mappings;
  bool *coveredp, oned_matrix_p;
  int *npositions, totalpositions;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
#ifndef USE_DIAGPOOL
  List_p;
  Diag_T diag;
#endif
#ifdef EXTRACT_GENOMICSEG
  Count_T *counts;
#endif

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (Chrpos_T **) CALLOC(querylength,sizeof(Chrpos_T *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
  Diagpool_reset(diagpool);
  diagonals = (List_T) NULL;
  while (source < Oligoindex_array_length(oligoindices) && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */
#ifdef PMAP
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#else

#ifdef EXTRACT_GENOMICSEG
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
    counts = Oligoindex_counts_copy(oligoindex);
#endif

    if (plusp == true) {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend,/*plusp*/true,
			  queryuc_ptr,querylength,/*chrpos*/chrstart,genestrand);
    } else {
      /* Need to add 1 to mappingend to cover same range as plusp */
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend+1,/*plusp*/false,
			  queryuc_ptr,querylength,/*chrpos*/(chrhigh-chroffset)-chrend,genestrand);
    }

#ifdef EXTRACT_GENOMICSEG
    assert(Oligoindex_counts_equal(oligoindex,counts));
    /* Oligoindex_counts_dump(oligoindex,counts); */
    FREE(counts);
#endif

#endif

    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindices,oligoindex,queryuc_ptr,
					querylength,chrstart,chrend,chroffset,chrhigh,plusp,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    debug(printf("Stage2_scan: source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage));

    source++;
  }
  *stage2_source = source;

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif

  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

#if 1
  for (source = 0; source < Oligoindex_array_length(oligoindices); source++) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    Oligoindex_untally(oligoindex,queryuc_ptr,querylength);
  }
#endif

  return ncovered;
}


List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	
#ifdef PMAP
		char *genomicuc_ptr,
#endif
		Chrpos_T chrstart, Chrpos_T chrend,
		Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp, int genestrand,
		Oligoindex_array_T oligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		int sufflookback, int nsufflookback, int maxintronlen, bool localp, bool skip_repetitive_p,
		bool favor_right_p, int max_nalignments, bool debug_graphic_p, bool diagnosticp,
		Stopwatch_T stopwatch, bool diag_debug) {
  List_T all_stage2results = NULL, all_paths, all_ends, all_starts, end_paths, start_paths, path, pairs, p, q;
  List_T middle;
  Pair_T firstpair, lastpair;
  int diag_querystart, diag_queryend;
  int indexsize, indexsize_nt;
  Oligoindex_T oligoindex;
  Chrpos_T **mappings;
  bool *coveredp, oned_matrix_p;
  int source;
  int *npositions, totalpositions;
  Chrpos_T *minactive, *maxactive;
  int ncovered;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
  int anchor_querypos, querystart, queryend;
  Chrpos_T anchor_position;

#ifndef USE_DIAGPOOL
  List_T p;
  Diag_T diag;
#endif
#ifdef DEBUG
  int nunique;
#endif
#ifdef DEBUG0
  int i;
#endif

#ifdef EXTRACT_GENOMICSEG
  Count_T *counts;
#endif

  debug(printf("Entered Stage2_compute with chrstart %u and chrend %u\n",chrstart,chrend));

  Stopwatch_start(stopwatch);

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (Chrpos_T **) CALLOC(querylength,sizeof(Chrpos_T *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
#ifdef USE_DIAGPOOL
  Diagpool_reset(diagpool);
#endif
  Cellpool_reset(cellpool);
  diagonals = (List_T) NULL;
  while (source < Oligoindex_array_length(oligoindices) && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */

#ifdef PMAP
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#else

#if 0
    /* Previously used this for user_genomicseg, but now creating genome_blocks on the fly */
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#endif

#ifdef EXTRACT_GENOMICSEG
    /* printf("indexsize = %d\n",indexsize); */
    /* printf("Query:  %.*s\n",querylength,queryuc_ptr); */
    /* printf("Genome: %s\n",genomicuc_ptr); */
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/mappingend-mappingstart,
		     queryuc_ptr,querylength,sequencepos);
    counts = Oligoindex_counts_copy(oligoindex);

    /* printf("plusp %d\n",plusp); */
    /* printf("genomicstart %u, genomicend %u, genomiclength %d\n",genomicstart,genomicend,genomiclength); */
    /* printf("mappingstart %u, mappingend %u\n",mappingstart,mappingend); */
#endif

    if (plusp == true) {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend,/*plusp*/true,
			  queryuc_ptr,querylength,/*chrpos*/chrstart,genestrand);
    } else {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend+1,/*plusp*/false,
			  queryuc_ptr,querylength,/*chrpos*/(chrhigh-chroffset)-chrend,genestrand);
    }

#ifdef EXTRACT_GENOMICSEG
    assert(Oligoindex_counts_equal(oligoindex,counts));
    /* Oligoindex_counts_dump(oligoindex,counts); */

    FREE(counts);
#endif

#endif

    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindices,oligoindex,queryuc_ptr,
					querylength,chrstart,chrend,chroffset,chrhigh,plusp,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    debug(printf("Stage2_compute: source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage));

    source++;
  }
  *stage2_source = source;
  *stage2_indexsize = indexsize;
#ifdef PMAP
  indexsize_nt = 3*indexsize;
#else
  indexsize_nt = indexsize;
#endif

  /* diag_runtime = */ Stopwatch_stop(stopwatch);

  minactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));
  maxactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));

  Stopwatch_start(stopwatch);

  if (diag_debug == true) {
    /* Do nothing */
    middle = (List_T) NULL;

  } else if (totalpositions == 0) {
    debug(printf("Quitting because totalpositions is zero\n"));
    middle = (List_T) NULL;

  } else if (querylength > 150 && pct_coverage < proceed_pctcoverage && ncovered < SUFF_NCOVERED) {
    /* Filter only on long queries */
    debug(printf("Quitting because querylength %d > 150, and pct_coverage is only %f < %f, and ncovered is only %d < %d, maxnconsecutive = %d\n",
		 querylength,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED,maxnconsecutive));
    middle = (List_T) NULL;

  } else {
    debug(printf("Proceeding because maxnconsecutive is %d and pct_coverage is %f > %f or ncovered = %d > %d\n",
		 maxnconsecutive,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED));

    debug(printf("Performing diag on genomiclength %u\n",chrend-chrstart));
    Diag_compute_bounds(&diag_querystart,&diag_queryend,minactive,maxactive,diagonals,querylength,
			debug_graphic_p,chrstart,chrend,chroffset,chrhigh,plusp);
    
    debug(
	  nunique = Diag_compute_bounds(&diag_querystart,&diag_queryend,minactive,maxactive,diagonals,querylength,
					debug_graphic_p,chrstart,chrend,chroffset,chrhigh,plusp);
	  fprintf(stderr,"%d diagonals (%d not dominated), maxnconsecutive = %d\n",
		  List_length(diagonals),nunique,maxnconsecutive);
	  );

    if (debug_graphic_p == true) {
      active_bounds_dump_R(minactive,maxactive,querylength);
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    all_paths = align_compute_lookback(mappings,npositions,totalpositions,
				       oned_matrix_p,minactive,maxactive,cellpool,
				       queryseq_ptr,queryuc_ptr,querylength,
				       /*querystart*/diag_querystart,/*queryend*/diag_queryend,
				       chroffset,chrhigh,plusp,
				       indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
				       /*anchoredp*/false,/*anchor_querypos*/0,/*anchor_position*/0,
				       localp,skip_repetitive_p,use_canonical_middle_p,NON_CANONICAL_PENALTY_MIDDLE,
				       favor_right_p,max_nalignments,debug_graphic_p);
    for (p = all_paths; p != NULL; p = List_next(p)) {
      path = (List_T) p->first;
      firstpair = path->first;
      pairs = List_reverse(path);
      lastpair = pairs->first;

      if (snps_p == true) {
	middle = convert_to_nucleotides_snps(pairs,
#ifndef PMAP
					     queryseq_ptr,queryuc_ptr,
#endif
					     chroffset,chrhigh,/*watsonp*/plusp,
					     query_offset,pairpool,indexsize_nt);
      } else {
	middle = convert_to_nucleotides(pairs,
#ifndef PMAP
					queryseq_ptr,queryuc_ptr,
#endif
					chroffset,chrhigh,/*watsonp*/plusp,
					query_offset,pairpool,indexsize_nt);
      }

      debug0(printf("For end, anchor querypos %d, anchor_position %u\n",lastpair->querypos,lastpair->genomepos));
      debug0(printf("For start, anchor querypos %d, anchor_position %u\n",firstpair->querypos,firstpair->genomepos));

      anchor_querypos = lastpair->querypos;
      anchor_position = lastpair->genomepos;
      querystart = anchor_querypos + 1;
      queryend = querylength - 1;

      all_ends = (List_T) NULL;
      end_paths = align_compute_lookback(mappings,npositions,totalpositions,
					 oned_matrix_p,minactive,maxactive,cellpool,
					 queryseq_ptr,queryuc_ptr,querylength,querystart,queryend,
					 chroffset,chrhigh,plusp,
					 indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
					 /*anchoredp*/true,anchor_querypos,anchor_position,
					 localp,skip_repetitive_p,use_canonical_ends_p,NON_CANONICAL_PENALTY_ENDS,
					 favor_right_p,max_nalignments,debug_graphic_p);

      /* fprintf(stderr,"%d ends\n",List_length(end_paths)); */
      if (List_length(end_paths) == 1) {
	pairs = (List_T) List_head(end_paths);
	path = List_reverse(pairs);
	if (snps_p == true) {
	  pairs = convert_to_nucleotides_snps(path,
#ifndef PMAP
					      queryseq_ptr,queryuc_ptr,
#endif
					      chroffset,chrhigh,/*watsonp*/plusp,
					      query_offset,pairpool,indexsize_nt);
	} else {
	  pairs = convert_to_nucleotides(path,
#ifndef PMAP
					 queryseq_ptr,queryuc_ptr,
#endif
					 chroffset,chrhigh,/*watsonp*/plusp,
					 query_offset,pairpool,indexsize_nt);
	}
	middle = List_reverse(Pairpool_join_end3(List_reverse(middle),pairs,pairpool,/*copy_end_p*/false));
	debug0(printf("ATTACHING SINGLE END TO MIDDLE\n"));
	debug0(Pair_dump_list(middle,true));

      } else {
	debug0(i = 0);
	for (q = end_paths; q != NULL; q = List_next(q)) {
	  pairs = (List_T) List_head(q);
	  path = List_reverse(pairs);
	  if (snps_p == true) {
	    pairs = convert_to_nucleotides_snps(path,
#ifndef PMAP
						queryseq_ptr,queryuc_ptr,
#endif
						chroffset,chrhigh,/*watsonp*/plusp,
						query_offset,pairpool,indexsize_nt);
	  } else {
	    pairs = convert_to_nucleotides(path,
#ifndef PMAP
					   queryseq_ptr,queryuc_ptr,
#endif
					   chroffset,chrhigh,/*watsonp*/plusp,
					   query_offset,pairpool,indexsize_nt);
	  }
	  debug0(printf("END %d/%d\n",i++,List_length(end_paths)));
	  debug0(Pair_dump_list(pairs,true));
	  all_ends = List_push(all_ends,(void *) pairs);
	}
      }
      List_free(&end_paths);

    
      anchor_querypos = firstpair->querypos;
      anchor_position = firstpair->genomepos;
      querystart = 0;
      queryend = anchor_querypos - 1;

      all_starts = (List_T) NULL;
      start_paths = align_compute_lookforward(mappings,npositions,totalpositions,
					      oned_matrix_p,minactive,maxactive,cellpool,
					      queryseq_ptr,queryuc_ptr,querylength,querystart,queryend,
					      chroffset,chrhigh,plusp,
					      indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
					      /*anchoredp*/true,anchor_querypos,anchor_position,
					      localp,skip_repetitive_p,use_canonical_ends_p,NON_CANONICAL_PENALTY_ENDS,
					      favor_right_p,max_nalignments,debug_graphic_p);

      /* fprintf(stderr,"%d starts\n",List_length(start_paths)); */
      if (List_length(start_paths) == 1) {
	path = (List_T) List_head(start_paths);
	if (snps_p == true) {
	  pairs = convert_to_nucleotides_snps(path,
#ifndef PMAP
					      queryseq_ptr,queryuc_ptr,
#endif
					      chroffset,chrhigh,/*watsonp*/plusp,
					      query_offset,pairpool,indexsize_nt);
	} else {
	  pairs = convert_to_nucleotides(path,
#ifndef PMAP
					 queryseq_ptr,queryuc_ptr,
#endif
					 chroffset,chrhigh,/*watsonp*/plusp,
					 query_offset,pairpool,indexsize_nt);
	}
	path = List_reverse(pairs);
	middle = Pairpool_join_end5(middle,path,pairpool,/*copy_end_p*/false);
	debug0(printf("ATTACHING SINGLE START TO MIDDLE\n"));
	debug0(Pair_dump_list(middle,true));
      } else {
	debug0(i = 0);
	for (q = start_paths; q != NULL; q = List_next(q)) {
	  path = (List_T) List_head(q);
	  if (snps_p == true) {
	    pairs = convert_to_nucleotides_snps(path,
#ifndef PMAP
						queryseq_ptr,queryuc_ptr,
#endif
						chroffset,chrhigh,/*watsonp*/plusp,
						query_offset,pairpool,indexsize_nt);
	  } else {
	    pairs = convert_to_nucleotides(path,
#ifndef PMAP
					   queryseq_ptr,queryuc_ptr,
#endif
					   chroffset,chrhigh,/*watsonp*/plusp,
					   query_offset,pairpool,indexsize_nt);
	  }
	  path = List_reverse(pairs);
	  debug0(printf("START %d/%d\n",i++,List_length(start_paths)));
	  debug0(Pair_dump_list(path,true));
	  all_starts = List_push(all_starts,(void *) path);
	}
      }
      List_free(&start_paths);

      all_stage2results = List_push(all_stage2results,(void *) Stage2_new(middle,all_starts,all_ends));
    }

    List_free(&all_paths);
  }


  FREE(maxactive);
  FREE(minactive);
  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

#if 1
  for (source = 0; source < Oligoindex_array_length(oligoindices); source++) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    Oligoindex_untally(oligoindex,queryuc_ptr,querylength);
  }
#endif

  Stopwatch_stop(stopwatch);

  if (diag_debug == true) {
    return diagonals;
  } else {

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif
  }

  return all_stage2results;
}



List_T
Stage2_compute_one (int *stage2_source, int *stage2_indexsize,
		    char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	
#ifdef PMAP
		    char *genomicuc_ptr,
#endif
		    Chrpos_T chrstart, Chrpos_T chrend,
		    Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp, int genestrand,

		    Oligoindex_array_T oligoindices, double proceed_pctcoverage,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    int sufflookback, int nsufflookback,  int maxintronlen, bool localp,
		    bool skip_repetitive_p, bool use_shifted_canonical_p,
		    bool favor_right_p, bool debug_graphic_p, bool diagnosticp) {
  List_T pairs, all_paths;
  List_T middle, path;
  int indexsize, indexsize_nt;
  Oligoindex_T oligoindex;
  Chrpos_T **mappings;
  bool *coveredp, oned_matrix_p;
  int source;
  int *npositions, totalpositions;
  Chrpos_T *minactive, *maxactive;
  int ncovered;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;


  debug(printf("Entered Stage2_compute_one with chrstart %u and chrend %u\n",chrstart,chrend));

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (Chrpos_T **) CALLOC(querylength,sizeof(Chrpos_T *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
#ifdef USE_DIAGPOOL
  Diagpool_reset(diagpool);
#endif
  Cellpool_reset(cellpool);
  diagonals = (List_T) NULL;
  while (source < Oligoindex_array_length(oligoindices) && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */

#ifdef PMAP
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#else

    if (plusp == true) {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend,/*plusp*/true,
			  queryuc_ptr,querylength,/*chrpos*/chrstart,genestrand);
    } else {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend+1,/*plusp*/false,
			  queryuc_ptr,querylength,/*chrpos*/(chrhigh-chroffset)-chrend,genestrand);
    }

#endif

    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindices,oligoindex,queryuc_ptr,
					querylength,chrstart,chrend,chroffset,chrhigh,plusp,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    debug(printf("Stage2_compute: source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage));

    source++;
  }
  *stage2_source = source;
  *stage2_indexsize = indexsize;
#ifdef PMAP
  indexsize_nt = 3*indexsize;
#else
  indexsize_nt = indexsize;
#endif


  minactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));
  maxactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));

  if (totalpositions == 0) {
    debug(printf("Quitting because totalpositions is zero\n"));
    middle = (List_T) NULL;

  } else {
    debug(printf("Proceeding because maxnconsecutive is %d and pct_coverage is %f > %f or ncovered = %d > %d\n",
		 maxnconsecutive,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED));

    debug(printf("Performing diag on genomiclength %u\n",chrend-chrstart));
    Diag_max_bounds(minactive,maxactive,querylength,chrstart,chrend,chroffset,chrhigh,plusp);
    
    if ((all_paths = align_compute_lookback(mappings,npositions,totalpositions,
					    oned_matrix_p,minactive,maxactive,cellpool,
					    queryseq_ptr,queryuc_ptr,querylength,
					    /*querystart*/0,/*queryend*/querylength-1,
					    chroffset,chrhigh,plusp,
					    indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
					    /*anchoredp*/false,/*anchor_querypos*/0,/*anchor_position*/0,
					    localp,skip_repetitive_p,use_canonical_middle_p,NON_CANONICAL_PENALTY_MIDDLE,
					    favor_right_p,/*max_nalignments*/1,debug_graphic_p)) == NULL) {
      middle = (List_T) NULL;
    } else if ((path = (List_T) List_head(all_paths)) == NULL) {
      middle = (List_T) NULL;
    } else if (snps_p == true) {
      pairs = List_reverse(path);
      middle = convert_to_nucleotides_snps(pairs,
#ifndef PMAP
					   queryseq_ptr,queryuc_ptr,
#endif
					   chroffset,chrhigh,/*watsonp*/plusp,
					   query_offset,pairpool,indexsize_nt);
    } else {
      pairs = List_reverse(path);
      middle = convert_to_nucleotides(pairs,
#ifndef PMAP
				      queryseq_ptr,queryuc_ptr,
#endif
				      chroffset,chrhigh,/*watsonp*/plusp,
				      query_offset,pairpool,indexsize_nt);
    }

    List_free(&all_paths);
  }


  FREE(maxactive);
  FREE(minactive);
  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

#if 1
  for (source = 0; source < Oligoindex_array_length(oligoindices); source++) {
    oligoindex = Oligoindex_array_elt(oligoindices,source);
    Oligoindex_untally(oligoindex,queryuc_ptr,querylength);
  }
#endif

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
  for (p = diagonals; p != NULL; p = List_next(p)) {
    diag = (Diag_T) List_head(p);
    Diag_free(&diag);
  }
  List_free(&diagonals);
#endif

  return List_reverse(middle);
}

