static char rcsid[] = "$Id: bamtally.c 68053 2012-07-03 23:22:34Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bamtally.h"

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
#include "iit-write.h"
#include "interval.h"
#include "genome.h"
#include "table.h"
#include "ucharlist.h"


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

/* Bamtally_iit */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Tally_T
typedef struct T *T;



typedef enum {FIRST, SECOND} Pairend_T;



typedef struct Match_T *Match_T;
struct Match_T {
  int shift;			/* Used to record shifts */
  int mapq;
  char quality;
  long int count;
};

static Match_T
Match_new (int shift, int mapq, char quality) {
  Match_T new = (Match_T) MALLOC(sizeof(*new));

  new->shift = shift;
  new->mapq = mapq;
  new->quality = quality;
  new->count = 1;

  return new;
}

static void
Match_free (Match_T *old) {
  FREE(*old);
  return;
}

static Match_T
find_match_byshift (List_T matches, int shift) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->shift == shift) {
      return match;
    }
  }
  return (Match_T) NULL;
}

static Match_T
find_match_byquality (List_T matches, char quality) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->quality == quality) {
      return match;
    }
  }
  return (Match_T) NULL;
}

static Match_T
find_match_bymapq (List_T matches, char mapq) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->mapq == mapq) {
      return match;
    }
  }
  return (Match_T) NULL;
}

/* Go -1 to -readlength, then +readlength to +1 */
static int
Match_shift_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
Match_quality_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->quality > y->quality) {
    return -1;
  } else if (x->quality < y->quality) {
    return +1;
  } else {
    return 0;
  }
}

static int
Match_mapq_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->mapq > y->mapq) {
    return -1;
  } else if (x->mapq < y->mapq) {
    return +1;
  } else {
    return 0;
  }
}


typedef struct Mismatch_T *Mismatch_T;
struct Mismatch_T {
  char nt;
  int shift;			/* Used to record shifts */
  char mapq;
  char quality;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Mismatch_T next;		/* Used for linking similar mismatches together */
};


static Mismatch_T
Mismatch_new (char nt, int shift, int mapq, char quality) {
  Mismatch_T new = (Mismatch_T) MALLOC(sizeof(*new));

  new->nt = nt;
  new->shift = shift;
  new->mapq = mapq;
  new->quality = quality;
  new->count = 1;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;
  return new;
}

static void
Mismatch_free (Mismatch_T *old) {
  FREE(*old);
  return;
}

static int
Mismatch_chain_length (Mismatch_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}


static int
Mismatch_count_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return 0;
  }
}


/* Go -1 to -readlength, then +readlength to +1 */
static int
Mismatch_shift_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
Mismatch_quality_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->quality > y->quality) {
    return -1;
  } else if (x->quality < y->quality) {
    return +1;
  } else {
    return 0;
  }
}

static int
Mismatch_mapq_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->mapq > y->mapq) {
    return -1;
  } else if (x->mapq < y->mapq) {
    return +1;
  } else {
    return 0;
  }
}


static Mismatch_T
find_mismatch_byshift (List_T mismatches, char nt, int shift) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->shift == shift) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_byquality (List_T mismatches, char nt, char quality) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->quality == quality) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_bymapq (List_T mismatches, char nt, char mapq) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->mapq == mapq) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_nt (List_T mismatches, char nt) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}



typedef struct Insertion_T *Insertion_T;
struct Insertion_T {
  char *segment;
  int mlength;
  int shift;			/* Used to record shifts */
  int mapq;
  char quality;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Insertion_T next;		/* Used for linking similar insertions together */
};

static Insertion_T
Insertion_new (char *query_insert, int mlength, int shift, int mapq, char quality) {
  Insertion_T new = (Insertion_T) MALLOC(sizeof(*new));

  new->segment = (char *) CALLOC(mlength+1,sizeof(char));
  strncpy(new->segment,query_insert,mlength);
  new->mlength = mlength;
  new->shift = shift;
  new->mapq = mapq;
  new->quality = quality;
  new->count = 1;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;

  return new;
}

static void
Insertion_free (Insertion_T *old) {
  FREE((*old)->segment);
  FREE(*old);
  return;
}

static int
Insertion_chain_length (Insertion_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}


static int
Insertion_count_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return 0;
  }
}


/* Go -1 to -readlength, then +readlength to +1 */
static int
Insertion_shift_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static Insertion_T
find_insertion_seg (List_T insertions, char *segment, int mlength) {
  List_T p;
  Insertion_T ins;

  for (p = insertions; p != NULL; p = List_next(p)) {
    ins = (Insertion_T) List_head(p);
    if (ins->mlength == mlength && !strncmp(ins->segment,segment,mlength)) {
      return ins;
    }
  }
  return (Insertion_T) NULL;
}

static Insertion_T
find_insertion_byshift (List_T insertions, char *segment, int mlength, int shift) {
  List_T p;
  Insertion_T ins;

  for (p = insertions; p != NULL; p = List_next(p)) {
    ins = (Insertion_T) List_head(p);
    if (ins->shift == shift && ins->mlength == mlength && !strncmp(ins->segment,segment,mlength)) {
      return ins;
    }
  }
  return (Insertion_T) NULL;
}


typedef struct Deletion_T *Deletion_T;
struct Deletion_T {
  char *segment;
  int mlength;
  int shift;			/* Used to record shifts */
  int mapq;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Deletion_T next;		/* Used for linking similar deletions together */
};

static Deletion_T
Deletion_new (char *deletion, int mlength, int shift, int mapq) {
  Deletion_T new = (Deletion_T) MALLOC(sizeof(*new));

  new->segment = (char *) CALLOC(mlength+1,sizeof(char));
  strncpy(new->segment,deletion,mlength);
  new->mlength = mlength;
  new->shift = shift;
  new->mapq = mapq;
  new->count = 1;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;

  return new;
}

static void
Deletion_free (Deletion_T *old) {
  FREE((*old)->segment);
  FREE(*old);
  return;
}

static int
Deletion_chain_length (Deletion_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}


static int
Deletion_count_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return 0;
  }
}


/* Go -1 to -readlength, then +readlength to +1 */
static int
Deletion_shift_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static Deletion_T
find_deletion_seg (List_T deletions, char *segment, int mlength) {
  List_T p;
  Deletion_T del;

  for (p = deletions; p != NULL; p = List_next(p)) {
    del = (Deletion_T) List_head(p);
    /* Not necessary to check segment, since it is defined by starting position and mlength */
    if (del->mlength == mlength /* && !strncmp(del->segment,segment,mlength)*/ ) {
      return del;
    }
  }
  return (Deletion_T) NULL;
}

static Deletion_T
find_deletion_byshift (List_T deletions, char *segment, int mlength, int shift) {
  List_T p;
  Deletion_T del;

  for (p = deletions; p != NULL; p = List_next(p)) {
    del = (Deletion_T) List_head(p);
    if (del->shift == shift && del->mlength == mlength /* && !strncmp(del->segment,segment,mlength)*/) {
      return del;
    }
  }
  return (Deletion_T) NULL;
}



/************************************************************************
 *   Tallies
 ************************************************************************/

struct T {
  char refnt;
  int nmatches;
  List_T matches_byshift;
  List_T matches_byquality;
  List_T matches_bymapq;

  List_T mismatches_byshift;
  List_T mismatches_byquality;
  List_T mismatches_bymapq;

  List_T insertions_byshift;
  List_T deletions_byshift;
};

static T
Tally_new () {
  T new = (T) MALLOC(sizeof(*new));

  new->refnt = ' ';
  new->nmatches = 0;
  new->matches_byshift = (List_T) NULL;
  new->matches_byquality = (List_T) NULL;
  new->matches_bymapq = (List_T) NULL;

  new->mismatches_byshift = (List_T) NULL;
  new->mismatches_byquality = (List_T) NULL;
  new->mismatches_bymapq = (List_T) NULL;

  new->insertions_byshift = (List_T) NULL;
  new->deletions_byshift = (List_T) NULL;

  return new;
}


static void
Tally_clear (T this) {
  List_T ptr;
  Match_T match;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;

  this->refnt = ' ';
  this->nmatches = 0;

  for (ptr = this->matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&(this->matches_byshift));
  this->matches_byshift = (List_T) NULL;

  for (ptr = this->matches_byquality; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&(this->matches_byquality));
  this->matches_byquality = (List_T) NULL;

  for (ptr = this->matches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&(this->matches_bymapq));
  this->matches_bymapq = (List_T) NULL;

  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_byshift));
  this->mismatches_byshift = (List_T) NULL;

  for (ptr = this->mismatches_byquality; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_byquality));
  this->mismatches_byquality = (List_T) NULL;

  for (ptr = this->mismatches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_bymapq));
  this->mismatches_bymapq = (List_T) NULL;

  for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&(this->insertions_byshift));
  this->insertions_byshift = (List_T) NULL;

  for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&(this->deletions_byshift));
  this->deletions_byshift = (List_T) NULL;

  return;
}


static void
Tally_transfer (T dest, T src) {

  dest->refnt = src->refnt;
  src->refnt = ' ';

  dest->nmatches = src->nmatches;
  src->nmatches = 0;

  dest->matches_byshift = src->matches_byshift;
  src->matches_byshift = (List_T) NULL;

  dest->matches_byquality = src->matches_byquality;
  src->matches_byquality = (List_T) NULL;

  dest->matches_bymapq = src->matches_bymapq;
  src->matches_bymapq = (List_T) NULL;

  dest->mismatches_byshift = src->mismatches_byshift;
  src->mismatches_byshift = (List_T) NULL;

  dest->mismatches_byquality = src->mismatches_byquality;
  src->mismatches_byquality = (List_T) NULL;

  dest->mismatches_bymapq = src->mismatches_bymapq;
  src->mismatches_bymapq = (List_T) NULL;

  dest->insertions_byshift = src->insertions_byshift;
  src->insertions_byshift = (List_T) NULL;

  dest->deletions_byshift = src->deletions_byshift;
  src->deletions_byshift = (List_T) NULL;

  return;
}



static void
Tally_free (T *old) {
  List_T ptr;
  Match_T match;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;

  (*old)->refnt = ' ';
  (*old)->nmatches = 0;

  for (ptr = (*old)->matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->matches_byshift));
  (*old)->matches_byshift = (List_T) NULL;

  for (ptr = (*old)->matches_byquality; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->matches_byquality));
  (*old)->matches_byquality = (List_T) NULL;

  for (ptr = (*old)->matches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->matches_bymapq));
  (*old)->matches_bymapq = (List_T) NULL;

  for (ptr = (*old)->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_byshift));
  (*old)->mismatches_byshift = (List_T) NULL;

  for (ptr = (*old)->mismatches_byquality; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_byquality));
  (*old)->mismatches_byquality = (List_T) NULL;

  for (ptr = (*old)->mismatches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_bymapq));
  (*old)->mismatches_bymapq = (List_T) NULL;

  for (ptr = (*old)->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&((*old)->insertions_byshift));
  (*old)->insertions_byshift = (List_T) NULL;

  for (ptr = (*old)->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&((*old)->deletions_byshift));
  (*old)->deletions_byshift = (List_T) NULL;

  FREE(*old);
  *old = (T) NULL;

  return;
}



/************************************************************************
 *   Probability calculations
 ************************************************************************/

#define MAX_QUALITY_SCORE 40

#define NGENOTYPES 10

static char *genotype_string[NGENOTYPES] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
static bool genotype_heterozygousp[NGENOTYPES] = {false,true,true,true,false,true,true,false,true,false};


static double
allele_logprob[3][MAX_QUALITY_SCORE+1] = {
  /* Not in genotype: log(1/3*10^(-Q/10)) */
  {-1.098612,
   -1.328871, -1.559129, -1.789388, -2.019646, -2.249905,
   -2.480163, -2.710422, -2.940680, -3.170939, -3.401197,
   -3.631456, -3.861714, -4.091973, -4.322231, -4.552490,
   -4.782748, -5.013007, -5.243265, -5.473524, -5.703782,
   -5.934041, -6.164299, -6.394558, -6.624817, -6.855075,
   -7.085334, -7.315592, -7.545851, -7.776109, -8.006368,
   -8.236626, -8.466885, -8.697143, -8.927402, -9.157660,
   -9.387919, -9.618177, -9.848436, -10.078694, -10.308953},

  /* In heterozygote: log(1/2 - 1/3*10^(-Q/10)) */
  {-1.7917595,
   -1.4472174, -1.2389754, -1.0998002, -1.0015828, -0.9299061,
   -0.8764201, -0.8358837, -0.8048159, -0.7808079, -0.7621401,
   -0.7475561, -0.7361213, -0.7271306, -0.7200462, -0.7144544,
   -0.7100349, -0.7065382, -0.7037694, -0.7015754, -0.6998362,
   -0.6984568, -0.6973624, -0.6964940, -0.6958048, -0.6952576,
   -0.6948232, -0.6944782, -0.6942043, -0.6939868, -0.6938141,
   -0.6936769, -0.6935679, -0.6934814, -0.6934126, -0.6933580,
   -0.6933147, -0.6932802, -0.6932528, -0.6932311, -0.6932138},

  /* In homozygote: log(1 - 10^(-Q/10)) */
  {/*-Inf*/-1.58,
   -1.5814737534, -0.9968430440, -0.6955244713, -0.5076758737, -0.3801304081,
   -0.2892681872, -0.2225515160, -0.1725565729, -0.1345519603, -0.1053605157,
   -0.0827653027, -0.0651741732, -0.0514182742, -0.0406248442, -0.0321335740,
   -0.0254397275, -0.0201543648, -0.0159758692, -0.0126691702, -0.0100503359,
   -0.0079749983, -0.0063295629, -0.0050244739, -0.0039890173, -0.0031672882,
   -0.0025150465, -0.0019972555, -0.0015861505, -0.0012597185, -0.0010005003,
   -0.0007946439, -0.0006311565, -0.0005013129, -0.0003981864, -0.0003162778,
   -0.0002512202, -0.0001995461, -0.0001585019, -0.0001259005, -0.0001000050}
};


#if 0
static void
make_quality_scores_constant (int quality_score_constant) {
  int i, Q;
  double value;

  for (i = 0; i < 3; i++) {
    value = allele_logprob[i][quality_score_constant];
    for (Q = 0; Q <= MAX_QUALITY_SCORE; Q++) {
      allele_logprob[i][Q] = value;
    }
  }

  return;
}
#endif


static int
genotype_compatibility[4][NGENOTYPES] = {
  
  {/*AA*/ 2, /*AC*/ 1, /*AG*/ 1, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 0,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* C */
  {/*AA*/ 0, /*AC*/ 1, /*AG*/ 0, /*AT*/ 0,
   /*CC*/ 2, /*CG*/ 1, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* G */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 1, /*AT*/ 0,
   /*CC*/ 0, /*CG*/ 1, /*CT*/ 0,
   /*GG*/ 2, /*GT*/ 1,
   /*TT*/ 0},

  /* T */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 0, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 1,
   /*TT*/ 2}
};


static int
compute_probs (double *probs, double *loglik, char refnt,
	       List_T matches_byquality, List_T mismatches_byquality,
	       int quality_score_adj) {
  int bestg;
  double total, adj_loglik[NGENOTYPES], maxlik;
  List_T p;
  Match_T match;
  Mismatch_T mismatch;
  int Q;
  int g;

  for (g = 0; g < NGENOTYPES; g++) {
    loglik[g] = 0.0;
  }

  for (p = matches_byquality; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    Q = match->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (refnt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }

  for (p = mismatches_byquality; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    Q = mismatch->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (mismatch->nt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }
  
  /* Raise all loglikelihoods to maximum, to avoid underflow when taking exp() */
  maxlik = loglik[0];
  bestg = 0;
  for (g = 1; g < NGENOTYPES; g++) {
    if (loglik[g] > maxlik) {
      maxlik = loglik[g];
      bestg = g;
    } else if (loglik[g] == maxlik) {
      bestg = -1;
    }
  }

  for (g = 0; g < NGENOTYPES; g++) {
    adj_loglik[g] = loglik[g] - maxlik;
  }

  total = 0.0;
  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] = exp(adj_loglik[g]);
    total += probs[g];
  }

  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] /= total;
  }


  return bestg;
}




static bool
pass_variant_filter_p (long int nmatches, List_T mismatches_by_shift,
		       int min_depth, int variant_strands) {
  long int depth;
  bool plus_strand_p[4], minus_strand_p[4];
  List_T ptr;
  Mismatch_T mismatch;

  depth = nmatches;
  ptr = mismatches_by_shift;
  while (ptr != NULL) {
    mismatch = (Mismatch_T) List_head(ptr);
    depth += mismatch->count;
    ptr = List_next(ptr);
  }
  if (depth < min_depth) {
    return false;
  }


  if (variant_strands == 0) {
    return true;
  } else if (variant_strands == 1) {
    if (mismatches_by_shift == NULL) {
      return false;
    } else {
      return true;
    }

  } else {
    plus_strand_p[0] = plus_strand_p[1] = plus_strand_p[2] = plus_strand_p[3] = false;
    minus_strand_p[0] = minus_strand_p[1] = minus_strand_p[2] = minus_strand_p[3] = false;

    /* Look for two strands on a given mismatch allele */
    for (ptr = mismatches_by_shift; ptr != NULL; ptr = List_next(ptr)) {
      mismatch = (Mismatch_T) List_head(ptr);
      if (mismatch->shift > 0) {
	switch (mismatch->nt) {
	case 'A':
	  if (minus_strand_p[0] == true) {
	    return true;
	  } else {
	    plus_strand_p[0] = true;
	  }
	  break;
	case 'C':
	  if (minus_strand_p[1] == true) {
	    return true;
	  } else {
	    plus_strand_p[1] = true;
	  }
	  break;
	case 'G': 
	  if (minus_strand_p[2] == true) {
	    return true;
	  } else {
	    plus_strand_p[2] = true;
	  }
	  break;
	case 'T':
	  if (minus_strand_p[3] == true) {
	    return true;
	  } else {
	    plus_strand_p[3] = true;
	  }
	  break;
	}
      } else if (mismatch->shift < 0) {
	switch (mismatch->nt) {
	case 'A':
	  if (plus_strand_p[0] == true) {
	    return true;
	  } else {
	    minus_strand_p[0] = true;
	  }
	  break;
	case 'C':
	  if (plus_strand_p[1] == true) {
	    return true;
	  } else {
	    minus_strand_p[1] = true;
	  }
	  break;
	case 'G':
	  if (plus_strand_p[2] == true) {
	    return true;
	  } else {
	    minus_strand_p[2] = true;
	  }
	  break;
	case 'T':
	  if (plus_strand_p[3] == true) {
	    return true;
	  } else {
	    minus_strand_p[3] = true;
	  }
	  break;
	}
      }
    }

    return false;
  }
}


static void
print_allele_counts_simple (FILE *fp, T this, Genome_T genome, Genomicpos_T chroffset,
			    Genomicpos_T chrpos) {
  List_T p;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;

  if (this->insertions_byshift != NULL) {
    fprintf(fp,"\n");
    fprintf(fp,"^");
    for (p = this->insertions_byshift; p != NULL; p = List_next(p)) {
      ins = (Insertion_T) List_head(p);
      fprintf(fp,"%s %ld ",ins->segment,ins->count);
    }
    fprintf(fp,"\n");
  }

  if (this->deletions_byshift != NULL) {
    fprintf(fp,"\n");
    fprintf(fp,"_");
    for (p = this->deletions_byshift; p != NULL; p = List_next(p)) {
      del = (Deletion_T) List_head(p);
      fprintf(fp,"%s %ld ",ins->segment,del->count);
    }
    fprintf(fp,"\n");
  }

  if (this->nmatches == 0) {
    fprintf(fp,"%c0",Genome_get_char(genome,chroffset+chrpos-1U));
  } else {
    fprintf(fp,"%c%d",this->refnt,this->nmatches);
  }

  for (p = this->mismatches_byshift; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    fprintf(fp," %c%ld",mismatch->nt,mismatch->count);
  }

  return;
}


static bool
pass_difference_filter_p (double *probs, double *loglik, T this,
			  Genome_T genome, char *chr, Genomicpos_T chroffset, Genomicpos_T chrpos,
			  int quality_score_adj, bool genomic_diff_p) {
  int bestg;

  if (genomic_diff_p == false) {
    return true;
  } else {
    bestg = compute_probs(probs,loglik,this->refnt,this->matches_byquality,this->mismatches_byquality,
			  quality_score_adj);
    if (bestg < 0) {
      return false;
    } else if (genotype_heterozygousp[bestg] == true) {
      if (this->refnt == 'A') {
	if (genotype_compatibility[0][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  chr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'C') {
	if (genotype_compatibility[1][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  chr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'G') {
	if (genotype_compatibility[2][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  chr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else if (this->refnt == 'T') {
	if (genotype_compatibility[3][bestg] == 0) {
	  fprintf(stderr,"At %s:%u, genotype %s inconsistent with genome %c: ",
		  chr,chrpos,genotype_string[bestg],this->refnt);
	  print_allele_counts_simple(stderr,this,genome,chroffset,chrpos);
	  fprintf(stderr,"\n");
	  return false;
	} else {
	  return true;
	}
      } else {
	fprintf(stderr,"Reference nt not ACGT\n");
	return false;
      }

    } else if (this->refnt == 'A') {
      return (bestg == 0) ? false : true;
    } else if (this->refnt == 'C') {
      return (bestg == 4) ? false : true;
    } else if (this->refnt == 'G') {
      return (bestg == 7) ? false : true;
    } else if (this->refnt == 'T') {
      return (bestg == 9) ? false : true;
    } else {
      fprintf(stderr,"Reference nt not ACGT\n");
      return false;
    }
  }
}



static Genomicpos_T
print_runlength (T *block_tallies, Genomicpos_T *exonstart, Genomicpos_T lastpos,
		 Genomicpos_T blockstart, Genomicpos_T blockptr, 
		 char *chr) {
  T this;
  int blocki, lasti;
  Genomicpos_T chrpos;

  lasti = blockptr - blockstart;

  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (this->nmatches > 0 || this->mismatches_byshift != NULL ||
	this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
      chrpos = blockstart + blocki;
      if (chrpos > lastpos + 1U) {
	if (*exonstart == 0U) {
	  /* Skip printing */
	} else {
	  printf(">1 %s:%u..%u\n",chr,*exonstart,lastpos);
	}
	*exonstart = chrpos;
      }
      lastpos = chrpos;
    }
  }

  return lastpos;
}


static List_T
make_mismatches_unique (List_T mismatches) {
  List_T unique_mismatches = NULL, ptr;
  Mismatch_T mismatch, mismatch0;

  for (ptr = mismatches; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    if ((mismatch0 = find_mismatch_nt(unique_mismatches,mismatch->nt)) != NULL) {
      mismatch0->count += mismatch->count;

      /* Insert mismatch into list */
      mismatch->next = mismatch0->next;
      mismatch0->next = mismatch;

      mismatch0->shift += 1; /* Used here as nshifts */
    } else {
      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,/*mapq*/0,' ');
      mismatch0->count = mismatch->count;
      mismatch0->next = mismatch;
      unique_mismatches = List_push(unique_mismatches,mismatch0);
    }
  }

  return unique_mismatches;
}

static List_T
make_mismatches_unique_signed (List_T mismatches) {
  List_T unique_mismatches = NULL, ptr;
  Mismatch_T mismatch, mismatch0;

  for (ptr = mismatches; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    if ((mismatch0 = find_mismatch_nt(unique_mismatches,mismatch->nt)) != NULL) {
      if (mismatch->shift > 0) {
	mismatch0->count_plus += mismatch->count;
      } else {
	mismatch0->count_minus += mismatch->count;
      }

      /* Insert mismatch into list */
      mismatch->next = mismatch0->next;
      mismatch0->next = mismatch;

      mismatch0->shift += 1; /* Used here as nshifts */
    } else {
      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,/*mapq*/0,' ');
      if (mismatch->shift > 0) {
	mismatch0->count_plus = mismatch->count;
	mismatch0->count_minus = 0;
      } else {
	mismatch0->count_minus = mismatch->count;
	mismatch0->count_plus = 0;
      }
      mismatch0->next = mismatch;
      unique_mismatches = List_push(unique_mismatches,mismatch0);
    }
  }

  return unique_mismatches;
}



static void
print_block (T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr, 
	     Genome_T genome, char *chr, Genomicpos_T chroffset,
	     bool blockp, int quality_score_adj, int min_depth, int variant_strands,
	     bool genomic_diff_p, bool signed_counts_p,
	     bool print_totals_p, bool print_cycles_p,
	     bool print_quality_scores_p, bool print_mapq_scores_p, bool want_genotypes_p) {
  T this;
  double probs[NGENOTYPES], loglik[NGENOTYPES];
  int length, i, j, g, bestg;
  Genomicpos_T chrpos;
  int blocki, lasti;
  Match_T *match_array, match;
  Mismatch_T mismatch, mismatch0, *mm_array, *mm_subarray;
  Insertion_T ins, ins0, *ins_array, *ins_subarray;
  Deletion_T del, del0, *del_array, *del_subarray;
  List_T unique_mismatches_byshift, unique_mismatches_byquality, unique_mismatches_bymapq, ptr;
  List_T unique_insertions_byshift, unique_deletions_byshift;
  long int total, total_matches_plus, total_matches_minus, total_mismatches_plus, total_mismatches_minus;


  lasti = blockptr - blockstart;
  debug1(printf("Printing blocki 0 to %d\n",lasti));

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    total += this->nmatches;
    for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
      mismatch = (Mismatch_T) List_head(ptr);
      total += mismatch->count;
    }
  }
  
  /* Total could be 0 if the block is outside chrstart..chrend */
  if (total > 0) {
    if (blockp == true) {
      printf(">%ld %s:%u..%u\n",total,chr,blockstart,blockptr-1U);
    }

    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;


      /* Handle insertions */
      if (this->insertions_byshift != NULL) {
	unique_insertions_byshift = NULL;
	for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byshift,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

	    ins0->shift += 1; /* Used here as nshifts */
	  } else {
	    ins0 = Insertion_new(ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*mapq*/0,' ');
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byshift = List_push(unique_insertions_byshift,ins0);
	  }
	}

	ins_array = (Insertion_T *) List_to_array(unique_insertions_byshift,NULL);
	qsort(ins_array,List_length(unique_insertions_byshift),sizeof(Insertion_T),Insertion_count_cmp);

	printf("^");
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",chr,chrpos);
	}

	for (i = 0; i < List_length(unique_insertions_byshift); i++) {
	  ins0 = ins_array[i];
	  if (i > 0) {
	    printf(" ");
	  }
	  if (signed_counts_p == false) {
	    printf("%s %ld",ins0->segment,ins0->count_plus+ins0->count_minus);
	  } else {
	    printf("%s %ld|%ld",ins0->segment,ins0->count_plus,ins0->count_minus);
	  }
	  if (print_cycles_p == true) {
	    printf("(");
	    length = Insertion_chain_length(ins0->next);
	    ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	    for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
	      ins_subarray[j] = ins;
	    }
	    qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_shift_cmp);

	    printf("%ld@%d",ins_subarray[0]->count,ins_subarray[0]->shift);
	    for (j = 1; j < length; j++) {
	      printf(",%ld@%d",ins_subarray[j]->count,ins_subarray[j]->shift);
	    }
	    FREE(ins_subarray);
	    printf(")");
	  }
	}

	printf("\n");

	FREE(ins_array);

	for (ptr = unique_insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byshift);
      }


      /* Handle deletions */
      if (this->deletions_byshift != NULL) {
	unique_deletions_byshift = NULL;
	for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byshift,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

	    del0->shift += 1; /* Used here as nshifts */
	  } else {
	    del0 = Deletion_new(del->segment,del->mlength,/*shift, used here as nshifts*/1,/*mapq*/0);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byshift = List_push(unique_deletions_byshift,del0);
	  }
	}

	del_array = (Deletion_T *) List_to_array(unique_deletions_byshift,NULL);
	qsort(del_array,List_length(unique_deletions_byshift),sizeof(Deletion_T),Deletion_count_cmp);

	printf("_");
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",chr,chrpos);
	}

	for (i = 0; i < List_length(unique_deletions_byshift); i++) {
	  del0 = del_array[i];
	  if (i > 0) {
	    printf(" ");
	  }
	  if (signed_counts_p == false) {
	    printf("%s %ld",del0->segment,del0->count_plus+del0->count_minus);
	  } else {
	    printf("%s %ld|%ld",del0->segment,del0->count_plus,del0->count_minus);
	  }
	  if (print_cycles_p == true) {
	    printf("(");
	    length = Deletion_chain_length(del0->next);
	    del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	    for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
	      del_subarray[j] = del;
	    }
	    qsort(del_subarray,length,sizeof(Deletion_T),Deletion_shift_cmp);

	    printf("%ld@%d",del_subarray[0]->count,del_subarray[0]->shift);
	    for (j = 1; j < length; j++) {
	      printf(",%ld@%d",del_subarray[j]->count,del_subarray[j]->shift);
	    }
	    FREE(del_subarray);
	    printf(")");
	  }
	}

	printf("\n");

	FREE(del_array);

	for (ptr = unique_deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byshift);
      }


      /* Handle mismatches */

      if (pass_variant_filter_p(this->nmatches,this->mismatches_byshift,
				min_depth,variant_strands) == false) {
	if (blockp == true) {
	  printf("\n");
	}

      } else if (pass_difference_filter_p(probs,loglik,this,genome,
					  chr,chroffset,chrpos,
					  quality_score_adj,genomic_diff_p) == false) {
	if (blockp == true) {
	  printf("\n");
	}

      } else {
	if (blockp == false) {
	  /* Genomic position */
	  printf("%s\t%u\t",chr,chrpos);
	}

	if (signed_counts_p == false) {
	  total = this->nmatches;
	  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    total += mismatch->count;
	  }
	} else {
	  total_matches_plus = total_matches_minus = 0;
	  for (ptr = this->matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    match = (Match_T) List_head(ptr);
	    if (match->shift > 0) {
	      total_matches_plus += match->count;
	    } else {
	      total_matches_minus += match->count;
	    }
	  }
	  assert(this->nmatches == total_matches_plus + total_matches_minus);

	  total_mismatches_plus = total_mismatches_minus = 0;
	  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    if (mismatch->shift > 0) {
	      total_mismatches_plus += mismatch->count;
	    } else {
	      total_mismatches_minus += mismatch->count;
	    }
	  }
	}

	if (print_totals_p == true) {
	  /* Positional total */
	  if (signed_counts_p == false) {
	    printf("%ld\t",total);
	  } else {
	    printf("%ld|%ld\t",total_matches_plus+total_mismatches_plus,total_matches_minus+total_mismatches_minus);
	  }
	}

	/* Details */
	if (this->nmatches == 0) {
	  /* Genome_T expects zero-based numbers */
	  if (signed_counts_p == false) {
	    printf("%c0",Genome_get_char(genome,chroffset+chrpos-1U));
	  } else {
	    printf("%c0|0",Genome_get_char(genome,chroffset+chrpos-1U));
	  }
	} else {
	  if (signed_counts_p == false) {
	    printf("%c%d",this->refnt,this->nmatches);
	  } else {
	    printf("%c%ld|%ld",this->refnt,total_matches_plus,total_matches_minus);
	  }

	  if (print_cycles_p == true) {
	    printf("(");
	    length = List_length(this->matches_byshift);
	    match_array = (Match_T *) List_to_array(this->matches_byshift,NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_shift_cmp);
	    printf("%ld@%d",match_array[0]->count,match_array[0]->shift);
	    for (j = 1; j < length; j++) {
	      printf(",%ld@%d",match_array[j]->count,match_array[j]->shift);
	    }
	    FREE(match_array);
	    printf(")");
	  }

	  if (print_quality_scores_p == true) {
	    printf("(");
	    length = List_length(this->matches_byquality);
	    match_array = (Match_T *) List_to_array(this->matches_byquality,NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_quality_cmp);
	    printf("%ldQ%d",match_array[0]->count,match_array[0]->quality - quality_score_adj);
	    for (j = 1; j < length; j++) {
	      printf(",%ldQ%d",match_array[j]->count,match_array[j]->quality - quality_score_adj);
	    }
	    FREE(match_array);
	    printf(")");
	  }

	  if (print_mapq_scores_p == true) {
	    printf("(");
	    length = List_length(this->matches_bymapq);
	    match_array = (Match_T *) List_to_array(this->matches_bymapq,NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_mapq_cmp);
	    printf("%ldMAPQ%d",match_array[0]->count,match_array[0]->mapq);
	    for (j = 1; j < length; j++) {
	      printf(",%ldMAPQ%d",match_array[j]->count,match_array[j]->mapq);
	    }
	    FREE(match_array);
	    printf(")");
	  }
	}

	if (this->mismatches_byshift != NULL) {
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift);
	  unique_mismatches_byquality = make_mismatches_unique(this->mismatches_byquality);
	  unique_mismatches_bymapq = make_mismatches_unique(this->mismatches_bymapq);

	  mm_array = (Mismatch_T *) List_to_array(unique_mismatches_byshift,NULL);
	  qsort(mm_array,List_length(unique_mismatches_byshift),sizeof(Mismatch_T),Mismatch_count_cmp);

	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    mismatch0 = mm_array[i];
	    if (signed_counts_p == false) {
	      printf(" %c%ld",mismatch0->nt,mismatch0->count_plus+mismatch0->count_minus);
	    } else {
	      printf(" %c%ld|%ld",mismatch0->nt,mismatch0->count_plus,mismatch0->count_minus);
	    }
	    if (print_cycles_p == true) {
	      printf("(");
	      length = Mismatch_chain_length(mismatch0->next);
	      mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		mm_subarray[j] = mismatch;
	      }
	      qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_shift_cmp);

	      printf("%ld@%d",mm_subarray[0]->count,mm_subarray[0]->shift);
	      for (j = 1; j < length; j++) {
		printf(",%ld@%d",mm_subarray[j]->count,mm_subarray[j]->shift);
	      }
	      FREE(mm_subarray);
	      printf(")");
	    }
	      
	    if (print_quality_scores_p == true) {
	      printf("(");
	      mismatch0 = find_mismatch_nt(unique_mismatches_byquality,mismatch0->nt);
	      length = Mismatch_chain_length(mismatch0->next);
	      mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		mm_subarray[j] = mismatch;
	      }
	      qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_quality_cmp);
	      
	      printf("%ldQ%d",mm_subarray[0]->count,mm_subarray[0]->quality - quality_score_adj);
	      for (j = 1; j < length; j++) {
		printf(",%ldQ%d",mm_subarray[j]->count,mm_subarray[j]->quality - quality_score_adj);
	      }
	      FREE(mm_subarray);
	      printf(")");
	    }

	    if (print_mapq_scores_p == true) {
	      printf("(");
	      mismatch0 = find_mismatch_nt(unique_mismatches_bymapq,mismatch0->nt);
	      length = Mismatch_chain_length(mismatch0->next);
	      mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		mm_subarray[j] = mismatch;
	      }
	      qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_mapq_cmp);
	      
	      printf("%ldMAPQ%d",mm_subarray[0]->count,mm_subarray[0]->mapq);
	      for (j = 1; j < length; j++) {
		printf(",%ldMAPQ%d",mm_subarray[j]->count,mm_subarray[j]->mapq);
	      }
	      FREE(mm_subarray);
	      printf(")");
	    }
	  }

	  FREE(mm_array);

	  for (ptr = unique_mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byshift);

	  for (ptr = unique_mismatches_byquality; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byquality);

	  for (ptr = unique_mismatches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_bymapq);
	}

	if (want_genotypes_p == true) {
	  bestg = compute_probs(probs,loglik,this->refnt,this->matches_byquality,
				this->mismatches_byquality,quality_score_adj);

	  /* Best genotype */
	  printf("\t");
	  if (bestg < 0) {
	    printf("NN");
	  } else {
	    printf("%s",genotype_string[bestg]);
	  }

	  /* Probabilities */
	  printf("\t");
	  printf("%.3g",1.0-probs[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",1.0-probs[g]);
	  }

	  /* Log-likelihoods */
	  printf("\t");
	  printf("%.3g",loglik[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",loglik[g]);
	  }
	}

	printf("\n");
      }

      /* Reset */
      Tally_clear(this);
    }
  }
    
  return;
}


static void
tally_block (long int *tally_matches, long int *tally_mismatches,
	     T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr, 
	     Genome_T genome, char *chr, Genomicpos_T chroffset,
	     int quality_score_adj, int min_depth, int variant_strands, bool genomic_diff_p) {
  T this;
  double probs[NGENOTYPES], loglik[NGENOTYPES];
  Genomicpos_T chrpos;
  int blocki, lasti;
  Mismatch_T mismatch;
  List_T ptr;
  long int total;


  lasti = blockptr - blockstart;
  debug1(printf("Printing blocki 0 to %d\n",lasti));

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    total += this->nmatches;
    for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
      mismatch = (Mismatch_T) List_head(ptr);
      total += mismatch->count;
    }
  }
  
  /* Total could be 0 if the block is outside chrstart..chrend */
  if (total > 0) {
    for (blocki = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;

      if (pass_variant_filter_p(this->nmatches,this->mismatches_byshift,
				min_depth,variant_strands) == false) {
	/* Skip */

      } else if (pass_difference_filter_p(probs,loglik,this,genome,
					  chr,chroffset,chrpos,
					  quality_score_adj,genomic_diff_p) == false) {
	/* Skip */

      } else {
#if 0
	bini = (double) (chrpos - mincoord)/binstep;
	if (bini < 0) {
	  bini = 0;
	} else if (bini >= nbins) {
	  bini = nbins - 1;
	}
	binx_matches[bini] += this->nmatches;

	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    binx_mismatches[bini] += mismatch->count;
	}
#else
	tally_matches[chrpos] += this->nmatches;

	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    tally_mismatches[chrpos] += mismatch->count;
	}
#endif
      }

      /* Reset */
      Tally_clear(this);
    }
  }
    
  return;
}


static Ucharlist_T
push_char (int *nbytes, Ucharlist_T list, unsigned char x) {
  debug2(printf("%d: Pushing char %c\n",*nbytes,x));
  list = Ucharlist_push(list,x);
  *nbytes += 1;
  return list;
}

static Ucharlist_T
push_string (int *nbytes, Ucharlist_T list, char *string) {
  int length, i;

  debug2(printf("%d: Pushing char %s\n",*nbytes,string));
  length = strlen(string);
  for (i = 0; i < length; i++) {
    list = Ucharlist_push(list,string[i]);
  }
  list = Ucharlist_push(list,'\0');
  *nbytes += (length + 1);
  return list;
}

static Ucharlist_T
push_int (int *nbytes, Ucharlist_T list, int value) {
  unsigned char x;

  debug2(printf("%d: Pushing int %d\n",*nbytes,value));
  x = value & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);

  x = (value >>= 8) & 0xff;
  list = Ucharlist_push(list,x);
  
  *nbytes += 4;
  return list;
}


static void
iit_block (List_T *intervallist, List_T *labellist, List_T *datalist,
	   T *block_tallies, Genomicpos_T blockstart, Genomicpos_T blockptr,
	   Genome_T genome, Genomicpos_T chroffset, int quality_score_adj,
	   int min_depth, int variant_strands) {
  T this;
  char Buffer[100], *label;
  int length, i, j, k;
  Genomicpos_T chrpos;
  int blocki, lasti;
  Match_T *match_array, match;
  Mismatch_T mismatch, mismatch0, *mm_array, *mm_subarray;
  Insertion_T ins, ins0, *ins_array, *ins_subarray;
  Deletion_T del, del0, *del_array, *del_subarray;
  List_T unique_mismatches_byshift, unique_mismatches_byquality, unique_mismatches_bymapq, ptr;
  List_T unique_insertions_byshift, unique_deletions_byshift;
  long int total, total_matches_plus, total_matches_minus, total_mismatches_plus, total_mismatches_minus;
  int ninsertions, ndeletions;

  unsigned char *bytearray;
  Ucharlist_T pointers = NULL, bytes = NULL;
  int ignore, nbytes, ntotal;


  lasti = blockptr - blockstart;

  /* Block total */
  total = 0;
  for (blocki = 0; blocki < lasti; blocki++) {
    this = block_tallies[blocki];
    if (pass_variant_filter_p(this->nmatches,this->mismatches_byshift,min_depth,variant_strands) == true) {
      total += this->nmatches;
      for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	mismatch = (Mismatch_T) List_head(ptr);
	total += mismatch->count;
      }
    }
  }
  
  /* Total could be 0 if the block is outside chrstart..chrend */
  if (total > 0) {
    sprintf(Buffer,"%ld",total);
    label = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
    strcpy(label,Buffer);

    *intervallist = List_push(*intervallist,(void *) Interval_new(blockstart,blockptr-1U,/*sign*/+1));
    *labellist = List_push(*labellist,(void *) label);

    ignore = 0;
    nbytes = 0;
    for (blocki = 0, k = 0; blocki < lasti; blocki++) {
      this = block_tallies[blocki];
      chrpos = blockstart + blocki;
      debug2(printf("data for chrpos %u:\n",chrpos));

      /* Handle insertions at this position */
      pointers = push_int(&ignore,pointers,nbytes);
      if (this->insertions_byshift != NULL) {
	unique_insertions_byshift = NULL;
	for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins = (Insertion_T) List_head(ptr);
	  if ((ins0 = find_insertion_seg(unique_insertions_byshift,ins->segment,ins->mlength)) != NULL) {
	    if (ins->shift > 0) {
	      ins0->count_plus += ins->count;
	    } else {
	      ins0->count_minus += ins->count;
	    }

	    /* Insert insertion into list */
	    ins->next = ins0->next;
	    ins0->next = ins;

	    ins0->shift += 1; /* Used here as nshifts */
	  } else {
	    ins0 = Insertion_new(ins->segment,ins->mlength,/*shift, used here as nshifts*/1,/*mapq*/0,' ');
	    if (ins->shift > 0) {
	      ins0->count_plus = ins->count;
	      ins0->count_minus = 0;
	    } else {
	      ins0->count_minus = ins->count;
	      ins0->count_plus = 0;
	    }
	    ins0->next = ins;
	    unique_insertions_byshift = List_push(unique_insertions_byshift,ins0);
	  }
	}

	ins_array = (Insertion_T *) List_to_array(unique_insertions_byshift,NULL);
	ninsertions = List_length(unique_insertions_byshift);
	qsort(ins_array,ninsertions,sizeof(Insertion_T),Insertion_count_cmp);

	/* Total number of different insertions at this position */
	bytes = push_int(&nbytes,bytes,ninsertions);
	for (i = 0; i < ninsertions; i++) {
	  ins0 = ins_array[i];
	  /* Counts and segment for insertion i */
	  bytes = push_int(&nbytes,bytes,ins0->count_plus);
	  bytes = push_int(&nbytes,bytes,ins0->count_minus);
	  bytes = push_string(&nbytes,bytes,ins0->segment);

	  /* Cycles for insertion i (no quality) */
	  length = Insertion_chain_length(ins0->next);
	  ins_subarray = (Insertion_T *) CALLOC(length,sizeof(Insertion_T));
	  for (ins = ins0->next, j = 0; ins != NULL; ins = ins->next, j++) {
	    ins_subarray[j] = ins;
	  }
	  qsort(ins_subarray,length,sizeof(Insertion_T),Insertion_shift_cmp);

	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->shift);
	    bytes = push_int(&nbytes,bytes,ins_subarray[j]->count);
	  }
	  FREE(ins_subarray);
	}

	FREE(ins_array);

	for (ptr = unique_insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  ins0 = List_head(ptr);
	  Insertion_free(&ins0);
	}
	List_free(&unique_insertions_byshift);
      }

      /* Handle deletions at this position */
      pointers = push_int(&ignore,pointers,nbytes);
      if (this->deletions_byshift != NULL) {
	unique_deletions_byshift = NULL;
	for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del = (Deletion_T) List_head(ptr);
	  if ((del0 = find_deletion_seg(unique_deletions_byshift,del->segment,del->mlength)) != NULL) {
	    if (del->shift > 0) {
	      del0->count_plus += del->count;
	    } else {
	      del0->count_minus += del->count;
	    }

	    /* Insert deletion into list */
	    del->next = del0->next;
	    del0->next = del;

	    del0->shift += 1; /* Used here as nshifts */
	  } else {
	    del0 = Deletion_new(del->segment,del->mlength,/*shift, used here as nshifts*/1,/*mapq*/0);
	    if (del->shift > 0) {
	      del0->count_plus = del->count;
	      del0->count_minus = 0;
	    } else {
	      del0->count_minus = del->count;
	      del0->count_plus = 0;
	    }
	    del0->next = del;
	    unique_deletions_byshift = List_push(unique_deletions_byshift,del0);
	  }
	}

	del_array = (Deletion_T *) List_to_array(unique_deletions_byshift,NULL);
	ndeletions = List_length(unique_deletions_byshift);
	qsort(del_array,ndeletions,sizeof(Deletion_T),Deletion_count_cmp);

	/* Total number of different deletions at this position */
	bytes = push_int(&nbytes,bytes,ndeletions);
	for (i = 0; i < ndeletions; i++) {
	  del0 = del_array[i];
	  /* Counts and segment for deletion i */
	  bytes = push_int(&nbytes,bytes,del0->count_plus);
	  bytes = push_int(&nbytes,bytes,del0->count_minus);
	  bytes = push_string(&nbytes,bytes,del0->segment);

	  /* Cycles for deletion i (no quality) */
	  length = Deletion_chain_length(del0->next);
	  del_subarray = (Deletion_T *) CALLOC(length,sizeof(Deletion_T));
	  for (del = del0->next, j = 0; del != NULL; del = del->next, j++) {
	    del_subarray[j] = del;
	  }
	  qsort(del_subarray,length,sizeof(Deletion_T),Deletion_shift_cmp);

	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->shift);
	    bytes = push_int(&nbytes,bytes,del_subarray[j]->count);
	  }
	  FREE(del_subarray);
	}

	FREE(del_array);

	for (ptr = unique_deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  del0 = List_head(ptr);
	  Deletion_free(&del0);
	}
	List_free(&unique_deletions_byshift);
      }
      

      /* Handle allele counts at this position */
      pointers = push_int(&ignore,pointers,nbytes);
      if (pass_variant_filter_p(this->nmatches,this->mismatches_byshift,min_depth,variant_strands) == true) {
	total_matches_plus = total_matches_minus = 0;
	for (ptr = this->matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  match = (Match_T) List_head(ptr);
	  if (match->shift > 0) {
	    total_matches_plus += match->count;
	  } else {
	    total_matches_minus += match->count;
	  }
	}
	assert(this->nmatches == total_matches_plus + total_matches_minus);

	total_mismatches_plus = total_mismatches_minus = 0;
	for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  if (mismatch->shift > 0) {
	    total_mismatches_plus += mismatch->count;
	  } else {
	    total_mismatches_minus += mismatch->count;
	  }
	}

	/* Total signed counts at this position */
	bytes = push_int(&nbytes,bytes,total_matches_plus + total_mismatches_plus);
	bytes = push_int(&nbytes,bytes,total_matches_minus + total_mismatches_minus);
	
	/* Reference nucleotide and counts */
	if (this->nmatches == 0) {
	  bytes = push_char(&nbytes,bytes,Genome_get_char(genome,chroffset+chrpos-1U));
	  bytes = push_int(&nbytes,bytes,/*plus*/0);
	  bytes = push_int(&nbytes,bytes,/*minus*/0);
	} else {
	  bytes = push_char(&nbytes,bytes,this->refnt);
	  bytes = push_int(&nbytes,bytes,total_matches_plus);
	  bytes = push_int(&nbytes,bytes,total_matches_minus);
	}

	/* Alternate nucleotide and counts */
	if (this->mismatches_byshift != NULL) {
	  unique_mismatches_byshift = make_mismatches_unique_signed(this->mismatches_byshift);
	  unique_mismatches_byquality = make_mismatches_unique(this->mismatches_byquality);
	  unique_mismatches_bymapq = make_mismatches_unique(this->mismatches_bymapq);

	  mm_array = (Mismatch_T *) List_to_array(unique_mismatches_byshift,NULL);
	  qsort(mm_array,List_length(unique_mismatches_byshift),sizeof(Mismatch_T),Mismatch_count_cmp);

	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    mismatch0 = mm_array[i];
	    bytes = push_char(&nbytes,bytes,mismatch0->nt);
	    bytes = push_int(&nbytes,bytes,mismatch0->count_plus);
	    bytes = push_int(&nbytes,bytes,mismatch0->count_minus);
	  }
	}

	/* Terminates nucleotides */
	bytes = push_char(&nbytes,bytes,'\0');
	
	/* Cycles and quality for reference */
	if (this->nmatches > 0) {
	  length = List_length(this->matches_byshift);
	  match_array = (Match_T *) List_to_array(this->matches_byshift,NULL);
	  qsort(match_array,length,sizeof(Match_T),Match_shift_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,match_array[j]->shift);
	    bytes = push_int(&nbytes,bytes,match_array[j]->count);
	  }

	  length = List_length(this->matches_byquality);
	  match_array = (Match_T *) List_to_array(this->matches_byquality,NULL);
	  qsort(match_array,length,sizeof(Match_T),Match_quality_cmp);
	  bytes = push_int(&nbytes,bytes,length);
	  for (j = 0; j < length; j++) {
	    bytes = push_int(&nbytes,bytes,match_array[j]->quality - quality_score_adj);
	    bytes = push_int(&nbytes,bytes,match_array[j]->count);
	  }
	}
	    
	if (this->mismatches_byshift != NULL) {
	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    /* Cycles for alternate allele i */
	    mismatch0 = mm_array[i];
	    length = Mismatch_chain_length(mismatch0->next);
	    mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	    for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
	      mm_subarray[j] = mismatch;
	    }
	    qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_shift_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->shift);
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->count);
	    }
	    FREE(mm_subarray);

	    /* Quality scores for alternate allele i */
	    mismatch0 = find_mismatch_nt(unique_mismatches_byquality,mismatch0->nt);
	    length = Mismatch_chain_length(mismatch0->next);
	    mm_subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	    for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
	      mm_subarray[j] = mismatch;
	    }
	    qsort(mm_subarray,length,sizeof(Mismatch_T),Mismatch_quality_cmp);
	    bytes = push_int(&nbytes,bytes,length);
	    for (j = 0; j < length; j++) {
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->quality - quality_score_adj);
	      bytes = push_int(&nbytes,bytes,mm_subarray[j]->count);
	    }
	    FREE(mm_subarray);
	  }

	  for (ptr = unique_mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byshift);

	  for (ptr = unique_mismatches_byquality; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byquality);

	  for (ptr = unique_mismatches_bymapq; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_bymapq);
	}
      }

      /* Reset */
      Tally_clear(this);
    }

    /* Final pointer for block */
    pointers = push_int(&ignore,pointers,nbytes);

    /* Reverse lists before appending */
    pointers = Ucharlist_append(Ucharlist_reverse(pointers),Ucharlist_reverse(bytes));
    bytearray = Ucharlist_to_array(&ntotal,pointers);
    Ucharlist_free(&pointers);
    
    *datalist = List_push(*datalist,(void *) bytearray);
  }

  return;
}



static Genomicpos_T
transfer_position (T *alloc_tallies, T *block_tallies, Genomicpos_T *exonstart, Genomicpos_T lastpos,
		   Genomicpos_T *blockptr, Genomicpos_T *blockstart, Genomicpos_T *blockend,
		   long int *tally_matches, long int *tally_mismatches, 
		   List_T *intervallist, List_T *labellist, List_T *datalist,
		   Genomicpos_T chrpos, int alloci,
		   Genome_T genome, char *chr, Genomicpos_T chroffset,
		   Genomicpos_T chrstart, Genomicpos_T chrend, Tally_outputtype_T output_type,
		   bool blockp, int blocksize,
		   int quality_score_adj, int min_depth, int variant_strands,
		   bool genomic_diff_p, bool signed_counts_p,
		   bool print_totals_p, bool print_cycles_p, bool print_quality_scores_p,
		   bool print_mapq_scores_p, bool want_genotypes_p) {
  int blocki;

  if (chrpos < chrstart) {
    debug0(printf("    excluding chrpos %u < chrstart %u\n",chrpos,chrstart));
    Tally_clear(alloc_tallies[alloci]);
  } else if (chrpos > chrend) {
    debug0(printf("    excluding chrpos %u > chrend %u\n",chrpos,chrend));
    Tally_clear(alloc_tallies[alloci]);
  } else {
    if (chrpos >= *blockend) {
      debug0(printf("    chrpos %u >= blockend %u\n",chrpos,*blockend));
      debug0(printf("      print block from blockstart %u to blockptr %u\n",*blockstart,*blockptr));
      
      if (output_type == OUTPUT_RUNLENGTHS) {
	lastpos = print_runlength(block_tallies,&(*exonstart),lastpos,*blockstart,*blockptr,chr);
      } else if (output_type == OUTPUT_TALLY) {
	tally_block(tally_matches,tally_mismatches,
		    block_tallies,*blockstart,*blockptr,genome,chr,chroffset,
		    quality_score_adj,min_depth,variant_strands,genomic_diff_p);
      } else if (output_type == OUTPUT_IIT) {
	iit_block(&(*intervallist),&(*labellist),&(*datalist),
		  block_tallies,*blockstart,*blockptr,
		  genome,chroffset,quality_score_adj,min_depth,variant_strands);
      } else {
	print_block(block_tallies,*blockstart,*blockptr,genome,chr,chroffset,blockp,
		    quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		    signed_counts_p,print_totals_p,print_cycles_p,
		    print_quality_scores_p,print_mapq_scores_p,want_genotypes_p);
      }
      debug0(printf("      set blockstart to chrpos, blockend to chrpos + blocksize\n"));
      *blockstart = chrpos;
      *blockend = chrpos + blocksize;
    }

    blocki = chrpos - (*blockstart);
#if 0
    debug0(printf("      transfer position from alloci %d to blocki %d\n",alloci,blocki));
#endif

    Tally_transfer(block_tallies[blocki],alloc_tallies[alloci]);

    /* Points just beyond maximum chrpos observed */
    if (chrpos + 1U > *blockptr) {
      *blockptr = chrpos + 1U;
#if 0
      debug0(printf("    advancing blockptr to %u\n",*blockptr));
#endif
    }
  }

  return lastpos;
}


static void
transfer_position_lh (T *alloc_tallies_low, T *alloc_tallies_high, T *block_tallies_low, T *block_tallies_high,
		      Genomicpos_T *blockptr, Genomicpos_T *blockstart, Genomicpos_T *blockend,
		      long int *tally_matches_low, long int *tally_mismatches_low,
		      long int *tally_matches_high, long int *tally_mismatches_high,
		      Genomicpos_T chrpos, int alloci,
		      Genome_T genome, char *chr, Genomicpos_T chroffset,
		      Genomicpos_T chrstart, Genomicpos_T chrend, int blocksize,
		      int quality_score_adj, int min_depth, int variant_strands,
		      bool genomic_diff_p) {
  int blocki;

  if (chrpos < chrstart) {
    debug0(printf("    excluding chrpos %u < chrstart %u\n",chrpos,chrstart));
    Tally_clear(alloc_tallies_low[alloci]);
    Tally_clear(alloc_tallies_high[alloci]);
  } else if (chrpos > chrend) {
    debug0(printf("    excluding chrpos %u > chrend %u\n",chrpos,chrend));
    Tally_clear(alloc_tallies_low[alloci]);
    Tally_clear(alloc_tallies_high[alloci]);
  } else {
    if (chrpos >= *blockend) {
      debug0(printf("    chrpos %u >= blockend %u\n",chrpos,*blockend));
      debug0(printf("      print block from blockstart %u to blockptr %u\n",*blockstart,*blockptr));
      
      tally_block(tally_matches_low,tally_mismatches_low,
		  block_tallies_low,*blockstart,*blockptr,genome,chr,chroffset,
		  quality_score_adj,min_depth,variant_strands,genomic_diff_p);

      tally_block(tally_matches_high,tally_mismatches_high,
		  block_tallies_high,*blockstart,*blockptr,genome,chr,chroffset,
		  quality_score_adj,min_depth,variant_strands,genomic_diff_p);

      debug0(printf("      set blockstart to chrpos, blockend to chrpos + blocksize\n"));
      *blockstart = chrpos;
      *blockend = chrpos + blocksize;
    }

    blocki = chrpos - (*blockstart);
#if 0
    debug0(printf("      transfer position from alloci %d to blocki %d\n",alloci,blocki));
#endif

    Tally_transfer(block_tallies_low[blocki],alloc_tallies_low[alloci]);
    Tally_transfer(block_tallies_high[blocki],alloc_tallies_high[alloci]);

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




#if 0
static void
print_chromosome_blocks_wig (Genomicpos_T allocoffset, Genomicpos_T chroffset,
			     Genomicpos_T chrstart, Genomicpos_T chrend,
			     char *refnts, int *nmatches, 
			     List_T *matches_byshift, List_T *matches_byquality,
			     List_T *mismatches_byshift, List_T *mismatches_byquality,
			     Genome_T genome, char *chr) {
  Genomicpos_T pos, posi, firsti, lasti;
  bool hitp;
  Mismatch_T mismatch;
  List_T ptr;
  long int total;
 
  for (pos = chrstart - 1U; pos < chrend; pos += blocksize) {
    hitp = false;
    for (posi = pos; posi < chrend && posi < pos + blocksize; posi += 1U) {
      if (nmatches[allocoffset+posi] > 0 || mismatches_byshift[allocoffset+posi] != NULL) {
	if (hitp == false) {
	  firsti = posi;
	  hitp = true;
	}
	lasti = posi;
      }
    }
      
    if (hitp == true) {
      /* old: printf(">%s_%d %u %u %s\n",chromosome,++segno,firsti+1U,lasti+1U,chromosome); */
      /* Block total */
      total = 0;
      for (posi = firsti; posi <= lasti; posi++) {
	total += nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
      }
      printf("variableStep chrom=%s start=%u\n",chr,firsti+1U);

      for (posi = firsti; posi <= lasti; posi++) {
	printf("%u",posi+1U);

	/* Positional total */
	total = nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
	printf(" %ld",total);

	printf("\n");
      }
    }
  }

  return;
}
#endif



static void
revise_position (char querynt, char genomicnt, int mapq, int quality, int signed_shift,
		 T this, int *quality_counts_match, int *quality_counts_mismatch,
		 bool ignore_query_Ns_p) {
  Match_T match;
  Mismatch_T mismatch;

  if (toupper(querynt) == toupper(genomicnt)) {
    /* Count matches, even if query quality score was low */
    this->nmatches += 1;
    if ((match = find_match_byshift(this->matches_byshift,signed_shift)) != NULL) {
      match->count += 1;
    } else {
      this->matches_byshift = List_push(this->matches_byshift,(void *) Match_new(signed_shift,mapq,quality));
    }

    if ((match = find_match_byquality(this->matches_byquality,quality)) != NULL) {
      match->count += 1;
    } else {
      this->matches_byquality = List_push(this->matches_byquality,(void *) Match_new(signed_shift,mapq,quality));
    }

    if ((match = find_match_bymapq(this->matches_bymapq,mapq)) != NULL) {
      match->count += 1;
    } else {
      this->matches_bymapq = List_push(this->matches_bymapq,(void *) Match_new(signed_shift,mapq,quality));
    }

    quality_counts_match[quality] += 1;

  } else if (toupper(querynt) == 'N' && ignore_query_Ns_p == true) {
    /* Skip query N's */

  } else {
    if ((mismatch = find_mismatch_byshift(this->mismatches_byshift,toupper(querynt),signed_shift)) != NULL) {
      mismatch->count += 1;
    } else {
      this->mismatches_byshift = List_push(this->mismatches_byshift,(void *) Mismatch_new(toupper(querynt),signed_shift,mapq,quality));
    }

    if ((mismatch = find_mismatch_byquality(this->mismatches_byquality,toupper(querynt),quality)) != NULL) {
      mismatch->count += 1;
    } else {
      this->mismatches_byquality = List_push(this->mismatches_byquality,(void *) Mismatch_new(toupper(querynt),signed_shift,mapq,quality));
    }

    if ((mismatch = find_mismatch_bymapq(this->mismatches_bymapq,toupper(querynt),mapq)) != NULL) {
      mismatch->count += 1;
    } else {
      this->mismatches_bymapq = List_push(this->mismatches_bymapq,(void *) Mismatch_new(toupper(querynt),signed_shift,mapq,quality));
    }

    quality_counts_mismatch[quality] += 1;
  }

  return;
}


static void
revise_insertion (char *query_insert, int mlength, int mapq, int quality, int signed_shift,
		  T this) {
  Insertion_T ins;

  if ((ins = find_insertion_byshift(this->insertions_byshift,query_insert,mlength,signed_shift)) != NULL) {
    ins->count += 1;
  } else {
    this->insertions_byshift = List_push(this->insertions_byshift,(void *) Insertion_new(query_insert,mlength,signed_shift,mapq,quality));
  }

  /* Not yet doing find_insertion_byquality or find_insertion_bymapq */

  return;
}

static void
revise_deletion (char *deletion, int mlength, int mapq, int signed_shift, T this) {
  Deletion_T del;

  if ((del = find_deletion_byshift(this->deletions_byshift,deletion,mlength,signed_shift)) != NULL) {
    del->count += 1;
  } else {
    this->deletions_byshift = List_push(this->deletions_byshift,(void *) Deletion_new(deletion,mlength,signed_shift,mapq));
  }

  /* Not yet doing find_deletion_byquality or find_deletion_bymapq */

  return;
}



static double
average (char *quality_scores, int n) {
  double total = 0.0;
  int i;

  for (i = 0; i < n; i++) {
    total += (double) quality_scores[i];
  }
  return total/(double) n;
}



static void
revise_read (T *alloc_tallies, Genomicpos_T chrpos_low,
	     unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
	     char *shortread, int mapq, char *quality_string, Genomicpos_T alloc_low,
	     int *quality_counts_match, int *quality_counts_mismatch,
	     Genome_T genome, Genomicpos_T chroffset, bool ignore_query_Ns_p,
	     bool print_indels_p) {
  T this;
  int alloci;
  int shift, signed_shift;
  char strand;
  Genomicpos_T pos;
  char *genomic = NULL, *p, *q;
  char *r;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;
  int quality_score;


  debug1(printf("Revising read at %s:%u..%u\n",chr,chrpos_low,chrpos_high));

  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
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
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'I') {

      if (print_indels_p == true) {
	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing insertion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = alloc_tallies[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	if (quality_string == NULL) {
	  quality_score = DEFAULT_QUALITY;
	} else {
	  quality_score = (int) average(r,mlength);
	}

	revise_insertion(/*query_insert*/p,mlength,mapq,quality_score,signed_shift,this);
      }

      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);

    } else if (type == 'D') {

      if (print_indels_p == true) {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  /* q = &(*p); */
	  fprintf(stderr,"Need genome to print deletions\n");
	  exit(9);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing deletion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = alloc_tallies[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	revise_deletion(/*deletion*/q,mlength,mapq,signed_shift,this);
      }

      pos += mlength;

    } else if (type == 'M') {
      if (0 /* mlength < min_mlength */) {
	debug1(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
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
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));
	  this = alloc_tallies[alloci];
	  if (genome == NULL) {
	    this->refnt = ' ';

	  } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	    this->refnt = toupper(*q);
	    
	  } else if (this->refnt != toupper(*q)) {
	    fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		    this->refnt,*q,pos+1U);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    this->nmatches,List_length(this->mismatches_byshift));
	    exit(9);

	  }

	  signed_shift = (strand == '+') ? shift : -shift;
	  if (quality_string == NULL) {
	    quality_score = DEFAULT_QUALITY;
	  } else {
	    quality_score = (int) *r;
	  }
	  revise_position(/*querynt*/*p,/*genomicnt*/*q,mapq,quality_score,signed_shift,
			  this,quality_counts_match,quality_counts_mismatch,ignore_query_Ns_p);

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
  
  FREE(genomic);
  return;
}


static int
count_nsplices (Intlist_T types) {
  int nsplices = 0;
  Intlist_T p;

  for (p = types; p != NULL; p = Intlist_next(p)) {
    if (Intlist_head(p) == 'N') {
      nsplices++;
    }
  }
  return nsplices;
}


#if 0
/* Considers low and high relative to splice */
static void
revise_read_lh_old (T *alloc_tallies_low, T *alloc_tallies_high, bool lowend_p, Genomicpos_T chrpos_low,
		    unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
		    char *shortread, int mapq, char *quality_string, Genomicpos_T alloc_low,
		    int *quality_counts_match, int *quality_counts_mismatch,
		    Genome_T genome, Genomicpos_T chroffset, bool ignore_query_Ns_p) {
  T this;
  bool *lowp, *highp;
  int nsplices, splicei;
  int alloci;
  int shift, signed_shift;
  char strand;
  Genomicpos_T pos;
  char genomic[1024], *p, *q;
  char *r;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;
  int quality_score;


  nsplices = count_nsplices(types);
  lowp = (bool *) CALLOC(nsplices+1,sizeof(bool));
  highp = (bool *) CALLOC(nsplices+1,sizeof(bool));
  for (splicei = 0; splicei <= nsplices; splicei++) {
    lowp[splicei] = true;
    highp[splicei] = true;
  }
  if (lowend_p == true) {
    highp[0] = false;
  } else {
    lowp[nsplices] = false;
  }

  debug1(printf("Revising read at %s:%u..%u\n",chr,chrpos_low,chrpos_high));

  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
  }

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */

  splicei = 0;
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
      splicei++;
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
	debug1(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
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
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));

	  signed_shift = (strand == '+') ? shift : -shift;
	  if (quality_string == NULL) {
	    quality_score = DEFAULT_QUALITY;
	  } else {
	    quality_score = (int) *r;
	  }

	  if (lowp[splicei] == true) {
	    this = alloc_tallies_low[alloci];
	    if (genome == NULL) {
	      this->refnt = ' ';

	    } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	      this->refnt = toupper(*q);
	    
	    } else if (this->refnt != toupper(*q)) {
	      fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		      this->refnt,*q,pos+1U);
	      fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		      this->nmatches,List_length(this->mismatches_byshift));
	      exit(9);
	    }

	    revise_position(/*querynt*/*p,/*genomicnt*/*q,mapq,quality_score,signed_shift,
			    this,quality_counts_match,quality_counts_mismatch,ignore_query_Ns_p);
	  }

	  if (highp[splicei] == true) {
	    this = alloc_tallies_high[alloci];
	    if (genome == NULL) {
	      this->refnt = ' ';

	    } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	      this->refnt = toupper(*q);
	    
	    } else if (this->refnt != toupper(*q)) {
	      fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		      this->refnt,*q,pos+1U);
	      fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		      this->nmatches,List_length(this->mismatches_byshift));
	      exit(9);
	    }

	    revise_position(/*querynt*/*p,/*genomicnt*/*q,mapq,quality_score,signed_shift,
			    this,quality_counts_match,quality_counts_mismatch,ignore_query_Ns_p);
	  }

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
  
  FREE(highp);
  FREE(lowp);

  return;
}
#endif


static void
revise_read_lh (T *alloc_tallies_low, T *alloc_tallies_high, bool lowend_p, Genomicpos_T chrpos_low,
		unsigned int flag, Intlist_T types, Uintlist_T npositions, int cigar_querylength,
		char *shortread, int mapq, char *quality_string, Genomicpos_T alloc_low,
		int *quality_counts_match, int *quality_counts_mismatch,
		Genome_T genome, Genomicpos_T chroffset, bool ignore_query_Ns_p, bool print_indels_p) {
  T this;
  int alloci;
  int shift, signed_shift;
  char strand;
  Genomicpos_T pos;
  char *genomic = NULL, *p, *q;
  char *r;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;
  int quality_score;


  debug1(printf("Revising read at %s:%u..%u\n",chr,chrpos_low,chrpos_high));

  if (flag & QUERY_MINUSP) {
    strand = '-';
    shift = cigar_querylength;
  } else {
    strand = '+';
    shift = 1;
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

      if (print_indels_p == true) {
	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing insertion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	if (quality_string == NULL) {
	  quality_score = DEFAULT_QUALITY;
	} else {
	  quality_score = (int) average(r,mlength);
	}

	revise_insertion(/*query_insert*/p,mlength,mapq,quality_score,signed_shift,this);
      }

      p += mlength;
      r += mlength;
      shift += ((strand == '+') ? +mlength : -mlength);

    } else if (type == 'D') {

      if (print_indels_p == true) {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	if (genome == NULL) {
	  /* q = &(*p); */
	  fprintf(stderr,"Need genome to print deletions\n");
	  exit(9);
	} else {
	  FREE(genomic);
	  genomic = (char *) CALLOC(mlength+1,sizeof(char));
	  Genome_fill_buffer_simple(genome,/*left*/chroffset + pos,mlength,genomic);
	  /* printf("After (+): %s\n",genomic); */
	  q = genomic;
	}

	alloci = (pos + 1U) - alloc_low;
	debug1(printf("Processing deletion of length %d at shift %d, pos %u, alloci %d\n",
		      mlength,shift,pos+1U,alloci));

	this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];

	signed_shift = (strand == '+') ? shift : -shift;
	revise_deletion(/*deletion*/q,mlength,mapq,signed_shift,this);
      }

      pos += mlength;

    } else if (type == 'M') {
      if (0 /* mlength < min_mlength */) {
	debug1(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
	p += mlength;
	r += mlength;
	pos += mlength;
	shift += ((strand == '+') ? +mlength : -mlength);

      } else {
	debug1(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
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
	debug1(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

	while (mlength-- > /* trimright */ 0) {
	  alloci = (pos + 1U) - alloc_low;
	  debug1(printf("Processing %c and %c at shift %d, pos %u, mlength %u, alloci %d\n",
			*p,*q,shift,pos+1U,mlength,alloci));

	  signed_shift = (strand == '+') ? shift : -shift;
	  if (quality_string == NULL) {
	    quality_score = DEFAULT_QUALITY;
	  } else {
	    quality_score = (int) *r;
	  }

	  this = (lowend_p == true) ? alloc_tallies_low[alloci] : alloc_tallies_high[alloci];
	  if (genome == NULL) {
	    this->refnt = ' ';

	  } else if (this->nmatches == 0 && this->mismatches_byshift == NULL) {
	    this->refnt = toupper(*q);
	    
	  } else if (this->refnt != toupper(*q)) {
	    fprintf(stderr,"Two different genomic chars %c and %c at position %u\n",
		    this->refnt,*q,pos+1U);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    this->nmatches,List_length(this->mismatches_byshift));
	    exit(9);
	  }

	  revise_position(/*querynt*/*p,/*genomicnt*/*q,mapq,quality_score,signed_shift,
			  this,quality_counts_match,quality_counts_mismatch,ignore_query_Ns_p);

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
  
  FREE(genomic);
  return;
}


/* Also appears in gstruct.c */
static bool
best_mapping_p (Tableuint_T resolve_low_table, Tableuint_T resolve_high_table, char *acc,
		Genomicpos_T chrpos_low, Genomicpos_T chroffset) {
  Genomicpos_T genomicpos, genomicpos_low, genomicpos_high;

  genomicpos_low = Tableuint_get(resolve_low_table,(void *) acc);
  genomicpos_high = Tableuint_get(resolve_high_table,(void *) acc);


  if (genomicpos_low == 0 && genomicpos_high == 0) {
    /* Was already unique.  Could also check nhits. */
    return true;
  } else {
    genomicpos = chroffset + chrpos_low;
    if (genomicpos == genomicpos_low || genomicpos == genomicpos_high) {
      return true;
    } else {
      return false;
    }
  }
}


void
Bamtally_run (long int **tally_matches, long int **tally_mismatches,
	      List_T *intervallist, List_T *labellist, List_T *datalist,
	      int *quality_counts_match, int *quality_counts_mismatch,
	      Bamreader_T bamreader, Genome_T genome, char *chr,
	      Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, int alloclength,
	      Tableuint_T resolve_low_table, Tableuint_T resolve_high_table,
	      int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p,
	      bool ignore_lowend_p, bool ignore_highend_p,
	      Tally_outputtype_T output_type, bool blockp, int blocksize,
	      int quality_score_adj, int min_depth, int variant_strands,
	      bool genomic_diff_p, bool signed_counts_p, bool ignore_query_Ns_p,
	      bool print_indels_p, bool print_totals_p, bool print_cycles_p,
	      bool print_quality_scores_p, bool print_mapq_scores_p, bool want_genotypes_p,
	      bool verbosep) {
  T this;
  Genomicpos_T alloc_ptr, alloc_low, alloc_high, chrpos_low, chrpos_high, chrpos;
  Genomicpos_T blockptr = 0U, blockstart = 0U, blockend = 0U;
  int delta, alloci;
  bool goodp;
  Bamline_T bamline;
  int i;

  Genomicpos_T exonstart = 0U, lastpos = 0U;

  T *alloc_tallies, *block_tallies;


  alloc_tallies = (T *) CALLOC(alloclength,sizeof(T));
  for (i = 0; i < alloclength; i++) {
    alloc_tallies[i] = Tally_new();
  }

  block_tallies = (T *) CALLOC(blocksize,sizeof(T));
  for (i = 0; i < blocksize; i++) {
    block_tallies[i] = Tally_new();
  }

  *tally_matches = (long int *) NULL;
  *tally_mismatches = (long int *) NULL;
  *intervallist = (List_T) NULL;
  *labellist = (List_T) NULL;
  *datalist = (List_T) NULL;

  if (output_type == OUTPUT_TALLY) {
    *tally_matches = (long int *) CALLOC(chrend+1,sizeof(long int));
    *tally_mismatches = (long int *) CALLOC(chrend+1,sizeof(long int));
  }
  
  alloc_ptr = 0U;
  alloc_low = 0U;
  alloc_high = alloc_low + alloclength;
  goodp = false;

  while (goodp == false &&
	 (bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,need_concordant_p)) != NULL) {
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */
    debug0(printf(">%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Bamread_print_cigar(bamline));
    debug0(printf("\n"));

    if (resolve_low_table != NULL &&
	best_mapping_p(resolve_low_table,resolve_high_table,
		       Bamline_acc(bamline),Bamline_chrpos_low(bamline),
		       chroffset) == false) {
      /* Skip */
    } else {
      chrpos_low = Bamline_chrpos_low(bamline);
      chrpos_high = Bamline_chrpos_high(bamline);

      alloc_low = chrpos_low;
      alloc_high = alloc_low + alloclength;
      alloc_ptr = chrpos_low;

      blockstart = chrpos_low;
      blockend = blockstart + blocksize;
      blockptr = chrpos_low;

      debug0(printf("    initialize alloc_low %u, alloc_high = %u, blockstart = %u, blockend = %u\n",
		   alloc_low,alloc_high,blockstart,blockend));

      if (chrpos_high >= alloc_high) {
	if (verbosep == true) {
	  fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		  Bamline_acc(bamline),chr,chrpos_low,chrpos_high,alloclength);
	}
      } else {
	if (chrpos_high + 1U > alloc_ptr) {
	  alloc_ptr = chrpos_high + 1U;
	  debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	}

	revise_read(alloc_tallies,chrpos_low,Bamline_flag(bamline),
		    Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		    Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		    Bamline_mapq(bamline),Bamline_quality_string(bamline),alloc_low,
		    quality_counts_match,quality_counts_mismatch,
		    genome,chroffset,ignore_query_Ns_p,print_indels_p);

	goodp = true;
      }
    }

    Bamline_free(&bamline);
  }

  while ((bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,need_concordant_p)) != NULL) {
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */

    debug0(printf("*  alloc: %u..%u..%u  block: %u..%u..%u  ",
		  alloc_low,alloc_ptr,alloc_high,blockstart,blockptr,blockend));
    debug0(printf("%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Bamread_print_cigar(bamline));
    debug0(printf("\n"));
    
    if (resolve_low_table != NULL &&
	best_mapping_p(resolve_low_table,resolve_high_table,
		       Bamline_acc(bamline),Bamline_chrpos_low(bamline),
		       chroffset) == false) {
      /* Skip */
    } else if (ignore_lowend_p == true && Bamline_lowend_p(bamline) == true) {
      /* Skip */
    } else if (ignore_highend_p == true && Bamline_lowend_p(bamline) == false) {
      /* Skip */
    } else {
      chrpos_low = Bamline_chrpos_low(bamline);
      chrpos_high = Bamline_chrpos_high(bamline);

      if (chrpos_low > alloc_ptr) {
	debug0(printf("Case 1: chrpos_low %u > alloc_ptr %u\n",chrpos_low,alloc_ptr));
	debug0(printf("    transfer from alloc_low %u to alloc_ptr %u\n",alloc_low,alloc_ptr));
	for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
	  alloci = chrpos - alloc_low;
	  this = alloc_tallies[alloci];
	  if (this->nmatches > 0 || this->mismatches_byshift != NULL ||
	      this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
	    lastpos = transfer_position(alloc_tallies,block_tallies,&exonstart,lastpos,
					&blockptr,&blockstart,&blockend,
					*tally_matches,*tally_mismatches,
					&(*intervallist),&(*labellist),&(*datalist),
					chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
					output_type,blockp,blocksize,
					quality_score_adj,min_depth,variant_strands,genomic_diff_p,
					signed_counts_p,print_totals_p,print_cycles_p,
					print_quality_scores_p,print_mapq_scores_p,want_genotypes_p);
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
	  this = alloc_tallies[alloci];
	  if (this->nmatches > 0 || this->mismatches_byshift != NULL ||
	      this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
	    lastpos = transfer_position(alloc_tallies,block_tallies,&exonstart,lastpos,
					&blockptr,&blockstart,&blockend,
					*tally_matches,*tally_mismatches,
					&(*intervallist),&(*labellist),&(*datalist),
					chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
					output_type,blockp,blocksize,
					quality_score_adj,min_depth,variant_strands,genomic_diff_p,
					signed_counts_p,print_totals_p,print_cycles_p,
					print_quality_scores_p,print_mapq_scores_p,want_genotypes_p);
	  }
	}

	delta = chrpos_low - alloc_low;
	debug0(printf("    shift alloc downward by %d so alloc_low = chrpos_low\n",delta));
	for (alloci = 0; alloci < (int) (alloc_ptr - alloc_low - delta); alloci++) {
	  Tally_transfer(alloc_tallies[alloci],alloc_tallies[alloci+delta]);
	}
	for ( ; alloci < (int) (alloc_ptr - alloc_low); alloci++) {
	  debug1(printf("Resetting alloci %d\n",alloci));
	  Tally_clear(alloc_tallies[alloci]);
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

      if (chrpos_high >= alloc_high) {
	if (verbosep == true) {
	  fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		  Bamline_acc(bamline),chr,chrpos_low,chrpos_high,alloc_high);
	}
      } else {
	if (chrpos_high + 1U > alloc_ptr) {
	  alloc_ptr = chrpos_high + 1U;
	  debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
	}

	revise_read(alloc_tallies,chrpos_low,Bamline_flag(bamline),
		    Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		    Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		    Bamline_mapq(bamline),Bamline_quality_string(bamline),alloc_low,
		    quality_counts_match,quality_counts_mismatch,
		    genome,chroffset,ignore_query_Ns_p,print_indels_p);
      }
    }

    Bamline_free(&bamline);
  }

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));
  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
    alloci = chrpos - alloc_low;
    this = alloc_tallies[alloci];
    if (this->nmatches > 0 || this->mismatches_byshift != NULL ||
	this->insertions_byshift != NULL || this->deletions_byshift != NULL) {
      lastpos = transfer_position(alloc_tallies,block_tallies,&exonstart,lastpos,
				  &blockptr,&blockstart,&blockend,
				  *tally_matches,*tally_mismatches,
				  &(*intervallist),&(*labellist),&(*datalist),
				  chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
				  output_type,blockp,blocksize,
				  quality_score_adj,min_depth,variant_strands,genomic_diff_p,
				  signed_counts_p,print_totals_p,print_cycles_p,
				  print_quality_scores_p,print_mapq_scores_p,want_genotypes_p);
    }
  }

  debug0(printf("print block from blockstart %u to blockptr %u\n",blockstart,blockptr));
  if (output_type == OUTPUT_RUNLENGTHS) {
    lastpos = print_runlength(block_tallies,&exonstart,lastpos,blockstart,blockptr,chr);
  } else if (output_type == OUTPUT_TALLY) {
    tally_block(*tally_matches,*tally_mismatches,
		block_tallies,blockstart,blockptr,genome,chr,chroffset,
		quality_score_adj,min_depth,variant_strands,genomic_diff_p);
  } else if (output_type == OUTPUT_IIT) {
    iit_block(&(*intervallist),&(*labellist),&(*datalist),
	      block_tallies,blockstart,blockptr,
	      genome,chroffset,quality_score_adj,min_depth,variant_strands);
  } else{
    print_block(block_tallies,blockstart,blockptr,genome,chr,chroffset,blockp,
		quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		signed_counts_p,print_totals_p,print_cycles_p,
		print_quality_scores_p,print_mapq_scores_p,want_genotypes_p);
  }


  for (i = 0; i < alloclength; i++) {
    Tally_free(&(alloc_tallies[i]));
  }
  FREE(alloc_tallies);

  for (i = 0; i < blocksize; i++) {
    Tally_free(&(block_tallies[i]));
  }
  FREE(block_tallies);


  return;
}


/* Assumes output_type is OUTPUT_TALLY */
void
Bamtally_run_lh (long int **tally_matches_low, long int **tally_mismatches_low,
		 long int **tally_matches_high, long int **tally_mismatches_high,
		 int *quality_counts_match, int *quality_counts_mismatch,
		 Bamreader_T bamreader, Genome_T genome, char *chr,
		 Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, int alloclength,
		 int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		 bool need_concordant_p, bool need_unique_p, bool need_primary_p,
		 int blocksize, int quality_score_adj, int min_depth, int variant_strands,
		 bool genomic_diff_p, bool ignore_query_Ns_p, bool verbosep) {
  T this_low, this_high;
  Genomicpos_T alloc_ptr, alloc_low, alloc_high, chrpos_low, chrpos_high, chrpos;
  Genomicpos_T blockptr = 0U, blockstart = 0U, blockend = 0U;
  int delta, alloci;
  bool goodp;
  Bamline_T bamline;
  int i;

  T *alloc_tallies_low, *alloc_tallies_high, *block_tallies_low, *block_tallies_high;


  alloc_tallies_low = (T *) CALLOC(alloclength,sizeof(T));
  for (i = 0; i < alloclength; i++) {
    alloc_tallies_low[i] = Tally_new();
  }
  alloc_tallies_high = (T *) CALLOC(alloclength,sizeof(T));
  for (i = 0; i < alloclength; i++) {
    alloc_tallies_high[i] = Tally_new();
  }

  block_tallies_low = (T *) CALLOC(blocksize,sizeof(T));
  for (i = 0; i < blocksize; i++) {
    block_tallies_low[i] = Tally_new();
  }
  block_tallies_high = (T *) CALLOC(blocksize,sizeof(T));
  for (i = 0; i < blocksize; i++) {
    block_tallies_high[i] = Tally_new();
  }

  *tally_matches_low = (long int *) CALLOC(chrend+1,sizeof(long int));
  *tally_mismatches_low = (long int *) CALLOC(chrend+1,sizeof(long int));
  *tally_matches_high = (long int *) CALLOC(chrend+1,sizeof(long int));
  *tally_mismatches_high = (long int *) CALLOC(chrend+1,sizeof(long int));

  alloc_ptr = 0U;
  alloc_low = 0U;
  alloc_high = alloc_low + alloclength;
  goodp = false;

  while (goodp == false &&
	 (bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,need_concordant_p)) != NULL) {
    chrpos_low = Bamline_chrpos_low(bamline);
    chrpos_high = Bamline_chrpos_high(bamline);
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */

    debug0(printf(">%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Bamread_print_cigar(bamline));
    debug0(printf("\n"));

    alloc_low = chrpos_low;
    alloc_high = alloc_low + alloclength;
    alloc_ptr = chrpos_low;

    blockstart = chrpos_low;
    blockend = blockstart + blocksize;
    blockptr = chrpos_low;

    debug0(printf("    initialize alloc_low %u, alloc_high = %u, blockstart = %u, blockend = %u\n",
		  alloc_low,alloc_high,blockstart,blockend));

    if (chrpos_high >= alloc_high) {
      if (verbosep == true) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		Bamline_acc(bamline),chr,chrpos_low,chrpos_high,alloclength);
      }
    } else {
      if (chrpos_high + 1U > alloc_ptr) {
	alloc_ptr = chrpos_high + 1U;
	debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
      }

      revise_read_lh(alloc_tallies_low,alloc_tallies_high,Bamline_lowend_p(bamline),
		     chrpos_low,Bamline_flag(bamline),
		     Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		     Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		     Bamline_mapq(bamline),Bamline_quality_string(bamline),alloc_low,
		     quality_counts_match,quality_counts_mismatch,
		     genome,chroffset,ignore_query_Ns_p,/*print_indels_p*/false);

      goodp = true;
    }

    Bamline_free(&bamline);
  }

  while ((bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,need_concordant_p)) != NULL) {
    /* chrpos_high = Samread_chrpos_high(types,npositions,chrpos_low); */

    debug0(printf("*  alloc: %u..%u..%u  block: %u..%u..%u  ",
		  alloc_low,alloc_ptr,alloc_high,blockstart,blockptr,blockend));
    debug0(printf("%s:%u..%u ",chr,chrpos_low,chrpos_high));
    debug0(Bamread_print_cigar(bamline));
    debug0(printf("\n"));
    
    chrpos_low = Bamline_chrpos_low(bamline);
    chrpos_high = Bamline_chrpos_high(bamline);

    if (chrpos_low > alloc_ptr) {
      debug0(printf("Case 1: chrpos_low %u > alloc_ptr %u\n",chrpos_low,alloc_ptr));
      debug0(printf("    transfer from alloc_low %u to alloc_ptr %u\n",alloc_low,alloc_ptr));
      for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
	alloci = chrpos - alloc_low;
	this_low = alloc_tallies_low[alloci];
	this_high = alloc_tallies_high[alloci];
	if (this_low->nmatches > 0 || this_low->mismatches_byshift != NULL ||
	    this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	    this_high->nmatches > 0 || this_high->mismatches_byshift != NULL ||
	    this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
	  transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			       &blockptr,&blockstart,&blockend,
			       *tally_matches_low,*tally_mismatches_low,
			       *tally_matches_high,*tally_mismatches_high,
			       chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
			       blocksize,
			       quality_score_adj,min_depth,variant_strands,genomic_diff_p);
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
	this_low = alloc_tallies_low[alloci];
	this_high = alloc_tallies_high[alloci];
	if (this_low->nmatches > 0 || this_low->mismatches_byshift != NULL ||
	    this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	    this_high->nmatches > 0 || this_high->mismatches_byshift != NULL ||
	    this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
	  transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			       &blockptr,&blockstart,&blockend,
			       *tally_matches_low,*tally_mismatches_low,
			       *tally_matches_high,*tally_mismatches_high,
			       chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
			       blocksize,
			       quality_score_adj,min_depth,variant_strands,genomic_diff_p);
	}
      }

      delta = chrpos_low - alloc_low;
      debug0(printf("    shift alloc downward by %d so alloc_low = chrpos_low\n",delta));
      for (alloci = 0; alloci < (int) (alloc_ptr - alloc_low - delta); alloci++) {
	Tally_transfer(alloc_tallies_low[alloci],alloc_tallies_low[alloci+delta]);
	Tally_transfer(alloc_tallies_high[alloci],alloc_tallies_high[alloci+delta]);
      }
      for ( ; alloci < (int) (alloc_ptr - alloc_low); alloci++) {
	debug1(printf("Resetting alloci %d\n",alloci));
	Tally_clear(alloc_tallies_low[alloci]);
	Tally_clear(alloc_tallies_high[alloci]);
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

    if (chrpos_high >= alloc_high) {
      if (verbosep == true) {
	fprintf(stderr,"read %s at %s:%u..%u is longer than allocated buffer ending at %u => skipping\n",
		Bamline_acc(bamline),chr,chrpos_low,chrpos_high,alloc_high);
      }
    } else {
      if (chrpos_high + 1U > alloc_ptr) {
	alloc_ptr = chrpos_high + 1U;
	debug0(printf("    revising alloc_ptr to be %u\n",alloc_ptr));
      }

      revise_read_lh(alloc_tallies_low,alloc_tallies_high,Bamline_lowend_p(bamline),
		     chrpos_low,Bamline_flag(bamline),
		     Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
		     Bamline_cigar_querylength(bamline),Bamline_read(bamline),
		     Bamline_mapq(bamline),Bamline_quality_string(bamline),alloc_low,
		     quality_counts_match,quality_counts_mismatch,
		     genome,chroffset,ignore_query_Ns_p,/*print_indels_p*/false);
    }

    Bamline_free(&bamline);
  }

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));

  debug0(printf("end of reads\n"));
  debug0(printf("    transfer from alloc_low %u up to alloc_ptr %u\n",alloc_low,alloc_ptr));
  for (chrpos = alloc_low; chrpos < alloc_ptr; chrpos++) {
    alloci = chrpos - alloc_low;
    this_low = alloc_tallies_low[alloci];
    this_high = alloc_tallies_high[alloci];
    if (this_low->nmatches > 0 || this_low->mismatches_byshift != NULL ||
	this_low->insertions_byshift != NULL || this_low->deletions_byshift != NULL ||
	this_high->nmatches > 0 || this_high->mismatches_byshift != NULL ||
	this_high->insertions_byshift != NULL || this_high->deletions_byshift != NULL) {
      transfer_position_lh(alloc_tallies_low,alloc_tallies_high,block_tallies_low,block_tallies_high,
			   &blockptr,&blockstart,&blockend,
			   *tally_matches_low,*tally_mismatches_low,
			   *tally_matches_high,*tally_mismatches_high,
			   chrpos,alloci,genome,chr,chroffset,chrstart,chrend,
			   blocksize,
			   quality_score_adj,min_depth,variant_strands,genomic_diff_p);
    }
  }

  debug0(printf("print block from blockstart %u to blockptr %u\n",blockstart,blockptr));
  tally_block(*tally_matches_low,*tally_mismatches_low,
	      block_tallies_low,blockstart,blockptr,genome,chr,chroffset,
	      quality_score_adj,min_depth,variant_strands,genomic_diff_p);
  tally_block(*tally_matches_high,*tally_mismatches_high,
	      block_tallies_high,blockstart,blockptr,genome,chr,chroffset,
	      quality_score_adj,min_depth,variant_strands,genomic_diff_p);


  for (i = 0; i < alloclength; i++) {
    Tally_free(&(alloc_tallies_low[i]));
    Tally_free(&(alloc_tallies_high[i]));
  }
  FREE(alloc_tallies_low);
  FREE(alloc_tallies_high);

  for (i = 0; i < blocksize; i++) {
    Tally_free(&(block_tallies_low[i]));
    Tally_free(&(block_tallies_high[i]));
  }
  FREE(block_tallies_low);
  FREE(block_tallies_high);


  return;
}



IIT_T
Bamtally_iit (Bamreader_T bamreader, char *desired_chr, Genomicpos_T chrstart, Genomicpos_T chrend,
	      Genome_T genome, IIT_T chromosome_iit, int alloclength,
	      int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p,
	      int min_depth, int variant_strands, bool ignore_query_Ns_p,
	      bool print_indels_p, int blocksize, bool verbosep) {
  IIT_T iit;

  long int *tally_matches, *tally_mismatches;
  List_T divlist = NULL, typelist = NULL;
  List_T intervallist, labellist, datalist, p;
  Interval_T interval;
  int quality_counts_match[256], quality_counts_mismatch[256];
  char *divstring, *typestring, *label;
  char *chr;

  int index;
  Genomicpos_T chroffset, chrlength;
  Table_T intervaltable, labeltable, datatable;
  bool allocp;

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  datatable = Table_new(65522,Table_string_compare,Table_string_hash);

  if (desired_chr == NULL) {
    /* Entire genome */
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      fprintf(stderr,"  processing chromosome %s...",chr);
      divlist = List_push(divlist,(void *) chr);

      chroffset = Interval_low(IIT_interval(chromosome_iit,index));
      chrlength = Interval_length(IIT_interval(chromosome_iit,index));

      Bamread_limit_region(bamreader,chr,/*chrstart*/1,/*chrend*/chrlength);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,chr,chroffset,/*chrstart*/1U,/*chrend*/chrlength,alloclength,
		   /*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_IIT,/*blockp*/false,blocksize,
		   /*quality_score_adj*/0,min_depth,variant_strands,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,ignore_query_Ns_p,
		   print_indels_p,/*print_totals_p*/false,/*print_cycles_p*/false,
		   /*print_quality_scores_p*/false,/*print_mapq_scores_p*/false,
		   /*want_genotypes_p*/false,verbosep);
      /* Reverse lists so we can specify presortedp == true */
      Table_put(intervaltable,(void *) chr,List_reverse(intervallist));
      Table_put(labeltable,(void *) chr,List_reverse(labellist));
      Table_put(datatable,(void *) chr,List_reverse(datalist));
      Bamread_unlimit_region(bamreader);
      /* Don't free chr yet, since divlist contains it */
      fprintf(stderr,"done\n");
    }

  } else if ((index = IIT_find_linear(chromosome_iit,desired_chr)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in genome\n",desired_chr);
    Table_free(&datatable);
    Table_free(&labeltable);
    Table_free(&intervaltable);

    return NULL;

  } else {
    /* Single chromosomal region */
    divlist = List_push(divlist,(void *) desired_chr);
    
    index = IIT_find_linear(chromosome_iit,desired_chr);
    chroffset = Interval_low(IIT_interval(chromosome_iit,index));
    
    Bamread_limit_region(bamreader,desired_chr,chrstart,chrend);
    Bamtally_run(&tally_matches,&tally_mismatches,
		 &intervallist,&labellist,&datalist,
		 quality_counts_match,quality_counts_mismatch,
		 bamreader,genome,chr,chroffset,/*chrstart*/chrstart,/*chrend*/chrend,alloclength,
		 /*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		 minimum_mapq,good_unique_mapq,maximum_nhits,
		 need_concordant_p,need_unique_p,need_primary_p,
		 /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		 /*output_type*/OUTPUT_IIT,/*blockp*/false,blocksize,
		 /*quality_score_adj*/0,min_depth,variant_strands,
		 /*genomic_diff_p*/false,/*signed_counts_p*/false,ignore_query_Ns_p,
		 print_indels_p,/*print_totals_p*/false,/*print_cycles_p*/false,
		 /*print_quality_scores_p*/false,/*print_mapq_scores_p*/false,
		 /*want_genotypes_p*/false,verbosep);
    /* Reverse lists so we can specify presortedp == true */
    Table_put(intervaltable,(void *) desired_chr,List_reverse(intervallist));
    Table_put(labeltable,(void *) desired_chr,List_reverse(labellist));
    Table_put(datatable,(void *) desired_chr,List_reverse(datalist));
    Bamread_unlimit_region(bamreader);
  }


  divlist = List_reverse(divlist);

  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(divlist,divstring);

  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(NULL,typestring);

  iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,intervaltable,
		   labeltable,datatable,/*divsort*/NO_SORT,/*version*/IIT_LATEST_VERSION,
		   /*presortedp*/true);

  FREE(typestring);
  List_free(&typelist);
  FREE(divstring);
  /* May need to free chr inside of divlist */
  List_free(&divlist);

  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      datalist = Table_get(datatable,(void *) chr);
      /* Do not free data yet, since tally_iit points to it */
      List_free(&datalist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    datalist = Table_get(datatable,(void *) desired_chr);
    /* Do not free data yet, since tally_iit points to it */
    List_free(&datalist);
  }
  Table_free(&datatable);


  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      labellist = Table_get(labeltable,(void *) chr);
      for (p = labellist; p != NULL; p = List_next(p)) {
	label = (char *) List_head(p);
	FREE(label);
      }
      List_free(&labellist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    labellist = Table_get(labeltable,(void *) desired_chr);
    for (p = labellist; p != NULL; p = List_next(p)) {
      label = (char *) List_head(p);
      FREE(label);
    }
    List_free(&labellist);
  }
  Table_free(&labeltable);
  

  if (desired_chr == NULL) {
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      chr = IIT_label(chromosome_iit,index,&allocp);
      intervallist = Table_get(intervaltable,(void *) chr);
      for (p = intervallist; p != NULL; p = List_next(p)) {
	interval = (Interval_T) List_head(p);
	Interval_free(&interval);
      }
      List_free(&intervallist);
      if (allocp) {
	FREE(chr);
      }
    }
  } else {
    intervallist = Table_get(intervaltable,(void *) desired_chr);
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&intervallist);
  }
  Table_free(&intervaltable);

  return iit;
}

