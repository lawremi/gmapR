static char rcsid[] = "$Id: tally.c 159976 2015-03-02 22:51:51Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tally.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"


#define INITIAL_READLENGTH 75


/* Matchpool and mismatchpool are buggy, resulting in a loop somewhere */
/* #define USE_MATCHPOOL 1 */
/* #define USE_MISMATCHPOOL 1 */


#ifndef USE_MATCHPOOL
Match_T
Match_new (int shift, int nm, int xs, int ncounts) {
  Match_T new = (Match_T) MALLOC(sizeof(*new));

  new->shift = shift;
  new->nm = nm;
  new->xs = xs;			/* intron strand (+1 for '+', +2 for '-', 0 for others) */
  new->count = ncounts;

  return new;
}

void
Match_free (Match_T *old) {
  FREE(*old);
  return;
}
#endif



#ifndef USE_MISMATCHPOOL
Mismatch_T
Mismatch_new (char nt, int shift, int nm, int xs, int ncounts) {
  Mismatch_T new = (Mismatch_T) MALLOC(sizeof(*new));

  new->nt = nt;
  new->shift = shift;
  new->nm = nm;
  new->xs = xs;
  new->count = ncounts;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;
  return new;
}

void
Mismatch_free (Mismatch_T *old) {
  FREE(*old);
  return;
}
#endif


Insertion_T
Insertion_new (Genomicpos_T chrpos, char *query_insert, int mlength, int shift, int nm, int ncounts) {
  Insertion_T new = (Insertion_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->segment = (char *) CALLOC(mlength+1,sizeof(char));
  strncpy(new->segment,query_insert,mlength);
  new->mlength = mlength;
  new->shift = shift;
  new->nm = nm;
  new->count = ncounts;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;

  return new;
}

void
Insertion_free (Insertion_T *old) {
  FREE((*old)->segment);
  FREE(*old);
  return;
}

int
Insertion_count_cmp (const void *a, const void *b) {
  Insertion_T x = * (Insertion_T *) a;
  Insertion_T y = * (Insertion_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return strcmp(x->segment,y->segment);
  }
}

Insertion_T
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

Insertion_T
find_insertion_bynm (List_T insertions, char *segment, int mlength, int nm) {
  List_T p;
  Insertion_T ins;

  for (p = insertions; p != NULL; p = List_next(p)) {
    ins = (Insertion_T) List_head(p);
    if (ins->nm == nm && ins->mlength == mlength && !strncmp(ins->segment,segment,mlength)) {
      return ins;
    }
  }
  return (Insertion_T) NULL;
}


Insertion_T
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

Deletion_T
Deletion_new (Genomicpos_T chrpos, char *deletion, int mlength, int shift, int nm, int ncounts) {
  Deletion_T new = (Deletion_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->segment = (char *) CALLOC(mlength+1,sizeof(char));
  strncpy(new->segment,deletion,mlength);
  new->mlength = mlength;
  new->shift = shift;
  new->nm = nm;
  new->count = ncounts;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;

  return new;
}

void
Deletion_free (Deletion_T *old) {
  FREE((*old)->segment);
  FREE(*old);
  return;
}


int
Deletion_count_cmp (const void *a, const void *b) {
  Deletion_T x = * (Deletion_T *) a;
  Deletion_T y = * (Deletion_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return strcmp(x->segment,y->segment);
  }
}

Deletion_T
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

Deletion_T
find_deletion_bynm (List_T deletions, char *segment, int mlength, int nm) {
  List_T p;
  Deletion_T del;

  for (p = deletions; p != NULL; p = List_next(p)) {
    del = (Deletion_T) List_head(p);
    if (del->nm == nm && del->mlength == mlength /* && !strncmp(del->segment,segment,mlength)*/) {
      return del;
    }
  }
  return (Deletion_T) NULL;
}



Deletion_T
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


/************************************************************************
 *   Readevid_T
 ************************************************************************/

struct Readevid_T {
  unsigned int linei;
  char nt;
  char nti;

  int shift;
  int nm;
  int xs;
};


Readevid_T
Readevid_new (unsigned int linei, char nt, int shift, int nm, int xs) {
  Readevid_T new = (Readevid_T) MALLOC(sizeof(*new));

  new->linei = linei;
  new->nt = nt;
  switch (nt) {
  case 'A': new->nti = 0; break;
  case 'C': new->nti = 1; break;
  case 'G': new->nti = 2; break;
  case 'T': new->nti = 3; break;
  default: new->nti = -1;
  }

  new->shift = shift;
  new->nm = nm;
  new->xs = xs;

  return new;
}


static void
Readevid_free (Readevid_T *old) {
  FREE(*old);
  return;
}

unsigned int
Readevid_linei (Readevid_T this) {
  return this->linei;
}

char
Readevid_nt (Readevid_T this) {
  return this->nt;
}


char
Readevid_codoni_plus (int *shift, int *nm, int *xs,
		      Readevid_T frame0, Readevid_T frame1, Readevid_T frame2) {
  if (frame0->nti < 0) {
    return -1;
  } else if (frame1->nti < 0) {
    return -1;
  } else if (frame2->nti < 0) {
    return -1;

  } else {
    *shift = frame0->shift;
    if (frame1->shift < *shift) {
      *shift = frame1->shift;
    }
    if (frame2->shift < *shift) {
      *shift = frame2->shift;
    }

#if 0
    /* MAPQ should be a read-level quantity */
    if (frame1->mapq < *mapq) {
      *mapq = frame1->mapq;
    }
    if (frame2->mapq < *mapq) {
      *mapq = frame2->mapq;
    }
#endif

    *nm = frame0->nm;
    *xs = frame0->xs;
#if 0
    /* XS should be a read-level quantity */
    if (frame1->xs != *xs) {
      *xs = 0;
    } else if (frame2->xs != *xs) {
      *xs = 0;
    }
#endif

    return 16*frame0->nti + 4*frame1->nti + frame2->nti;
  }
}

char
Readevid_codoni_minus (int *shift, int *nm, int *xs,
		       Readevid_T frame0, Readevid_T frame1, Readevid_T frame2) {
  if (frame0->nti < 0) {
    return -1;
  } else if (frame1->nti < 0) {
    return -1;
  } else if (frame2->nti < 0) {
    return -1;
  } else {
    *shift = frame0->shift;
    if (frame1->shift < *shift) {
      *shift = frame1->shift;
    }
    if (frame2->shift < *shift) {
      *shift = frame2->shift;
    }

#if 0
    /* MAPQ should be a read-level quantity */
    if (frame1->mapq < *mapq) {
      *mapq = frame1->mapq;
    }
    if (frame2->mapq < *mapq) {
      *mapq = frame2->mapq;
    }
#endif

    *nm = frame0->nm;
    *xs = frame0->xs;
#if 0
    /* XS should be a read-level quantity */
    if (frame1->xs != *xs) {
      *xs = 0;
    } else if (frame2->xs != *xs) {
      *xs = 0;
    }
#endif

    return 16*(3-frame2->nti) + 4*(3-frame1->nti) + (3-frame0->nti);
  }
}


int
Readevid_cmp (const void *a, const void *b) {
  Readevid_T x = * (Readevid_T *) a;
  Readevid_T y = * (Readevid_T *) b;

  if (x->linei < y->linei) {
    return -1;
  } else if (y->linei < x->linei) {
    return +1;
  } else {
    return 0;
  }
}




/************************************************************************
 *   Tally_T
 ************************************************************************/

#define T Tally_T

T
Tally_new () {
  T new = (T) MALLOC(sizeof(*new));

  new->refnt = ' ';
  new->nmatches = 0;
  new->delcounts_plus = 0;
  new->delcounts_minus = 0;
  new->n_fromleft_plus = 0;
  new->n_fromleft_minus = 0;
  
#ifdef USE_MATCHPOOL
  new->matchpool = Matchpool_new();
#endif

  new->use_array_p = false;
  new->list_matches_byshift = (List_T) NULL;
  new->list_matches_bynm = (List_T) NULL;
  new->list_matches_byxs = (List_T) NULL;

  new->n_matches_byshift_plus = INITIAL_READLENGTH+1;
  new->n_matches_byshift_minus = INITIAL_READLENGTH+1;
  new->matches_byshift_plus = (int *) CALLOC(new->n_matches_byshift_plus,sizeof(int));
  new->matches_byshift_minus = (int *) CALLOC(new->n_matches_byshift_minus,sizeof(int));

  new->n_matches_bynm = INITIAL_READLENGTH+1;
  new->matches_bynm = (int *) CALLOC(new->n_matches_bynm,sizeof(int));

  new->n_matches_byxs = 3+1;	/* for 0, 1, and 2 */
  new->matches_byxs = (int *) CALLOC(new->n_matches_byxs,sizeof(int));

#ifdef USE_MISMATCHPOOL
  new->mismatchpool = Mismatchpool_new();
#endif

  new->mismatches_byshift = (List_T) NULL;
  new->mismatches_bynm = (List_T) NULL;
  new->mismatches_byxs = (List_T) NULL;

  new->insertions_byshift = (List_T) NULL;
  new->insertions_bynm = (List_T) NULL;

  new->deletions_byshift = (List_T) NULL;
  new->deletions_bynm = (List_T) NULL;

  new->readevidence = (List_T) NULL;

  return new;
}


void
Tally_clear (T this) {
  List_T ptr;
  Match_T match;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;
  Readevid_T readevid;


  this->refnt = ' ';
  this->nmatches = 0;
  this->delcounts_plus = 0;
  this->delcounts_minus = 0;
  this->n_fromleft_plus = 0;
  this->n_fromleft_minus = 0;

  if (this->use_array_p == true) {
#if 1
    /* Note: these memset instructions are necessary to get correct results */
    memset((void *) this->matches_byshift_plus,0,this->n_matches_byshift_plus * sizeof(int));
    memset((void *) this->matches_byshift_minus,0,this->n_matches_byshift_minus * sizeof(int));
    memset((void *) this->matches_bynm,0,this->n_matches_bynm * sizeof(int));
    memset((void *) this->matches_byxs,0,this->n_matches_byxs * sizeof(int));
#endif
    this->use_array_p = false;
  } else {

#ifdef USE_MATCHPOOL
    Matchpool_reset(this->matchpool);
#else
    for (ptr = this->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
      match = (Match_T) List_head(ptr);
      Match_free(&match);
    }
    List_free(&(this->list_matches_byshift));
    this->list_matches_byshift = (List_T) NULL;

    for (ptr = this->list_matches_bynm; ptr != NULL; ptr = List_next(ptr)) {
      match = (Match_T) List_head(ptr);
      Match_free(&match);
    }
    List_free(&(this->list_matches_bynm));
    this->list_matches_bynm = (List_T) NULL;

    for (ptr = this->list_matches_byxs; ptr != NULL; ptr = List_next(ptr)) {
      match = (Match_T) List_head(ptr);
      Match_free(&match);
    }
    List_free(&(this->list_matches_byxs));
    this->list_matches_byxs = (List_T) NULL;
#endif
  }

#ifdef USE_MISMATCHPOOL
  Mismatchpool_reset(this->mismatchpool);
#else
  for (ptr = this->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_byshift));
  this->mismatches_byshift = (List_T) NULL;

  for (ptr = this->mismatches_bynm; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_bynm));
  this->mismatches_bynm = (List_T) NULL;

  for (ptr = this->mismatches_byxs; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&(this->mismatches_byxs));
  this->mismatches_byxs = (List_T) NULL;
#endif


  for (ptr = this->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&(this->insertions_byshift));
  this->insertions_byshift = (List_T) NULL;

  for (ptr = this->insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&(this->insertions_bynm));
  this->insertions_bynm = (List_T) NULL;

  for (ptr = this->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&(this->deletions_byshift));
  this->deletions_byshift = (List_T) NULL;

  for (ptr = this->deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&(this->deletions_bynm));
  this->deletions_bynm = (List_T) NULL;


  for (ptr = this->readevidence; ptr != NULL; ptr = List_next(ptr)) {
    readevid = (Readevid_T) List_head(ptr);
    Readevid_free(&readevid);
  }
  List_free(&(this->readevidence));
  this->readevidence = (List_T) NULL;


  return;
}


void
Tally_transfer (T *dest, T *src) {
  T temp;


  temp = *dest;
  *dest = *src;

  temp->refnt = ' ';
  temp->nmatches = 0;
  temp->delcounts_plus = 0;
  temp->delcounts_minus = 0;
  temp->n_fromleft_plus = 0;
  temp->n_fromleft_minus = 0;

  if (temp->use_array_p == true) {
    memset((void *) temp->matches_byshift_plus,0,temp->n_matches_byshift_plus * sizeof(int));
    memset((void *) temp->matches_byshift_minus,0,temp->n_matches_byshift_minus * sizeof(int));
    memset((void *) temp->matches_bynm,0,temp->n_matches_bynm * sizeof(int));
    memset((void *) temp->matches_byxs,0,temp->n_matches_byxs * sizeof(int));
    temp->use_array_p = false;
  }
  temp->list_matches_byshift = (List_T) NULL;
  temp->list_matches_bynm = (List_T) NULL;
  temp->list_matches_byxs = (List_T) NULL;

  temp->mismatches_byshift = (List_T) NULL;
  temp->mismatches_bynm = (List_T) NULL;
  temp->mismatches_byxs = (List_T) NULL;

  temp->insertions_byshift = (List_T) NULL;
  temp->insertions_bynm = (List_T) NULL;
  temp->deletions_byshift = (List_T) NULL;
  temp->deletions_bynm = (List_T) NULL;

  temp->readevidence = (List_T) NULL;

  *src = temp;

  return;
}



void
Tally_free (T *old) {
  List_T ptr;
  Match_T match;
  Mismatch_T mismatch;
  Insertion_T ins;
  Deletion_T del;
  Readevid_T readevid;


#if 0
  (*old)->refnt = ' ';
  (*old)->nmatches = 0;
  (*old)->delcounts_plus = 0;
  (*old)->delcounts_minus = 0;
  (*old)->n_fromleft_plus = 0;
  (*old)->n_fromleft_minus = 0;
#endif

  FREE((*old)->matches_byshift_plus);
  FREE((*old)->matches_byshift_minus);
  FREE((*old)->matches_bynm);
  FREE((*old)->matches_byxs);

#ifdef USE_MATCHPOOL
  Matchpool_free(&((*old)->matchpool));
#else
  for (ptr = (*old)->list_matches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->list_matches_byshift));
  (*old)->list_matches_byshift = (List_T) NULL;

  for (ptr = (*old)->list_matches_bynm; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->list_matches_bynm));
  (*old)->list_matches_bynm = (List_T) NULL;

  for (ptr = (*old)->list_matches_byxs; ptr != NULL; ptr = List_next(ptr)) {
    match = (Match_T) List_head(ptr);
    Match_free(&match);
  }
  List_free(&((*old)->list_matches_byxs));
  (*old)->list_matches_byxs = (List_T) NULL;
#endif

#ifdef USE_MISMATCHPOOL
  Mismatchpool_free(&(*old)->mismatchpool);
#else
  for (ptr = (*old)->mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_byshift));
  (*old)->mismatches_byshift = (List_T) NULL;

  for (ptr = (*old)->mismatches_bynm; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_bynm));
  (*old)->mismatches_bynm = (List_T) NULL;

  for (ptr = (*old)->mismatches_byxs; ptr != NULL; ptr = List_next(ptr)) {
    mismatch = (Mismatch_T) List_head(ptr);
    Mismatch_free(&mismatch);
  }
  List_free(&((*old)->mismatches_byxs));
  (*old)->mismatches_byxs = (List_T) NULL;
#endif


  for (ptr = (*old)->insertions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&((*old)->insertions_byshift));
  /* (*old)->insertions_byshift = (List_T) NULL; */

  for (ptr = (*old)->insertions_bynm; ptr != NULL; ptr = List_next(ptr)) {
    ins = (Insertion_T) List_head(ptr);
    Insertion_free(&ins);
  }
  List_free(&((*old)->insertions_bynm));
  /* (*old)->insertions_bynm = (List_T) NULL; */

  for (ptr = (*old)->deletions_byshift; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&((*old)->deletions_byshift));
  /* (*old)->deletions_byshift = (List_T) NULL; */

  for (ptr = (*old)->deletions_bynm; ptr != NULL; ptr = List_next(ptr)) {
    del = (Deletion_T) List_head(ptr);
    Deletion_free(&del);
  }
  List_free(&((*old)->deletions_bynm));
  /* (*old)->deletions_bynm = (List_T) NULL; */


  for (ptr = (*old)->readevidence; ptr != NULL; ptr = List_next(ptr)) {
    readevid = (Readevid_T) List_head(ptr);
    Readevid_free(&readevid);
  }
  List_free(&((*old)->readevidence));
  /* (*old)->readevidence = (List_T) NULL; */


  FREE(*old);
  *old = (T) NULL;

  return;
}


char
Tally_codoni_plus (Tally_T tally0, Tally_T tally1, Tally_T tally2,
		   Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		   Genome_T genome, Genomicpos_T chroffset) {
  char nti0, nti1, nti2;
  char refnt;

  switch (tally0->refnt) {
  case 'A': nti0 = 0; break;
  case 'C': nti0 = 1; break;
  case 'G': nti0 = 2; break;
  case 'T': nti0 = 3; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos0-1U);
    switch (refnt) {
    case 'A': nti0 = 0; break;
    case 'C': nti0 = 1; break;
    case 'G': nti0 = 2; break;
    case 'T': nti0 = 3; break;
    default: return -1;
    }
  }

  switch (tally1->refnt) {
  case 'A': nti1 = 0; break;
  case 'C': nti1 = 1; break;
  case 'G': nti1 = 2; break;
  case 'T': nti1 = 3; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos1-1U);
    switch (refnt) {
    case 'A': nti1 = 0; break;
    case 'C': nti1 = 1; break;
    case 'G': nti1 = 2; break;
    case 'T': nti1 = 3; break;
    default: return -1;
    }
  }

  switch (tally2->refnt) {
  case 'A': nti2 = 0; break;
  case 'C': nti2 = 1; break;
  case 'G': nti2 = 2; break;
  case 'T': nti2 = 3; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos2-1U);
    switch (refnt) {
    case 'A': nti2 = 0; break;
    case 'C': nti2 = 1; break;
    case 'G': nti2 = 2; break;
    case 'T': nti2 = 3; break;
    default: return -1;
    }
  }

  return 16*nti0 + 4*nti1 + nti2;
}

char
Tally_codoni_minus (Tally_T tally0, Tally_T tally1, Tally_T tally2,
		    Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		    Genome_T genome, Genomicpos_T chroffset) {
  char nti0, nti1, nti2;
  char refnt;

  switch (tally0->refnt) {
  case 'A': nti0 = 3; break;
  case 'C': nti0 = 2; break;
  case 'G': nti0 = 1; break;
  case 'T': nti0 = 0; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos0-1U);
    switch (refnt) {
    case 'A': nti0 = 3; break;
    case 'C': nti0 = 2; break;
    case 'G': nti0 = 1; break;
    case 'T': nti0 = 0; break;
    default: return -1;
    }
  }

  switch (tally1->refnt) {
  case 'A': nti1 = 3; break;
  case 'C': nti1 = 2; break;
  case 'G': nti1 = 1; break;
  case 'T': nti1 = 0; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos1-1U);
    switch (refnt) {
    case 'A': nti1 = 3; break;
    case 'C': nti1 = 2; break;
    case 'G': nti1 = 1; break;
    case 'T': nti1 = 0; break;
    default: return -1;
    }
  }

  switch (tally2->refnt) {
  case 'A': nti2 = 3; break;
  case 'C': nti2 = 2; break;
  case 'G': nti2 = 1; break;
  case 'T': nti2 = 0; break;
  default:
    refnt = Genome_get_char(genome,chroffset+chrpos2-1U);
    switch (refnt) {
    case 'A': nti2 = 3; break;
    case 'C': nti2 = 2; break;
    case 'G': nti2 = 1; break;
    case 'T': nti2 = 0; break;
    default: return -1;
    }
  }

  return 16*nti2 + 4*nti1 + nti0;
}
