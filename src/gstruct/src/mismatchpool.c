static char rcsid[] = "$Id: mismatchpool.c 135603 2014-05-08 20:10:58Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mismatchpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "mem.h"
#include "mismatchdef.h"
#include "listdef.h"


#define CHUNKSIZE 1000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For mechanics of memory allocation and deallocation */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* For popping */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Mismatchpool_T
struct T {
  int nobjects;
  int objectctr;
  struct Mismatch_T *objectptr;
  List_T objectchunks;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;
};

void
Mismatchpool_free (T *old) {
  List_T p;
  struct Mismatch_T *objectptr;
  struct List_T *listcellptr;

  if (*old) {
    for (p = (*old)->objectchunks; p != NULL; p = List_next(p)) {
      objectptr = (struct Mismatch_T *) List_head(p);
      FREE(objectptr);
    }
    List_free(&(*old)->objectchunks);
    for (p = (*old)->listcellchunks; p != NULL; p = List_next(p)) {
      listcellptr = (struct List_T *) List_head(p);
      FREE(listcellptr);
    }
    List_free(&(*old)->listcellchunks);
    FREE(*old);
  }
  return;
}

void
Mismatchpool_free_memory (T this) {
  List_T p;
  struct Mismatch_T *objectptr;
  struct List_T *listcellptr;

  for (p = this->objectchunks; p != NULL; p = List_next(p)) {
    objectptr = (struct Mismatch_T *) List_head(p);
    FREE(objectptr);
  }
  List_free(&this->objectchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE(listcellptr);
  }
  List_free(&this->listcellchunks);

  this->nobjects = 0;
  this->objectctr = 0;
  this->objectchunks = NULL;
  /* this->objectptr = add_new_objectchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}


void
Mismatchpool_report_memory (T this) {
  printf("Mismatchpool has %d pairchunks and %d listcellchunks\n",
	 List_length(this->objectchunks),List_length(this->listcellchunks));
  return;
}


static struct Mismatch_T *
add_new_objectchunk (T this) {
  struct Mismatch_T *chunk;

  chunk = (struct Mismatch_T *) MALLOC(CHUNKSIZE*sizeof(struct Mismatch_T));
  this->objectchunks = List_push_keep(this->objectchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of objects.  Ptr for object %d is %p\n",
		this->nobjects,chunk));

  this->nobjects += CHUNKSIZE;

  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push_keep(this->listcellchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;

  return chunk;
}

T
Mismatchpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->nobjects = 0;
  new->objectctr = 0;
  new->objectchunks = NULL;
  /* new->objectptr = add_new_objectchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Mismatchpool_reset (T this) {
  this->objectctr = 0;
  this->listcellctr = 0;
  return;
}

List_T
Mismatchpool_push (List_T list, T this, char nt, int shift, int mapq, char quality, int xs) {
  List_T listcell;
  Mismatch_T new;
  List_T p;
  int n;

  if (this->objectctr >= this->nobjects) {
    this->objectptr = add_new_objectchunk(this);
  } else if ((this->objectctr % CHUNKSIZE) == 0) {
    for (n = this->nobjects - CHUNKSIZE, p = this->objectchunks;
	 n > this->objectctr; p = p->rest, n -= CHUNKSIZE) ;
    this->objectptr = (struct Mismatch_T *) p->first;
    debug1(printf("Located object %d at %p\n",this->objectctr,this->objectptr));
  }    
  new = this->objectptr++;
  this->objectctr++;


  new->nt = nt;
  new->shift = shift;
  new->mapq = mapq;
  new->quality = quality;
  new->xs = xs;
  new->count = 1;

  /* Assigned when mismatch added to unique list */
  /* new->count_plus = 0; */
  /* new->count_minus = 0; */
  
  new->next = NULL;

  debug(printf("Creating %p: %d %d %c\n",new,shift,mapq,quality));


  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) new;
  listcell->rest = list;

  return listcell;
}



/* Note: this does not free the list cell */
List_T
Mismatchpool_pop (List_T list, Mismatch_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Mismatch_T) list->first;
    return head;
  } else {
    return list;
  }
}


