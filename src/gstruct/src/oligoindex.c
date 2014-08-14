#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bool.h"
#include "mem.h"
#include "oligoindex.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



typedef unsigned int Shortoligomer_T;
typedef int Chrpos_T;

#define T Oligoindex_T
struct T {
  Chrpos_T **pointers;
  Chrpos_T **positions;
  int *counts;
};

static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static char *
shortoligo_nt (Shortoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Shortoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}

static int
allocate_positions (Chrpos_T **pointers, Chrpos_T **positions,
		    int *counts, int oligospace, int indexsize,
		    Shortoligomer_T mask, char *sequence, int seqlength) {
  int i = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
  Chrpos_T *ptr;
  int totalcounts;
  int sequencepos = 0;

  sequencepos -= indexsize;
  for (i = 0, p = sequence; i < seqlength; i++, p++) {
    in_counter++;
    sequencepos++;

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0;
    }

    if (in_counter == indexsize) {
      masked = oligo & mask;
      counts[masked] += 1;
      in_counter--;
    }
  }

  totalcounts = 0;
  for (i = 0; i < oligospace; i++) {
    totalcounts += counts[i];
  }

  if (totalcounts == 0) {
    positions[0] = (Chrpos_T *) NULL;
  } else {
    ptr = (Chrpos_T *) CALLOC(totalcounts,sizeof(Chrpos_T));
    for (i = 0; i < oligospace; i++) {
      positions[i] = ptr;
      ptr += counts[i];
    }
    positions[i] = ptr;		/* For positions[oligospace], used for indicating if pointer hits next position */
    /* Does not copy positions[oligospace] */
    memcpy((void *) pointers,positions,oligospace*sizeof(Chrpos_T *));
  }

  return totalcounts;
}


static int
store_positions (Chrpos_T **pointers, int indexsize,
		 Shortoligomer_T mask, char *sequence, int seqlength) {
  int nstored = 0;
  int i = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
  int sequencepos = 0;

  sequencepos -= indexsize;
  for (i = 0, p = sequence; i < seqlength; i++, p++) {
    in_counter++;
    sequencepos++;

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0;
    }

    if (in_counter == indexsize) {
      masked = oligo & mask;
      pointers[masked][0] = (Chrpos_T) sequencepos;
      pointers[masked]++;
      nstored++;
      in_counter--;
    }
  }

  return nstored;
}


static T
Oligoindex_new (char *sequence, int seqlength) {
  T new = MALLOC(sizeof(*new));
  int nallocated, nstored;
  int indexsize = 8;
  int oligospace;

  oligospace = power(4,indexsize);
  new->positions = (Chrpos_T **) CALLOC(oligospace+1,sizeof(Chrpos_T *));
  new->pointers = (Chrpos_T **) CALLOC(oligospace,sizeof(Chrpos_T *));
  new->counts = (int *) CALLOC(oligospace,sizeof(int));

  if ((nallocated = allocate_positions(new->pointers,new->positions,new->counts,oligospace,
				       /*indexsize*/8,/*mask*/STRAIGHT_MASK_8,
				       sequence,seqlength)) > 0) {
    nstored = store_positions(new->pointers,/*indexsize*/8,/*mask*/STRAIGHT_MASK_8,
			      sequence,seqlength);

    if (nstored != nallocated) {
      fprintf(stderr,"Bug in create_tally: %d allocated, but %d stored\n",
	      nallocated,nstored);
      abort();
    }
  }

  return new;
}


static void
Oligoindex_free (T *old) {
  FREE((*old)->pointers);
  FREE((*old)->positions[0]);
  FREE((*old)->positions);
  FREE((*old)->counts);
  FREE(*old);

  return;
}


static int
compute_nmatches (int *maxpos, char *x, int posi, int xlen, char *y, int posj, int ylen,
		  int indexsize, int width) {
  int max_nmatches = 0, nmatches = 0, i, j, k;
  int startposi, startposj, endposi, endposj;
  
  if (posi < posj) {
    if (posi - width < 0) {
      startposi = 0;
      startposj = posj - posi;
    } else {
      startposi = posi - width;
      startposj = posj - width;
    }
  } else {
    if (posj - width < 0) {
      startposi = posi - posj;
      startposj = 0;
    } else {
      startposi = posi - width;
      startposj = posj - width;
    }
  }
    
  if (xlen - posi < ylen - posj) {
    if (posi + indexsize + width >= xlen) {
      endposi = xlen;
      endposj = posj - posi + xlen;
    } else {
      endposi = posi + indexsize + width;
      endposj = posj + indexsize + width;
    }
  } else {
    if (posj + indexsize + width >= ylen) {
      endposi = posi - posj + ylen;
      endposj = ylen;
    } else {
      endposi = posi + indexsize + width;
      endposj = posj + indexsize + width;
    }
  }

  debug(printf("startposi %d, startposj %d, endposi %d, endposj %d\n",
	       startposi,startposj,endposi,endposj));

  for (k = 0, i = startposi, j = startposj;
       k < width && i < endposi && j < endposj;
       k++, i++, j++) {
    if (x[i] == y[j]) {
      nmatches++;
    }
  }

  if (k < width) {
    *maxpos = -1;
    return 0;
  } else {
    *maxpos = startposi - posi;
    max_nmatches = nmatches;

    while (i < endposi && j < endposj) {
      if (x[i] == y[j]) {
	nmatches++;
      }
      if (x[i-width] == y[j-width]) {
	nmatches--;
      }
      if (nmatches > max_nmatches) {
	*maxpos = (i - width + 1) - posi;
	max_nmatches = nmatches;
      }
      i++;
      j++;
    }

    debug(printf("max_nmatches %d at maxpos %d\n",max_nmatches,*maxpos));

    return max_nmatches;
  }
}


int
Oligoindex_compare_seqs (int *maxposi, int *maxposj, char *donor_sequence, char *acceptor_sequence,
			 int width) {
  T this;
  int max_nmatches = 0, nmatches;
  int donor_len = strlen(donor_sequence);
  int acceptor_len = strlen(acceptor_sequence);
  int maxpos;
  int nhits, i, k;
  int indexsize = 8;
  Shortoligomer_T mask = STRAIGHT_MASK_8;

  Shortoligomer_T oligo = 0U, masked;
  char *p;
  int in_counter = 0;
  int sequencepos = 0;

  debug(printf("Entered Oligoindex_compare_seqs\n"));

  sequencepos -= indexsize;
  if (donor_len > acceptor_len) {
    debug(printf("Donor %d is longer than acceptor %d\n",donor_len,acceptor_len));
    this = Oligoindex_new(donor_sequence,donor_len);
    for (i = 0, p = acceptor_sequence; i < acceptor_len; i++, p++) {
      in_counter++;
      sequencepos++;

      switch (*p) {
      case 'A': oligo = (oligo << 2); break;
      case 'C': oligo = (oligo << 2) | 1; break;
      case 'G': oligo = (oligo << 2) | 2; break;
      case 'T': oligo = (oligo << 2) | 3; break;
      default: oligo = 0U; in_counter = 0;
      }

      if (in_counter == indexsize) {
	masked = oligo & mask;
	nhits = this->counts[masked];
	debug(printf("acceptor pos %d: oligo %u => %d counts\n",sequencepos,masked,nhits));
	for (k = 0; k < nhits; k++) {
	  /* Search in donor_sequence around this->positions[masked]; */
	  debug(printf("Found match between donor pos %d and acceptor pos %d\n",
		       this->positions[masked][k],sequencepos));
	  if ((nmatches = compute_nmatches(&maxpos,donor_sequence,/*donorpos*/this->positions[masked][k],donor_len,
					   acceptor_sequence,/*acceptorpos*/sequencepos,acceptor_len,
					   indexsize,width)) > max_nmatches) {
	    *maxposi = this->positions[masked][k] + maxpos;
	    *maxposj = sequencepos + maxpos;
	    max_nmatches = nmatches;
	  }
	}
	in_counter--;
      }
    }

  } else {
    debug(printf("Acceptor %d is longer than donor %d\n",acceptor_len,donor_len));
    this = Oligoindex_new(acceptor_sequence,acceptor_len);
    for (i = 0, p = donor_sequence; i < donor_len; i++, p++) {
      in_counter++;
      sequencepos++;

      switch (*p) {
      case 'A': oligo = (oligo << 2); break;
      case 'C': oligo = (oligo << 2) | 1; break;
      case 'G': oligo = (oligo << 2) | 2; break;
      case 'T': oligo = (oligo << 2) | 3; break;
      default: oligo = 0U; in_counter = 0;
      }

      if (in_counter == indexsize) {
	masked = oligo & mask;
	nhits = this->counts[masked];
	debug(printf("donor pos %d: oligo %u => %d counts\n",sequencepos,masked,nhits));
	for (k = 0; k < nhits; k++) {
	  /* Search in acceptor_sequence around this->positions[masked]; */
	  debug(printf("Found match between donor pos %d and acceptor pos %d\n",
		       sequencepos,this->positions[masked][k]));
	  if ((nmatches = compute_nmatches(&maxpos,donor_sequence,/*donorpos*/sequencepos,donor_len,
					   acceptor_sequence,/*acceptorpos*/this->positions[masked][k],
					   acceptor_len,indexsize,width)) > max_nmatches) {
	    *maxposi = sequencepos + maxpos;
	    *maxposj = this->positions[masked][k] + maxpos;
	    max_nmatches = nmatches;
	  }
	}
	in_counter--;
      }
    }
  }

  Oligoindex_free(&this);

  return max_nmatches;
}

