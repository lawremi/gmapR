#ifndef TALLY_INCLUDED
#define TALLY_INCLUDED
#include "list.h"
#include "intlist.h"
#include "bool.h"

#include "genomicpos.h"
#include "genome.h"

#include "matchdef.h"
#include "matchpool.h"
#include "mismatchdef.h"
#include "mismatchpool.h"


#define MAX_QUALITY_SCORE 40	/* Also in bamtally.c */
#define MAX_MAPQ_SCORE 40	/* Also in bamtally.c */


extern Match_T
Match_new (int shift, int nm, int xs, int ncounts);
extern void
Match_free (Match_T *old);
extern Mismatch_T
Mismatch_new (char nt, int shift, int nm, int xs, int ncounts);
extern void
Mismatch_free (Mismatch_T *old);

typedef struct Softclip_jcn_T *Softclip_jcn_T;
struct Softclip_jcn_T {
  Genomicpos_T chrpos;
  int shift;			/* Used to record shifts */
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */
};


Softclip_jcn_T
Softclip_jcn_new (Genomicpos_T chrpos, int shift, int ncounts);
extern void
Softclip_jcn_free (Softclip_jcn_T *old);
extern int
Softclip_jcn_count_cmp (const void *a, const void *b);



typedef struct Insertion_T *Insertion_T;
struct Insertion_T {
  Genomicpos_T chrpos;
  char *segment;
  int mlength;
  int shift;			/* Used to record shifts */
  int nm;
  int xs;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Insertion_T next;		/* Used for linking similar insertions together */
};

extern Insertion_T
Insertion_new (Genomicpos_T chrpos, char *query_insert, int mlength, int shift, int nm, int xs, int ncounts);
extern void
Insertion_free (Insertion_T *old);
extern int
Insertion_count_cmp (const void *a, const void *b);
extern Insertion_T
find_insertion_byshift (List_T insertions, char *segment, int mlength, int shift);
extern Insertion_T
find_insertion_bynm (List_T insertions, char *segment, int mlength, int nm);
extern Insertion_T
find_insertion_byxs (List_T insertions, char *segment, int mlength, int xs);
extern Insertion_T
find_insertion_seg (List_T insertions, char *segment, int mlength);


typedef struct Deletion_T *Deletion_T;
struct Deletion_T {
  Genomicpos_T chrpos;
  char *segment;
  int mlength;
  int shift;			/* Used to record shifts */
  int nm;
  int xs;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Deletion_T next;		/* Used for linking similar deletions together */
};

extern Deletion_T
Deletion_new (Genomicpos_T chrpos, char *deletion, int mlength, int shift, int nm, int xs, int ncounts);
extern void
Deletion_free (Deletion_T *old);
extern int
Deletion_count_cmp (const void *a, const void *b);
extern Deletion_T
find_deletion_byshift (List_T deletions, char *segment, int mlength, int shift);
extern Deletion_T
find_deletion_bynm (List_T deletions, char *segment, int mlength, int nm);
extern Deletion_T
find_deletion_byxs (List_T deletions, char *segment, int mlength, int xs);
extern Deletion_T
find_deletion_seg (List_T deletions, char *segment, int mlength);


typedef struct Microinv_T *Microinv_T;
struct Microinv_T {
  Genomicpos_T chrpos;
  char *segment;
  int mlength;
  int shift;			/* Used to record shifts */
  int nm;
  int xs;
  long int count;

  long int count_plus;		/* Used by unique elements */
  long int count_minus;		/* Used by unique elements */

  Microinv_T next;		/* Used for linking similar microinversions together */
};


extern Microinv_T
Microinv_new (Genomicpos_T chrpos, char *microinv, int mlength, int shift, int nm, int xs, int ncounts);
extern void
Microinv_free (Microinv_T *old);
extern int
Microinv_count_cmp (const void *a, const void *b);
extern Microinv_T
find_microinv_byshift (List_T microinvs, char *segment, int mlength, int shift);
extern Microinv_T
find_microinv_bynm (List_T microinvs, char *segment, int mlength, int nm);
extern Microinv_T
find_microinv_byxs (List_T microinvs, char *segment, int mlength, int xs);
extern Microinv_T
find_microinv_seg (List_T microinvs, char *segment, int mlength);



typedef struct Readevid_T *Readevid_T;
extern Readevid_T
Readevid_new (unsigned int linei, int nreps, char nt, int shift, int nm, int xs);
extern unsigned int
Readevid_linei (Readevid_T this);
extern int
Readevid_nreps (Readevid_T this);
extern char
Readevid_nt (Readevid_T this);
extern int
Readevid_quality_score (Readevid_T this);

/* Needs to be signed char, because can return -1 for non-ACGT */
extern char
Readevid_codoni_plus (int *shift, int *nm, int *xs,
		      Readevid_T frame0, Readevid_T frame1, Readevid_T frame2);
/* Needs to be signed char, because can return -1 for non-ACGT */
extern char
Readevid_codoni_minus (int *shift, int *nm, int *xs,
		       Readevid_T frame0, Readevid_T frame1, Readevid_T frame2);
extern int
Readevid_cmp (const void *a, const void *b);



typedef struct Tally_T *Tally_T;
struct Tally_T {
  char refnt;
  int depth_passing;
  int depth_total;

  int nmatches_passing;
  int nmismatches_passing;
  int nsoftclip_nts_low_passing;
  int nsoftclip_nts_high_passing;

  int nmatches_total;		/* Kept for informational purposes */
  int nmismatches_total;	/* Kept for informational purposes */
  int nsoftclip_nts_low_total;	/* Kept for informational purposes */
  int nsoftclip_nts_high_total;	/* Kept for informational purposes */

  int delcounts_plus;
  int delcounts_minus;

#if 0
  long int n_fromleft_plus; /* Used for reference count for insertions */
  long int n_fromleft_minus; /* Used for reference count for insertions */
#endif

#ifdef USE_MATCHPOOL
  Matchpool_T matchpool;
#endif
#ifdef USE_MISMATCHPOOL
  Mismatchpool_T mismatchpool;
#endif

  bool use_array_p;
  List_T list_matches_byshift;
  List_T list_matches_bynm;
  List_T list_matches_byxs;

#ifdef REUSE_ARRAYS
  int avail_matches_byshift_plus;
  int avail_matches_byshift_minus;
#endif

  /* Includes only passing quality scores */
  int max_byshift_plus;
  int *matches_byshift_plus;
  int max_byshift_minus;
  int *matches_byshift_minus;

  int max_nm;
  int *matches_bynm;

  int *matches_byxs;

  /* Includes only passing quality scores */
  List_T mismatches_byshift;
  List_T mismatches_bynm;
  List_T mismatches_byxs;

  List_T softclip_nts_low_byshift;
  List_T softclip_nts_low_bynm;
  List_T softclip_nts_low_byxs;

  List_T softclip_nts_high_byshift;
  List_T softclip_nts_high_bynm;
  List_T softclip_nts_high_byxs;

  List_T softclip_jcns_low_byshift;
  List_T softclip_jcns_high_byshift;

  List_T insertions_byshift;
  List_T insertions_bynm;
  List_T insertions_byxs;

  List_T deletions_byshift;
  List_T deletions_bynm;
  List_T deletions_byxs;

  List_T microinvs_byshift;
  List_T microinvs_bynm;
  List_T microinvs_byxs;

  List_T readevidence;

  List_T bamlines;		/* Stored at chrpos_high */
  Intlist_T bamline_nreps_plus;
  Intlist_T bamline_nreps_minus;
};


extern Tally_T
Tally_new ();
extern void
Tally_clear (Tally_T this);
extern void
Tally_transfer (Tally_T *dest, Tally_T *src);
extern void
Tally_free (Tally_T *old);
extern char
Tally_codoni_plus (Tally_T tally0, Tally_T tally1, Tally_T tally2,
		   Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		   Genome_T genome, Genomicpos_T chroffset);
extern char
Tally_codoni_minus (Tally_T tally0, Tally_T tally1, Tally_T tally2,
		    Genomicpos_T chrpos0, Genomicpos_T chrpos1, Genomicpos_T chrpos2,
		    Genome_T genome, Genomicpos_T chroffset);


#endif


