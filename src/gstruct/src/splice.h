/* $Id: splice.h 107778 2013-09-12 20:12:22Z twu $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "list.h"
#include "iit-read.h"


typedef enum {VALID, MINCOUNT, MINSUPPORT, BADPROB, NONCANONICAL, WRONGDIR_CMP, WRONGDIR_EXTENT} Splicestatus_T;
typedef enum {NOVEL, KNOWN_MAIN, KNOWN_ALT} Spliceknownp_T;

extern char *
Splicestatus_string (Splicestatus_T status);


#define T Splice_T
typedef struct T *T;

extern void
Splice_free (T *old);

extern void
Splice_gc (List_T *splices);

extern bool
Splice_intron_usedp (T this);

extern void
Splice_set_intron_usedp (T this);

extern int
Splice_count (T this);

extern int
Splice_known_count (T this);

extern int
Splice_primary_count (T this);

extern int
Splice_crosshyb_count (T this);

extern int
Splice_maxminsupport (T this);

extern int
Splice_npositions (T this, Genomicpos_T genestart, Genomicpos_T geneend, int insertlength, int readlength,
		   int min_overhang);

extern double
Splice_density (T this, Genomicpos_T genestart, Genomicpos_T geneend, int insertlength, int readlength,
		int min_overhang);

extern int
Splice_sign (T this);

extern Genomicpos_T
Splice_donorpos (T this);

extern Genomicpos_T
Splice_acceptorpos (T this);

extern Genomicpos_T
Splice_low (T this);

extern Genomicpos_T
Splice_high (T this);

extern bool
Splice_canonicalp (T this);

extern double
Splice_donorprob (T this);

extern double
Splice_acceptorprob (T this);

extern Splicestatus_T
Splice_status (T this);

extern Spliceknownp_T
Splice_knownp (T this);

extern bool
Splice_validp (T this);

extern bool
Splice_valid_end_p (T this);

extern int
Splice_level (T this);

extern T *
Splice_array_copy (T *array, int nsplices);

extern List_T
Splice_valid_list (List_T splices);

extern int
Splice_cmp (const void *a, const void *b);

extern int
Splice_low_cmp (const void *a, const void *b);

extern int
Splice_high_cmp (const void *a, const void *b);

extern int
Splice_bycount_cmp (const void *a, const void *b);

extern int
Splice_compute_levels (List_T splices, Genomicpos_T mincoord,
		       Genomicpos_T maxcoord, int max_allowed_levels,
		       double xfactor, int mincount);

extern T
Splice_new_known (Genomicpos_T donorpos, Genomicpos_T acceptorpos, int sign,
		  bool mainpath_p);

extern void
Splice_transfer_info (T knownsplice, T obssplice);

extern List_T
Splice_add_at_donor (List_T list, int sign, Genomicpos_T donorpos, Genomicpos_T acceptorpos,
		     bool crosshybp, int nhits, int support1, int support2,
		     int overhang1, int overhang2, int querylength, bool low_end_p, bool canonicalp,
		     double donorprob, double acceptorprob,
		     char donor1, char donor2, char acceptor1, char acceptor2, bool knownp);

extern void
Splice_filter_list (List_T splices, int mincount_end_alt, int minsupport, bool need_canonical_p);

extern void
Splice_filter_array (T *array, int nsplices, int mincount_end_alt, int minsupport, bool need_canonical_p);

extern void
Splice_resolve (T *array, int nsplices, long int *fwd_runlengths, long int *rev_runlengths,
		Genomicpos_T chrlength, Genomicpos_T max_exonlength, char *chr);

extern void
Splice_turn (T *array, int nsplices, Genomicpos_T max_exonlength, char *chr);

extern void
Splice_print (T splice, char *chr, IIT_T genes_iit, bool show_invalid_p);

#undef T

#endif

