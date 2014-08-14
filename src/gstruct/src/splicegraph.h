/* $Id: splicegraph.h 58843 2012-03-02 01:04:39Z twu $ */
#ifndef SPLICEGRAPH_INCLUDED
#define SPLICEGRAPH_INCLUDED
#include "bool.h"
#include "list.h"
#include "genome.h"
#include "genomicpos.h"


typedef enum {OUTPUT_SPLICES, OUTPUT_PATHS, OUTPUT_GENES, ACCUMULATE_GENES, OUTPUT_DIFFS} Gene_outputtype_T;

typedef struct Gene_T *Gene_T;

extern void
Gene_free (Gene_T *old);

extern void
Splicegraph_genes_gc (List_T *genes);

#define T Splicegraph_T
typedef struct T *T;

extern void
Splicegraph_free (T *old);

extern T
Splicegraph_new ();

extern List_T
Splicegraph_solve_obs (T this, List_T obs_splices, List_T knowngenes, char *chr,
		       long int *tally_matches_low, long int *tally_mismatches_low,
		       long int *tally_matches_high, long int *tally_mismatches_high,
		       long int *primary_extents, long int *crosshyb_extents,
		       Genome_T genome, Genomicpos_T chroffset, Genomicpos_T chrlength,
		       int insertlength, int readlength, int min_overhang,
		       int mincount_alt, int minsupport, bool need_canonical_p,
		       int min_exonlength, int max_exonlength, int auto_exonlength,
		       int min_intronlength, int max_intronlength, bool altpaths_p);

extern List_T
Splicegraph_solve_known (T this, List_T splices, IIT_T introns_iit,
			 IIT_T middle_exons_iit, IIT_T end_exons_iit,
			 char *chr, Genomicpos_T chrlength);

extern void
Splicegraph_print_genes (FILE *fp, List_T genes, int npairs, IIT_T knowngenes_iit, char *chr,
			 long int *tally_matches_low, long int *tally_mismatches_low,
			 long int *tally_matches_high, long int *tally_mismatches_high,
			 int insertlength, int readlength, int min_overhang);

extern List_T
Splicegraph_splices (List_T genes);

#undef T
#endif

