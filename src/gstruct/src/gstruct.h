/* $Id: gstruct.h 127943 2014-02-19 23:49:02Z twu $ */
#ifndef GSTRUCT_INCLUDED
#define GSTRUCT_INCLUDED
#include "bool.h"
#include "table.h"
#include "tableint.h"
#include "genomicpos.h"
#include "genome.h"
#include "iit-read.h"
#include "uinttable.h"
#include "splice.h"
#include "bamread.h"


#ifdef BAM_INPUT
#include "bamread.h"
#endif

#define T Gstruct_T
typedef struct T *T;


extern void
Gstruct_free (T *old);

extern int
Gstruct_npairs (T this);


#ifdef BAM_INPUT

extern T
Gstruct_bam_input (int *mean_readlength, int *mean_insertlength, List_T bamfiles,
		   char *chromosome, Genomicpos_T chrstart, Genomicpos_T chrend,
		   Tableint_T ngoodhits_low_table, Tableint_T ngoodhits_high_table,
		   IIT_T genes_iit, Genomicpos_T shortsplicedist, Genomicpos_T max_pairlength,
		   Genome_T genome, IIT_T chromosome_iit,
		   char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		   bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		   bool trust_sam_p, bool need_canonical_p, char *bam_lacks_chr, int bam_lacks_chr_length);

extern T
Gstruct_knowngenes_input (IIT_T *introns_iit, IIT_T *middle_exons_iit, IIT_T *end_exons_iit,
			  IIT_T knowngenes_iit, Genome_T genome, IIT_T chromosome_iit);


extern void
Gstruct_recount_multimappers (Tableint_T *ngoodhits_low_table, Tableint_T *ngoodhits_high_table, Bamreader_T bamreader,
			      IIT_T genes_iit, int maximum_nhits);

#else

extern T
Gstruct_sam_input (int minimum_mapq, Genomicpos_T max_pairlength,
		   Genome_T genome, IIT_T chromosome_iit, Genomicpos_T shortsplicedist,
		   bool need_unique_p, bool need_primary_p, bool trust_sam_p, bool need_canonical_p);

#endif

extern char *
Gstruct_next_chr (List_T *splices, long int **fwd_extents, long int **rev_extents, long int **null_extents,
		  long int **primary_extents, long int **crosshyb_extents,
		  Genomicpos_T *chrlength, T this,
		  int max_exonlength, int mincount_end_alt, int minsupport, bool need_canonical_p,
		  char *bam_lacks_chr, int bam_lacks_chr_length);

extern List_T
Gstruct_get_known (T this, char *chr);


#undef T

#endif

