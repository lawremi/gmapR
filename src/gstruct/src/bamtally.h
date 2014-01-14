/* $Id: bamtally.h 123358 2014-01-14 00:34:54Z twu $ */
#ifndef BAMTALLY_INCLUDED
#define BAMTALLY_INCLUDED
#include "bamread.h"
#include "genome.h"
#include "genomicpos.h"
#include "list.h"
#include "iit-read.h"
#include "tableuint.h"


typedef enum {OUTPUT_BLOCKS, OUTPUT_RUNLENGTHS, OUTPUT_TALLY, OUTPUT_IIT, OUTPUT_TOTAL} Tally_outputtype_T;
#define DEFAULT_QUALITY 40  /* quality_score_adj + 40 */

extern long int
Bamtally_run (long int **tally_matches, long int **tally_mismatches,
	      List_T *intervallist, List_T *labellist, List_T *datalist,
	      int *quality_counts_match, int *quality_counts_mismatch,
	      Bamreader_T bamreader, Genome_T genome, char *printchr,
	      Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, int alloclength,
	      Tableuint_T resolve_low_table, Tableuint_T resolve_high_table,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      bool ignore_lowend_p, bool ignore_highend_p,
	      Tally_outputtype_T output_type, bool blockp, int blocksize,
	      int quality_score_adj, int min_depth, int variant_strands,
	      bool genomic_diff_p, bool signed_counts_p, bool ignore_query_Ns_p,
	      bool print_indels_p, bool print_totals_p, bool print_cycles_p,
	      bool print_quality_scores_p, bool print_mapq_scores_p, bool want_genotypes_p,
	      bool verbosep, bool readlevel_p);

/* Version that keeps separate tallies for the low and high ends of paired-end reads */
extern void
Bamtally_run_lh (long int *tally_matches_low, long int *tally_mismatches_low,
		 long int *tally_matches_high, long int *tally_mismatches_high,
		 int *quality_counts_match, int *quality_counts_mismatch,
		 Bamreader_T bamreader, Genome_T genome, char *printchr,
		 Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend, int alloclength,
		 char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		 bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		 int blocksize, int quality_score_adj, int min_depth, int variant_strands,
		 bool genomic_diff_p, bool ignore_query_Ns_p, bool verbosep, bool readlevel_p);

extern IIT_T
Bamtally_iit (Bamreader_T bamreader, char *desired_chr, char *bam_lacks_chr,
	      Genomicpos_T chrstart, Genomicpos_T chrend,
	      Genome_T genome, IIT_T chromosome_iit, int alloclength,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      int min_depth, int variant_strands, bool ignore_query_Ns_p,
	      bool print_indels_p, int blocksize, bool verbosep, bool readlevel_p);

#endif
