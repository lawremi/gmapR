/* $Id: bamread.h 178960 2015-11-16 19:52:26Z twu $ */
#ifndef BAMREAD_INCLUDED
#define BAMREAD_INCLUDED
/* Cannot use bool, since it appears to conflict with samtools */
#include <stdio.h>
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "genome.h"
#include "iit-read.h"

#include "table.h"
#include "uinttable.h"
#include "chrom.h"
#include "list.h"


#define T Bamreader_T
typedef struct T *T;

extern void
Bamread_free (T *old);

extern T
Bamread_new (char *filename);

extern void
Bamread_write_header (T this);

extern Genomicpos_T
Bamread_chrlength (T this, char *chr);

extern bool
Bamread_limit_region (T this, char *chr, Genomicpos_T chrstart, Genomicpos_T chrend);

extern void
Bamread_unlimit_region (T this);

extern int
Bamread_nreads (int *npositions, T this, char *chr, Genomicpos_T chrpos1, Genomicpos_T chrpos2);

extern int
Bamread_next_line (T this, char **acc, unsigned int *flag, int *mapq, char **chr, Genomicpos_T *chrpos,
		   char **mate_chr, Genomicpos_T *mate_chrpos,
		   Intlist_T *cigartypes, Uintlist_T *cigarlengths, int *cigarlength,
		   int *readlength, char **read, char **quality_string, char **hardclip, char **hardclip_quality,
		   char **read_group, bool *terminalp);

typedef struct Bamline_T *Bamline_T;

extern char *
Bamline_acc (Bamline_T this);
extern unsigned int
Bamline_flag (Bamline_T this);
extern int
Bamline_concordantp (Bamline_T this);
extern int
Bamline_lowend_p (Bamline_T this);
extern int
Bamline_paired_read_p (Bamline_T this);
extern int
Bamline_firstend_p (Bamline_T this);
extern int
Bamline_hiti (Bamline_T this);
extern int
Bamline_nhits (Bamline_T this);
extern bool
Bamline_good_unique_p (Bamline_T this);
extern bool
Bamline_perfect_match_p (Bamline_T this);
extern int
Bamline_mapq (Bamline_T this);
extern char *
Bamline_chr (Bamline_T this);
extern Genomicpos_T
Bamline_chrpos_low (Bamline_T this);
extern Genomicpos_T
Bamline_chrpos_low_noclip (Bamline_T this);
extern char *
Bamline_mate_chr (Bamline_T this);
extern Genomicpos_T
Bamline_mate_chrpos_low (Bamline_T this);
extern int
Bamline_insert_length (Bamline_T this);
extern Intlist_T
Bamline_cigar_types (Bamline_T this);
extern Uintlist_T
Bamline_cigar_npositions (Bamline_T this);
extern Intlist_T
Bamline_diffcigar (int *min_overhang, Uintlist_T *npositions, Uintlist_T *chrpositions, Bamline_T this);
extern int
Bamline_cigar_querylength (Bamline_T this);
extern void
Bamread_print_cigar (FILE *fp, Bamline_T this);
extern char *
Bamline_cigar_string (Bamline_T this);
extern int
Bamline_readlength (Bamline_T this);
extern char *
Bamline_read (Bamline_T this);
extern char *
Bamline_quality_string (Bamline_T this);
extern char *
Bamline_hardclip (Bamline_T this);
extern char *
Bamline_hardclip_quality (Bamline_T this);
extern bool
Bamline_terminalp (Bamline_T this);
extern char *
Bamline_read_group (Bamline_T this);
extern void
Bamline_print (FILE *fp, Bamline_T this, unsigned int newflag, int quality_score_adj);
extern void
Bamline_print_new_cigar (FILE *fp, Bamline_T this, Genomicpos_T chrpos_low, char *new_cigar,
			 char *new_md_string, int quality_score_adj);
extern void
Bamline_print_new_mate (FILE *fp, Bamline_T this, char *mate_chr, Genomicpos_T mate_chrpos_low,
			int insert_length);

extern int
Bamline_nm (Bamline_T this);
extern char
Bamline_splice_strand (Bamline_T this);

extern char
Bamline_strand (Bamline_T this, Genome_T genome, IIT_T chromosome_iit);
extern Genomicpos_T
Bamline_chrpos_high (Bamline_T this);
extern Genomicpos_T
Bamline_chrpos_high_noclip (Bamline_T this);
extern Genomicpos_T
Bamline_total_ins (Bamline_T this);
extern int
Bamline_nmismatches (Bamline_T this);


extern void
Bamline_free (Bamline_T *old);
extern Bamline_T
Bamread_next_bamline (T this, char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		      bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		      bool need_concordant_p);
extern Bamline_T
Bamread_next_indel_bamline (T this, char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
			    bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
			    bool need_concordant_p);

extern Bamline_T *
Bamread_next_bamline_set (int *nlines, Bamline_T *prev_bamline,
			  T this, char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
			  bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
			  bool need_concordant_p);

extern Bamline_T **
Bamread_block (int **nlines, Genomicpos_T chrstart, Genomicpos_T chrend,
	       T this, char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	       bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	       bool need_concordant_p);

extern Bamline_T
Bamread_get_acc (T this, char *desired_chr, Genomicpos_T desired_chrpos, char *desired_acc);


typedef struct Bamstore_T *Bamstore_T;

extern void
Bamstore_free (Bamstore_T *old);

extern Bamstore_T
Bamstore_new (Genomicpos_T chrpos);

extern Bamline_T
Bamstore_get (Table_T bamstore_chrtable, char *chr, Genomicpos_T low, char *acc,
	      Genomicpos_T mate_low, int hiti);

extern void
Bamstore_add_at_low (Table_T bamstore_chrtable, char *chr, Genomicpos_T low,
		     Bamline_T bamline);

extern void
Bamstore_table_free (Uinttable_T *bamstore_table);


typedef struct Bampair_T *Bampair_T;

extern char *
Bampair_acc (Bampair_T this);
extern Bamline_T
Bampair_bamline_low (Bampair_T this);
extern Bamline_T
Bampair_bamline_high (Bampair_T this);
extern Genomicpos_T
Bampair_chrpos_low (Bampair_T this);
extern Genomicpos_T
Bampair_chrpos_high (Bampair_T this);
extern Genomicpos_T
Bampair_chrpos_low_noclip (Bampair_T this);
extern Genomicpos_T
Bampair_chrpos_high_noclip (Bampair_T this);

extern int
Bampair_level (Bampair_T this);
extern bool
Bampair_plusp (Bampair_T this);
extern bool
Bampair_good_unique_p (Bampair_T this);
extern bool
Bampair_uniquep (Bampair_T this);
extern bool
Bampair_primaryp (Bampair_T this);

extern void
Bampair_free (Bampair_T *old);
extern void
Bampair_print (FILE *fp, Bampair_T this, int quality_score_adj);
extern void
Bampair_details (Uintlist_T *chrpos_first_lows, Uintlist_T *chrpos_first_highs,
		 Uintlist_T *chrpos_second_lows, Uintlist_T *chrpos_second_highs,
		 Uintlist_T *chrpos_overlap_lows, Uintlist_T *chrpos_overlap_highs,
		 Uintlist_T *splice_lows, Uintlist_T *splice_highs, Intlist_T *splice_signs,
		 Bampair_T this);

extern List_T
Bamread_all_pairs (T bamreader, char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		   bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
		   bool need_concordant_p);

extern int
Bampair_compute_levels (List_T bampairs, Genomicpos_T mincoord,
			Genomicpos_T maxcoord, int max_allowed_levels,
			double xfactor, Genomicpos_T min_pairlength, bool only_internal_p);


#undef T
#endif

