#ifndef GMAPR_H
#define GMAPR_H

#include <Rinternals.h>

/* .Call entry points for the gmapR package */

SEXP
R_Bamread_new (SEXP bamfile_R);

SEXP
R_Bamtally_iit (SEXP bamreader_R, SEXP genome_dir_R, SEXP db_R,
                SEXP which_R, SEXP desired_read_group_R,
                SEXP alloclength_R,
                SEXP minimum_mapq_R, SEXP good_unique_mapq_R,
                SEXP maximum_nhits_R,
                SEXP need_concordant_p_R, SEXP need_unique_p_R,
                SEXP need_primary_p_R, SEXP ignore_duplicates_p_R,
                SEXP min_depth_R, SEXP variant_strands_R,
                SEXP ignore_query_Ns_p,
                SEXP print_indels_p_R,
                SEXP blocksize_R, 
                SEXP verbosep_R, SEXP max_softclip_R,
                SEXP genome_iit_file_R,
                SEXP print_xs_scores_p_R, SEXP print_cycles_p_R,
                SEXP minimum_quality_score_R, SEXP nonconvered_R,
                SEXP print_nm_scores_p_R);

SEXP
R_tally_iit_parse(SEXP tally_iit_R, SEXP cycle_breaks_R,
                  SEXP high_base_quality, SEXP which_R);

SEXP
R_Genome_getSeq (SEXP genome_dir_R, SEXP db_R,
                 SEXP seqnames_R, SEXP start_R, SEXP width_R, SEXP strand_R);

#endif
