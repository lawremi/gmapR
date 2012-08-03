#ifndef GMAPR_H
#define GMAPR_H

#include <Rinternals.h>

SEXP
R_Bamread_new (SEXP bamfile_R);

SEXP
R_Bamtally_iit (SEXP bamreader_R, SEXP genome_dir_R, SEXP db_R,
                SEXP which_R,
                SEXP cycle_breaks_R, SEXP high_quality_cutoff_R,
                SEXP alloclength_R,
                SEXP minimum_mapq_R, SEXP good_unique_mapq_R,
                SEXP maximum_nhits_R,
                SEXP need_concordant_p_R, SEXP need_unique_p_R,
                SEXP need_primary_p_R,
                SEXP min_depth_R, SEXP variant_strands_R,
                SEXP ignore_query_Ns_p,
                SEXP print_indels_p_R,
                SEXP blocksize_R, 
                SEXP verbosep_R);

#endif