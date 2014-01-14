#include <stdlib.h>

#include <gstruct/bool.h>
#include <gstruct/iit-read.h>
#include <gstruct/genome.h>
#include <gstruct/genomicpos.h>
#include <gstruct/bamread.h>
#include <gstruct/bamtally.h>

#include "gmapR.h"
#include "iit.h"
#include "genome.h"

SEXP
R_Bamtally_iit (SEXP bamreader_R, SEXP genome_dir_R, SEXP db_R,
                SEXP which_R, SEXP desired_read_group_R,
                SEXP alloclength_R,
                SEXP minimum_mapq_R, SEXP good_unique_mapq_R,
                SEXP maximum_nhits_R,
                SEXP need_concordant_p_R, SEXP need_unique_p_R,
                SEXP need_primary_p_R, SEXP ignore_duplicates_p_R,
                SEXP min_depth_R, SEXP variant_strands_R,
                SEXP ignore_query_Ns_p_R,
                SEXP print_indels_p_R,
                SEXP blocksize_R, 
                SEXP verbosep_R)
{
  Bamreader_T bamreader = (Bamreader_T) R_ExternalPtrAddr(bamreader_R);
  const char *genome_dir =
    genome_dir_R == R_NilValue ? NULL : CHAR(asChar(genome_dir_R));
  const char *db = CHAR(asChar(db_R));
  int alloclength = asInteger(alloclength_R);
  const char *desired_read_group =
    desired_read_group_R == R_NilValue ? NULL :
    CHAR(asChar(desired_read_group_R));
  int minimum_mapq = asInteger(minimum_mapq_R);
  int good_unique_mapq = asInteger(good_unique_mapq_R);
  int maximum_nhits = asInteger(maximum_nhits_R);
  bool need_concordant_p = asLogical(need_concordant_p_R);
  bool need_unique_p = asLogical(need_unique_p_R);
  bool need_primary_p = asLogical(need_primary_p_R);
  bool ignore_duplicates_p = asLogical(ignore_duplicates_p_R);
  int min_depth = asInteger(min_depth_R);
  int variant_strands = asInteger(variant_strands_R);
  bool ignore_query_Ns_p = asLogical(ignore_query_Ns_p_R);
  bool print_indels_p = asLogical(print_indels_p_R);
  int blocksize = asInteger(blocksize_R);
  int verbosep = asLogical(verbosep_R);

  Genome_T genome = createGenome(genome_dir, db);
  IIT_T chromosome_iit = readChromosomeIIT(genome_dir, db);
  
  const char *chr = NULL;
  Genomicpos_T start = 0;
  Genomicpos_T end = 0;

  if (which_R != R_NilValue) {
    chr = CHAR(asChar(VECTOR_ELT(which_R, 0)));
    start = asInteger(VECTOR_ELT(which_R, 1));
    end = asInteger(VECTOR_ELT(which_R, 2));
  }

  IIT_T tally_iit = Bamtally_iit(bamreader, (char *)chr, 
                                 /* TODO: bam_lacks_chr */ NULL,
                                 start, end,
                                 genome, chromosome_iit, alloclength,
                                 (char *)desired_read_group,
                                 minimum_mapq, good_unique_mapq,
                                 maximum_nhits, need_concordant_p,
                                 need_unique_p, need_primary_p,
                                 ignore_duplicates_p,
                                 min_depth, variant_strands, ignore_query_Ns_p,
                                 print_indels_p, blocksize, verbosep,
                                 /*readlevel_p*/false);
  IIT_free(&chromosome_iit);
  Genome_free(&genome);

  if (tally_iit == NULL) {
    error("Could not create tally\n");
  }

  return R_IIT_new(tally_iit);
}
