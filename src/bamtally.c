#include <stdlib.h>
#include <string.h>

#include <gstruct/bool.h>
#include <gstruct/datadir.h>
#include <gstruct/iit-read.h>
#include <gstruct/genome.h>
#include <gstruct/interval.h>
#include <gstruct/genomicpos.h>
#include <gstruct/bamread.h>
#include <gstruct/bamtally.h>

#include "gmapR.h"

enum { SEQNAMES, POS, REF, READ, N_CYCLES, N_CYCLES_REF, COUNT, COUNT_REF,
         COUNT_TOTAL, HIGH_QUALITY, HIGH_QUALITY_REF, MEAN_QUALITY,
         MEAN_QUALITY_REF, COUNT_PLUS, COUNT_PLUS_REF, COUNT_MINUS,
         COUNT_MINUS_REF, N_BASE_COLS };

static inline int
read_int (unsigned char **bytes) {
  int x = ((int *)*bytes)[0];
  *bytes += sizeof(int);
  return x;
}

static inline char
read_char (unsigned char **bytes) {
  char x = (*bytes)[0];
  (*bytes)++;
  return x;
}

static inline const char *
read_string (unsigned char **bytes) {
  const char *string = *bytes;
  (*bytes) += strlen(string) + 1;
  return string;
}

static int
parse_indel_count(unsigned char *bytes) {
  int count = read_int(&bytes);
  return count;
}

static int
parse_allele_count(unsigned char *bytes) {
  int n_alleles = 1; /* always have a reference */
  bytes += sizeof(int) * 4 + 1; /* skip total and reference */
  while(bytes[0] != '\0') {
    bytes += sizeof(int) * 2 + 1;
    n_alleles++;
  }
  return n_alleles;
}

static void
read_total_counts(unsigned char **bytes, int row, int *count_total) {
  int count_total_plus = read_int(bytes);
  int count_total_minus = read_int(bytes);
  count_total[row] = count_total_plus + count_total_minus;
}

static void
read_cycle_counts(unsigned char **bytes, int row, int *n_cycles,
                  SEXP tally_R, SEXP cycle_breaks_R)
{
  n_cycles[row] = read_int(bytes);
  if (cycle_breaks_R != R_NilValue) {
    for (int index = 0; index < n_cycles[row]; index++) {
      int cycle = abs(read_int(bytes));
      int count = read_int(bytes);
      int bin = 0;
      while(length(cycle_breaks_R) > bin &&
            cycle > INTEGER(cycle_breaks_R)[bin])
        bin++;
      if (bin > 0 && bin < length(cycle_breaks_R)) {
        SEXP bin_vector = VECTOR_ELT(tally_R, N_BASE_COLS + bin - 1);
        INTEGER(bin_vector)[row] += count;
      }
    }
  } else {
    (*bytes) += n_cycles[row] * 2 * sizeof(int);
  }
}

static int
parse_indels(unsigned char *bytes, int row,
             SEXP read_R, int *n_cycles, int *count,
             int *high_quality, double *mean_quality,
             int *count_plus, int *count_minus,
             SEXP cycle_breaks_R, SEXP tally_R)
{
  int indel_count = read_int(&bytes);
  for (int indel = 0; indel < indel_count; indel++) {
    count_plus[row] = read_int(&bytes);
    count_minus[row] = read_int(&bytes);
    count[row] = count_plus[row] + count_minus[row];
    SEXP seq_R = mkChar(read_string(&bytes));
    SET_STRING_ELT(read_R, row, seq_R);
    read_cycle_counts(&bytes, row, n_cycles, tally_R, cycle_breaks_R);
    /* stuff not relevant for indels */
    mean_quality[row] = NA_REAL;
    high_quality[row] = NA_INTEGER;
  }
  return indel_count;
}

static void
read_quality_counts(unsigned char **bytes, int row, int *high_quality,
                    double *mean_quality, int high_quality_cutoff)
{
  int n_qualities = read_int(bytes);
  int total_quality = 0;
  int total_quality_weight = 0;
  high_quality[row] = 0;
  for (int index = 0; index < n_qualities; index++) {
    int quality = read_int(bytes);
    int count = read_int(bytes);
    if (quality > high_quality_cutoff) {
      total_quality += quality * count;
      total_quality_weight += count;
      high_quality[row] += count;
    }
  }
  mean_quality[row] = total_quality_weight > 0 ?
    (double)total_quality / total_quality_weight : R_NaN;
}

static int
read_allele_counts(unsigned char **bytes, int row, SEXP read_R,
                   int *count_plus, int *count_minus, int *count)
{
  int n_alleles = 0;
  char allele;
  while((allele = read_char(bytes)) != '\0') {
    SET_STRING_ELT(read_R, row, mkCharLen(&allele, 1));
    count_plus[row] = read_int(bytes);
    count_minus[row] = read_int(bytes);
    count[row] = count_plus[row] + count_minus[row];
    row++;
    n_alleles++;
  }
  return n_alleles;
}

static int
parse_alleles(unsigned char *bytes, int row,
              SEXP read_R, int *n_cycles, int *count, int *count_total,
              int *high_quality, double *mean_quality, int *count_plus,
              int *count_minus, SEXP cycle_breaks_R, SEXP tally_R,
              int high_quality_cutoff)
{
  read_total_counts(&bytes, row, count_total);
  int n_alleles = read_allele_counts(&bytes, row, read_R,
                                     count_plus, count_minus, count);
  for (int allele = 0; allele < n_alleles; allele++, row++) {
    n_cycles[row] = 0;
    for (int b = 0; b < length(cycle_breaks_R) - 1; b++) {
      INTEGER(VECTOR_ELT(tally_R, N_BASE_COLS + b))[row] = 0;
    }
    high_quality[row] = 0;
    mean_quality[row] = R_NaN;
    if (count[row] > 0) {
      read_cycle_counts(&bytes, row, n_cycles, tally_R, cycle_breaks_R);
      read_quality_counts(&bytes, row, high_quality, mean_quality,
                          high_quality_cutoff);
    }
  }
  return n_alleles;
}


/* FORMAT

   block header:
   [0][ins][del][mm]

   mismatches for each position:
   [t+][t-][rn][r+][r-]([an][a+][a-])*[\0]([c#]([cv][cc])*)*([q#]([qv][qc])*)*

   insertions:
   [i#]([seq][t+][t-][c#]([cv][cc])*)*

   deletions:
   [d#]([seq][t+][t-][c#]([cv][cc])*)*
*/

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
  int minimum_mapq = asInteger(minimum_mapq_R);
  int good_unique_mapq = asInteger(good_unique_mapq_R);
  int maximum_nhits = asInteger(maximum_nhits_R);
  bool need_concordant_p = asLogical(need_concordant_p_R);
  bool need_unique_p = asLogical(need_unique_p_R);
  bool need_primary_p = asLogical(need_primary_p_R);
  int min_depth = asInteger(min_depth_R);
  int variant_strands = asInteger(variant_strands_R);
  bool ignore_query_Ns_p = asLogical(ignore_query_Ns_p_R);
  bool print_indels_p = asLogical(print_indels_p_R);
  int blocksize = asInteger(blocksize_R);
  int verbosep = asLogical(verbosep_R);
  int high_quality_cutoff = asInteger(high_quality_cutoff_R);
  
  char *genomesubdir, *fileroot, *dbversion, *iitfile;
  genomesubdir = Datadir_find_genomesubdir(&fileroot, &dbversion,
                                           (char *)genome_dir, (char *)db);
  
  iitfile = (char *) calloc(strlen(genomesubdir) + strlen("/") +
                            strlen(fileroot) + strlen(".chromosome.iit") + 1,
                            sizeof(char));
  sprintf(iitfile, "%s/%s.chromosome.iit", genomesubdir, fileroot);
  IIT_T chromosome_iit = IIT_read(iitfile, /*name*/NULL, /*readonlyp*/true,
                                  /*divread*/READ_ALL, /*divstring*/NULL,
                                  /*add_iit_p*/false, /*labels_read_p*/true);
  Genome_T genome = Genome_new(genomesubdir, fileroot, /*snps_root*/NULL,
                               /*uncompressedp*/false, /*access*/USE_MMAP_ONLY);
  free(iitfile);
  free(fileroot);
  free(dbversion);
  free(genomesubdir);
  
  const char *chr = NULL;
  Genomicpos_T start = 0;
  Genomicpos_T end = 0;

  if (which_R != R_NilValue) {
    chr = CHAR(asChar(VECTOR_ELT(which_R, 0)));
    start = asInteger(VECTOR_ELT(which_R, 1));
    end = asInteger(VECTOR_ELT(which_R, 2));
  }

  IIT_T tally_iit = Bamtally_iit(bamreader, (char *)chr, start, end,
                                 genome, chromosome_iit,
                                 alloclength, minimum_mapq, good_unique_mapq,
                                 maximum_nhits, need_concordant_p,
                                 need_unique_p, need_primary_p,
                                 min_depth, variant_strands, ignore_query_Ns_p,
                                 print_indels_p, blocksize, verbosep);
  IIT_free(&chromosome_iit);
  Genome_free(&genome);

  if (tally_iit == NULL) {
    error("Could not create tally\n");
  }

  int n_cycle_bins =
    length(cycle_breaks_R) == 0 ? 0 : length(cycle_breaks_R) - 1;

  int n_rows = 0;
  /* loop over the IIT, getting the total number of rows
     this is (num alts + 1) for every position. */
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    unsigned char *bytes = IIT_data(tally_iit, index);
    int width = IIT_length(tally_iit, index);
    unsigned char *base = bytes + (3 * width + 1) * sizeof(int);
    for (int pos = 0; pos < width; pos++) {
      int insertion_offset = read_int(&bytes);
      int deletion_offset = read_int(&bytes);
      int allele_offset = read_int(&bytes);
      int next_offset = read_int(&bytes);
      bytes -= 4; /* rewind from read-ahead */
      if (deletion_offset - insertion_offset > 0) {
        n_rows += parse_indel_count(base + insertion_offset);
      }
      if (allele_offset - deletion_offset > 0) {
        n_rows += parse_indel_count(base + deletion_offset);
      }
      if (next_offset - allele_offset > 0) {
        n_rows += parse_allele_count(base + allele_offset);
      }
    }
  }

  SEXP tally_R; /* the result list */
  PROTECT(tally_R =
          allocVector(VECSXP, N_BASE_COLS + n_cycle_bins));

  SET_VECTOR_ELT(tally_R, SEQNAMES, allocVector(STRSXP, n_rows));
  SEXP seqnames_R = VECTOR_ELT(tally_R, SEQNAMES);
  SET_VECTOR_ELT(tally_R, POS, allocVector(INTSXP, n_rows));
  int *pos = INTEGER(VECTOR_ELT(tally_R, POS));
  SET_VECTOR_ELT(tally_R, REF, allocVector(STRSXP, n_rows));
  SEXP ref_R = VECTOR_ELT(tally_R, REF);
  SET_VECTOR_ELT(tally_R, READ, allocVector(STRSXP, n_rows));
  SEXP read_R = VECTOR_ELT(tally_R, READ);
  SET_VECTOR_ELT(tally_R, N_CYCLES, allocVector(INTSXP, n_rows));
  int *n_cycles = INTEGER(VECTOR_ELT(tally_R, N_CYCLES));
  SET_VECTOR_ELT(tally_R, N_CYCLES_REF, allocVector(INTSXP, n_rows));
  int *n_cycles_ref = INTEGER(VECTOR_ELT(tally_R, N_CYCLES_REF));
  SET_VECTOR_ELT(tally_R, COUNT, allocVector(INTSXP, n_rows));
  int *count = INTEGER(VECTOR_ELT(tally_R, COUNT));
  SET_VECTOR_ELT(tally_R, COUNT_REF, allocVector(INTSXP, n_rows));
  int *count_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_REF));
  SET_VECTOR_ELT(tally_R, COUNT_TOTAL, allocVector(INTSXP, n_rows));
  int *count_total = INTEGER(VECTOR_ELT(tally_R, COUNT_TOTAL));
  SET_VECTOR_ELT(tally_R, HIGH_QUALITY, allocVector(INTSXP, n_rows));
  int *high_quality = INTEGER(VECTOR_ELT(tally_R, HIGH_QUALITY));
  SET_VECTOR_ELT(tally_R, HIGH_QUALITY_REF, allocVector(INTSXP, n_rows));
  int *high_quality_ref = INTEGER(VECTOR_ELT(tally_R, HIGH_QUALITY_REF));
  SET_VECTOR_ELT(tally_R, MEAN_QUALITY, allocVector(REALSXP, n_rows));
  double *mean_quality = REAL(VECTOR_ELT(tally_R, MEAN_QUALITY));
  SET_VECTOR_ELT(tally_R, MEAN_QUALITY_REF, allocVector(REALSXP, n_rows));
  double *mean_quality_ref = REAL(VECTOR_ELT(tally_R, MEAN_QUALITY_REF));
  SET_VECTOR_ELT(tally_R, COUNT_PLUS, allocVector(INTSXP, n_rows));
  int *count_plus = INTEGER(VECTOR_ELT(tally_R, COUNT_PLUS));
  SET_VECTOR_ELT(tally_R, COUNT_PLUS_REF, allocVector(INTSXP, n_rows));
  int *count_plus_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_PLUS_REF));
  SET_VECTOR_ELT(tally_R, COUNT_MINUS, allocVector(INTSXP, n_rows));
  int *count_minus = INTEGER(VECTOR_ELT(tally_R, COUNT_MINUS));
  SET_VECTOR_ELT(tally_R, COUNT_MINUS_REF, allocVector(INTSXP, n_rows));
  int *count_minus_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_MINUS_REF));
  for (int bin = 0; bin < n_cycle_bins; bin++)    
    SET_VECTOR_ELT(tally_R, bin + N_BASE_COLS, allocVector(INTSXP, n_rows));

  int row = 0;
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    unsigned char *bytes = IIT_data(tally_iit, index);
    int width = IIT_length(tally_iit, index);
    unsigned char *base = bytes + (3 * width + 1) * sizeof(int);
    int start = IIT_interval_low(tally_iit, index);
    SEXP divstring_R;
    PROTECT(divstring_R = mkChar(IIT_divstring_from_index(tally_iit, index)));
    
    for (int position = 0; position < width; position++) {
      int insertion_offset = read_int(&bytes);
      int deletion_offset = read_int(&bytes);
      int allele_offset = read_int(&bytes);
      int next_offset = read_int(&bytes);
      int ref_row = row;
      bytes -= 4; /* rewind from read-ahead */
      if (next_offset - allele_offset > 0)
        row += parse_alleles(base + allele_offset, row, read_R, n_cycles,
                             count, count_total,
                             high_quality, mean_quality, count_plus,
                             count_minus, cycle_breaks_R, tally_R,
                             high_quality_cutoff);
      if (deletion_offset - insertion_offset > 0)
        row += parse_indels(base + insertion_offset, row, read_R, n_cycles,
                            count,
                            high_quality, mean_quality,
                            count_plus, count_minus,
                            cycle_breaks_R, tally_R);
      if (allele_offset - deletion_offset > 0)
        row += parse_indels(base + deletion_offset, row, read_R, n_cycles,
                            count,
                            high_quality, mean_quality,
                            count_plus, count_minus,
                            cycle_breaks_R, tally_R);
      /* fill in position and reference information */
      for (int r = ref_row; r < row; r++) {
        SET_STRING_ELT(seqnames_R, r, divstring_R);
        pos[r] = start + position;
        n_cycles_ref[r] = n_cycles[ref_row];
        count_total[r] = count_total[ref_row];
        mean_quality_ref[r] = mean_quality[ref_row];
        high_quality_ref[r] = high_quality[ref_row];
        count_plus_ref[r] = count_plus[ref_row];
        count_minus_ref[r] = count_minus[ref_row];
        count_ref[r] = count_plus_ref[r] + count_minus_ref[r];
        SET_STRING_ELT(ref_R, r, STRING_ELT(read_R, ref_row));  
      }
      bool have_ref_row = row > ref_row;
      if (have_ref_row) {
        /* clear the 'alt' columns for the 'ref' row with NAs */
        SET_STRING_ELT(read_R, ref_row, NA_STRING);
        n_cycles[ref_row] = NA_INTEGER;
        mean_quality[ref_row] = NA_REAL;
        high_quality[ref_row] = NA_REAL;
        count_plus[ref_row] = NA_INTEGER;
        count_minus[ref_row] = NA_INTEGER;
        count[ref_row] = NA_INTEGER;
      }
    }
    UNPROTECT(1);
  }

  IIT_free(&tally_iit);
  UNPROTECT(1);
  return tally_R;
}
