#include <stdlib.h> /* for abs */
#include <string.h> /* for strlen */

#include "iit.h"
#include "bytestream.h"

enum { SEQNAMES, POS, REF, READ, N_CYCLES, N_CYCLES_REF, COUNT, COUNT_REF,
       COUNT_TOTAL, HIGH_QUALITY, HIGH_QUALITY_REF, HIGH_QUALITY_TOTAL,
       MEAN_QUALITY, MEAN_QUALITY_REF, COUNT_PLUS, COUNT_PLUS_REF, COUNT_MINUS,
       COUNT_MINUS_REF, N_BASE_COLS };

typedef struct TallyTable {
  SEXP seqnames_R;
  int *pos;
  SEXP ref_R;
  SEXP read_R;
  int *n_cycles;
  int *n_cycles_ref;
  int *count;
  int *count_ref;
  int *count_total;
  int *high_quality;
  int *high_quality_ref;
  int *high_quality_total;
  double *mean_quality;
  double *mean_quality_ref;
  int *count_plus;
  int *count_plus_ref;
  int *count_minus;
  int *count_minus_ref;
  int **cycle_bins;
} TallyTable;

static SEXP R_TallyTable_new(int n_rows, int n_cycle_bins) {
  SEXP tally_R; /* the result list */
  PROTECT(tally_R = allocVector(VECSXP, N_BASE_COLS + n_cycle_bins));
  
  SET_VECTOR_ELT(tally_R, SEQNAMES, allocVector(STRSXP, n_rows));
  SET_VECTOR_ELT(tally_R, POS, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, REF, allocVector(STRSXP, n_rows));
  SET_VECTOR_ELT(tally_R, READ, allocVector(STRSXP, n_rows));
  SET_VECTOR_ELT(tally_R, N_CYCLES, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, N_CYCLES_REF, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_REF, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_TOTAL, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, HIGH_QUALITY, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, HIGH_QUALITY_REF, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, HIGH_QUALITY_TOTAL, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, MEAN_QUALITY, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, MEAN_QUALITY_REF, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_PLUS, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_PLUS_REF, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_MINUS, allocVector(INTSXP, n_rows));
  SET_VECTOR_ELT(tally_R, COUNT_MINUS_REF, allocVector(INTSXP, n_rows));
  for (int bin = 0; bin < n_cycle_bins; bin++) {
    SEXP cycle_bin_R = allocVector(INTSXP, n_rows);
    SET_VECTOR_ELT(tally_R, bin + N_BASE_COLS, cycle_bin_R);
  }

  UNPROTECT(1);
  return tally_R;
}

static TallyTable *TallyTable_new(SEXP tally_R) {
  TallyTable *tally = (TallyTable *) R_alloc(sizeof(TallyTable), 1);
  int n_cycle_bins = length(tally_R) - N_BASE_COLS;
  
  tally->seqnames_R = VECTOR_ELT(tally_R, SEQNAMES);
  tally->pos = INTEGER(VECTOR_ELT(tally_R, POS));
  tally->ref_R = VECTOR_ELT(tally_R, REF);
  tally->read_R = VECTOR_ELT(tally_R, READ);
  tally->n_cycles = INTEGER(VECTOR_ELT(tally_R, N_CYCLES));
  tally->n_cycles_ref = INTEGER(VECTOR_ELT(tally_R, N_CYCLES_REF));
  tally->count = INTEGER(VECTOR_ELT(tally_R, COUNT));
  tally->count_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_REF));
  tally->count_total = INTEGER(VECTOR_ELT(tally_R, COUNT_TOTAL));
  tally->high_quality = INTEGER(VECTOR_ELT(tally_R, HIGH_QUALITY));
  tally->high_quality_ref = INTEGER(VECTOR_ELT(tally_R, HIGH_QUALITY_REF));
  tally->high_quality_total = INTEGER(VECTOR_ELT(tally_R, HIGH_QUALITY_TOTAL));
  tally->mean_quality = REAL(VECTOR_ELT(tally_R, MEAN_QUALITY));
  tally->mean_quality_ref = REAL(VECTOR_ELT(tally_R, MEAN_QUALITY_REF));
  tally->count_plus = INTEGER(VECTOR_ELT(tally_R, COUNT_PLUS));
  tally->count_plus_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_PLUS_REF));
  tally->count_minus = INTEGER(VECTOR_ELT(tally_R, COUNT_MINUS));
  tally->count_minus_ref = INTEGER(VECTOR_ELT(tally_R, COUNT_MINUS_REF));
  tally->cycle_bins = (int **) R_alloc(sizeof(int*), n_cycle_bins);
  for (int bin = 0; bin < n_cycle_bins; bin++) {
    tally->cycle_bins[bin] = INTEGER(VECTOR_ELT(tally_R, bin + N_BASE_COLS));
  }

  return tally;
}

static void
read_total_counts(unsigned char **bytes, int row, int *count_total) {
  int count_total_plus = read_int(bytes);
  int count_total_minus = read_int(bytes);
  count_total[row] = count_total_plus + count_total_minus;
}

static void
read_cycle_counts(unsigned char **bytes, int row, int *n_cycles,
                  int **cycle_bins, int *cycle_breaks, int n_cycle_bins)
{
  int n_cycle_breaks = n_cycle_bins + 1;
  n_cycles[row] = read_int(bytes);
  if (cycle_breaks != NULL) {
    for (int index = 0; index < n_cycles[row]; index++) {
      int cycle = abs(read_int(bytes));
      int count = read_int(bytes);
      int bin = 0;
      while(n_cycle_breaks > bin &&
            cycle > cycle_breaks[bin])
        bin++;
      if (bin > 0 && bin < n_cycle_breaks) {
        cycle_bins[bin-1][row] += count;
      }
    }
  } else {
    (*bytes) += n_cycles[row] * 2 * sizeof(int);
  }
}

static int
parse_indels(unsigned char *bytes, int row,
             int *cycle_breaks, int n_cycle_bins,
             int high_base_quality, TallyTable *tally, bool insertion)
{
  int indel_count = read_int(&bytes);
  for (int indel = 0; indel < indel_count; indel++) {
    /* TODO: add ref counts when we have them from Tom. */
    /* TODO: sum the ref counts with the indel counts for total. */
    tally->count_plus[row] = read_int(&bytes);
    tally->count_minus[row] = read_int(&bytes);
    tally->count[row] = tally->count_plus[row] + tally->count_minus[row];
    SEXP seq_R = mkChar(read_string(&bytes));
    if (insertion) {
      SET_STRING_ELT(tally->read_R, row, seq_R);
      SET_STRING_ELT(tally->ref_R, row, R_BlankString);
    }
    else {
      SET_STRING_ELT(tally->read_R, row, R_BlankString);
      SET_STRING_ELT(tally->ref_R, row, seq_R);
    }
    read_cycle_counts(&bytes, row, tally->n_cycles, tally->cycle_bins,
                      cycle_breaks, n_cycle_bins);
    /* quality is per-base and so not relevant for indels */
    tally->mean_quality[row] = NA_REAL;
    tally->mean_quality_ref[row] = NA_REAL;
    tally->high_quality[row] = NA_INTEGER;
    tally->high_quality_ref[row] = NA_INTEGER;
    tally->high_quality_total[row] = NA_INTEGER;
  }
  return indel_count;
}

static void
read_quality_counts(unsigned char **bytes, int row, int *high_quality,
                    double *mean_quality, int high_base_quality)
{
  int n_qualities = read_int(bytes);
  int total_quality = 0;
  int total_quality_weight = 0;
  high_quality[row] = 0;
  for (int index = 0; index < n_qualities; index++) {
    int quality = read_int(bytes);
    int count = read_int(bytes);
    if (quality > high_base_quality) {
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
parse_alleles(unsigned char *bytes, int row, int ref_row, int *cycle_breaks,
              int n_cycle_bins, int high_base_quality, TallyTable *tally)
{
  read_total_counts(&bytes, row, tally->count_total);
  int n_alleles = read_allele_counts(&bytes, row, tally->read_R,
                                     tally->count_plus, tally->count_minus,
                                     tally->count);
  for (int allele = 0; allele < n_alleles; allele++, row++) {
    tally->n_cycles[row] = 0;
    for (int b = 0; b < n_cycle_bins; b++) {
      tally->cycle_bins[b][row] = 0;
    }
    tally->high_quality[row] = 0;
    tally->mean_quality[row] = R_NaN;
    if (tally->count[row] > 0) {
      read_cycle_counts(&bytes, row, tally->n_cycles, tally->cycle_bins,
                        cycle_breaks, n_cycle_bins);
      read_quality_counts(&bytes, row, tally->high_quality, tally->mean_quality,
                          high_base_quality);
    }
  }
  int high_quality_total = 0;
  for (int r = ref_row; r < row; r++) {
    high_quality_total += tally->high_quality[r];
  }
  for (int r = ref_row; r < row; r++) {
    tally->n_cycles_ref[r] = tally->n_cycles[ref_row];
    tally->count_total[r] = tally->count_total[ref_row];
    tally->mean_quality_ref[r] = tally->mean_quality[ref_row];
    tally->high_quality_ref[r] = tally->high_quality[ref_row];
    tally->high_quality_total[r] = high_quality_total;
    tally->count_plus_ref[r] = tally->count_plus[ref_row];
    tally->count_minus_ref[r] = tally->count_minus[ref_row];
    tally->count_ref[r] = tally->count_plus_ref[r] + tally->count_minus_ref[r];
    SET_STRING_ELT(tally->ref_R, r, STRING_ELT(tally->read_R, ref_row));  
  }
  bool have_ref_row = row > ref_row;
  if (have_ref_row) {
    /* clear the 'alt' columns for the 'ref' row with NAs */
    SET_STRING_ELT(tally->read_R, ref_row, NA_STRING);
    tally->n_cycles[ref_row] = NA_INTEGER;
    tally->mean_quality[ref_row] = NA_REAL;
    tally->high_quality[ref_row] = NA_REAL;
    tally->count_plus[ref_row] = NA_INTEGER;
    tally->count_minus[ref_row] = NA_INTEGER;
    tally->count[ref_row] = NA_INTEGER;
  }
  return n_alleles;
}
    
static int parse_indel_count(unsigned char *bytes) {
  int count = read_int(&bytes);
  return count;
}

static int parse_allele_count(unsigned char *bytes) {
  int n_alleles = 1; /* always have a reference */
  bytes += sizeof(int) * 4 + 1; /* skip total and reference */
  while(bytes[0] != '\0') {
    bytes += sizeof(int) * 2 + 1;
    n_alleles++;
  }
  return n_alleles;
}

static int count_rows_for_interval(IIT_T tally_iit, int index) {
  int n_rows;
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
  return n_rows;
}

static int parse_interval(IIT_T tally_iit, int index,
                          TallyTable *tally, int row,
                          int *cycle_breaks, int n_cycle_bins,
                          int high_base_quality)
{
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
      row += parse_alleles(base + allele_offset, row, ref_row,
                           cycle_breaks, n_cycle_bins,
                           high_base_quality, tally);
    if (deletion_offset - insertion_offset > 0)
      row += parse_indels(base + insertion_offset, row,
                          cycle_breaks, n_cycle_bins,
                          high_base_quality, tally, true);
    if (allele_offset - deletion_offset > 0)
      row += parse_indels(base + deletion_offset, row,
                          cycle_breaks, n_cycle_bins,
                          high_base_quality, tally, false);
    /* fill in position information */
    for (int r = ref_row; r < row; r++) {
      SET_STRING_ELT(tally->seqnames_R, r, divstring_R);
      tally->pos[r] = start + position;
    }      
  }
  UNPROTECT(1);
  return row;
}

static SEXP parse_all(IIT_T tally_iit, int *cycle_breaks,
                      int n_cycle_bins,
                      int high_base_quality)
{
  int n_rows = 0;
  /* loop over the IIT, getting the total number of rows
     this is (num alts + 1) for every position. */
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    n_rows += count_rows_for_interval(tally_iit, index);
  }
  
  SEXP tally_R;
  PROTECT(tally_R = R_TallyTable_new(n_rows, n_cycle_bins));
  TallyTable *tally = TallyTable_new(tally_R);
  
  int row = 0;
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    row = parse_interval(tally_iit, index, tally, row,
                         cycle_breaks, n_cycle_bins,
                         high_base_quality);
  }

  UNPROTECT(1);
  return tally_R;
}

static SEXP parse_some(IIT_T tally_iit, int *cycle_breaks,
                       int n_cycle_bins,
                       int high_base_quality,
                       SEXP chr_R, int *start, int *end)
{
  int n_rows = 0;
  int *nmatches = (int *) R_alloc(sizeof(int), length(chr_R));
  int **indexes = (int **) R_alloc(sizeof(int*), length(chr_R));
  
  for (int i = 0; i < length(chr_R); i++) {
    indexes[i] = IIT_get(nmatches + i, tally_iit,
                         (char *) CHAR(STRING_ELT(chr_R, i)),
                         start[i], end[i], false);
    for (int j = 0; j < nmatches[i]; j++) {
      n_rows += count_rows_for_interval(tally_iit, indexes[i][j]);
    }
  }
  
  SEXP tally_R;
  PROTECT(tally_R = R_TallyTable_new(n_rows, n_cycle_bins));
  TallyTable *tally = TallyTable_new(tally_R);

  int row = 0;
  for (int i = 0; i < length(chr_R); i++) {
    for (int j = 0; j < nmatches[i]; j++) {
      row = parse_interval(tally_iit, indexes[i][j],
                           tally, row,
                           cycle_breaks, n_cycle_bins,
                           high_base_quality);
    }
  }
  
  UNPROTECT(1);
  return tally_R;
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

SEXP R_tally_iit_parse(SEXP tally_iit_R, SEXP cycle_breaks_R,
                       SEXP high_base_quality_R, SEXP which_R)
{
  int high_base_quality = asInteger(high_base_quality_R);
  int *cycle_breaks =
    cycle_breaks_R == R_NilValue ? NULL : INTEGER(cycle_breaks_R);
  int n_cycle_bins =
    length(cycle_breaks_R) == 0 ? 0 : length(cycle_breaks_R) - 1;
  IIT_T tally_iit = (IIT_T) R_ExternalPtrAddr(tally_iit_R);
  SEXP tally_R;
  
  if (which_R == R_NilValue) {
    tally_R = parse_all(tally_iit, cycle_breaks, n_cycle_bins,
                        high_base_quality);
  } else {
    SEXP chr_R = VECTOR_ELT(which_R, 0);
    int *start = INTEGER(VECTOR_ELT(which_R, 1));
    int *end = INTEGER(VECTOR_ELT(which_R, 2));
    tally_R = parse_some(tally_iit, cycle_breaks, n_cycle_bins,
                         high_base_quality, chr_R, start, end);
  }
  
  return tally_R;
}
