#include <stdlib.h> /* for abs */
#include <string.h> /* for strlen */

#include "iit.h"
#include "bytestream.h"

enum { SEQNAMES, POS, REF, READ, N_CYCLES, N_CYCLES_REF, COUNT, COUNT_REF,
       COUNT_TOTAL, HIGH_QUALITY, HIGH_QUALITY_REF, HIGH_QUALITY_TOTAL,
       MEAN_QUALITY, MEAN_QUALITY_REF, COUNT_PLUS, COUNT_PLUS_REF, COUNT_MINUS,
       COUNT_MINUS_REF, READ_POS_MEAN, READ_POS_MEAN_REF, READ_POS_VAR,
       READ_POS_VAR_REF, MDFNE, MDFNE_REF, N_BASE_COLS };

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
  double *read_pos_mean;
  double *read_pos_mean_ref;
  double *read_pos_var;
  double *read_pos_var_ref;
  double *mdfne;
  double *mdfne_ref;
  int **cycle_bins;
} TallyTable;

typedef struct TallyParam {
  int *cycle_breaks;
  int n_cycle_bins;
  int high_base_quality;
  int read_length;
  double *mdfne_buf;
} TallyParam;

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
  SET_VECTOR_ELT(tally_R, READ_POS_MEAN, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, READ_POS_MEAN_REF, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, READ_POS_VAR, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, READ_POS_VAR_REF, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, MDFNE, allocVector(REALSXP, n_rows));
  SET_VECTOR_ELT(tally_R, MDFNE_REF, allocVector(REALSXP, n_rows));

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
  tally->read_pos_mean = REAL(VECTOR_ELT(tally_R, READ_POS_MEAN));
  tally->read_pos_mean_ref = REAL(VECTOR_ELT(tally_R, READ_POS_MEAN_REF));
  tally->read_pos_var = REAL(VECTOR_ELT(tally_R, READ_POS_VAR));
  tally->read_pos_var_ref = REAL(VECTOR_ELT(tally_R, READ_POS_VAR_REF));
  tally->mdfne = REAL(VECTOR_ELT(tally_R, MDFNE));
  tally->mdfne_ref = REAL(VECTOR_ELT(tally_R, MDFNE_REF));
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
read_cycle_counts(unsigned char **bytes, int row, TallyParam param,
                  int *n_cycles, double *read_pos_mean, double *read_pos_var,
                  double *mdfne, int **cycle_bins)
{
  int n_cycle_breaks = param.n_cycle_bins + 1;
  int count_sum = 0, weighted_sum = 0, weighted_sum_sq = 0;
  int midpoint = param.read_length / 2.0 + 0.5;

  if (param.mdfne_buf != NULL) {
    memset(param.mdfne_buf, 0, sizeof(int) * midpoint);
  }
  n_cycles[row] = read_int(bytes);
  
  for (int index = 0; index < n_cycles[row]; index++) {
    int cycle = abs(read_int(bytes));
    int count = read_int(bytes);
    int bin = 0;
    if (param.mdfne_buf != NULL && cycle <= param.read_length) {
      param.mdfne_buf[midpoint-abs(cycle-midpoint)-1] = count;
    }
    count_sum += count;
    weighted_sum += cycle * count;
    weighted_sum_sq += cycle * cycle * count;
    if (param.n_cycle_bins > 0) {
      while(n_cycle_breaks > bin &&
            cycle > param.cycle_breaks[bin])
        bin++;
      if (bin > 0 && bin < n_cycle_breaks) {
        cycle_bins[bin-1][row] += count;
      }
    }
  }

  read_pos_mean[row] = ((double)weighted_sum) / count_sum;
  if (count_sum > 1) {
    read_pos_var[row] = ((double)weighted_sum_sq) / (count_sum - 1) -
      (count_sum / (count_sum - 1)) * read_pos_mean[row] * read_pos_mean[row];
  } else {
    read_pos_var[row] = NA_REAL;
  }
  
  mdfne[row] = NA_REAL;
  if (param.mdfne_buf != NULL) {
    int prev_dist = 0;
    int total_count = count_sum;
    count_sum = 0;
    for (int dist = 0; dist < midpoint; dist++) {
      count_sum += param.mdfne_buf[dist];
      if (count_sum > total_count / 2) {
        if (total_count % 2 == 0) {
          mdfne[row] = prev_dist + (dist - prev_dist) / 2.0;
        } else {
          mdfne[row] = dist;
        }
        break;
      }
      if (param.mdfne_buf[dist] > 0) {
        prev_dist = dist;
      }
    }
  }
}

static int
parse_indels(unsigned char *bytes, int row,
             TallyParam param, TallyTable *tally, bool insertion)
{
  int indel_count = read_int(&bytes);
  for (int indel = 0; indel < indel_count; indel++, row++) {
    for (int b = 0; b < param.n_cycle_bins; b++) {
      tally->cycle_bins[b][row] = 0;
    }
    tally->count_plus[row] = read_int(&bytes);
    tally->count_minus[row] = read_int(&bytes);
    tally->count[row] = tally->count_plus[row] + tally->count_minus[row];
    tally->count_plus_ref[row] = read_int(&bytes);
    tally->count_minus_ref[row] = read_int(&bytes);
    tally->count_ref[row] = tally->count_plus_ref[row] +
      tally->count_minus_ref[row];
    tally->count_total[row] = tally->count_ref[row] + tally->count[row];
    SEXP seq_R = mkChar(read_string(&bytes));
    if (insertion) {
      SET_STRING_ELT(tally->read_R, row, seq_R);
      SET_STRING_ELT(tally->ref_R, row, R_BlankString);
    }
    else {
      SET_STRING_ELT(tally->read_R, row, R_BlankString);
      SET_STRING_ELT(tally->ref_R, row, seq_R);
    }
    read_cycle_counts(&bytes, row, param, tally->n_cycles,
                      tally->read_pos_mean, tally->read_pos_var,
                      tally->mdfne, tally->cycle_bins);
    /* quality is per-base and so not relevant for indels */
    tally->mean_quality[row] = NA_REAL;
    tally->mean_quality_ref[row] = NA_REAL;
    tally->high_quality[row] = NA_INTEGER;
    tally->high_quality_ref[row] = NA_INTEGER;
    tally->high_quality_total[row] = NA_INTEGER;
    /* no position from which to tabulate the cycles */
    tally->n_cycles_ref[row] = NA_INTEGER;
    tally->read_pos_mean_ref[row] = NA_REAL;
    tally->read_pos_var_ref[row] = NA_REAL;
    tally->mdfne_ref[row] = NA_REAL;
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
    if (quality >= high_base_quality) {
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
parse_alleles(unsigned char *bytes, int row, int ref_row,
              TallyParam param, TallyTable *tally)
{
  read_total_counts(&bytes, row, tally->count_total);
  int n_alleles = read_allele_counts(&bytes, row, tally->read_R,
                                     tally->count_plus, tally->count_minus,
                                     tally->count);
  for (int allele = 0; allele < n_alleles; allele++, row++) {
    tally->n_cycles[row] = 0;
    for (int b = 0; b < param.n_cycle_bins; b++) {
      tally->cycle_bins[b][row] = 0;
    }
    tally->high_quality[row] = 0;
    tally->mean_quality[row] = R_NaN;
    if (tally->count[row] > 0) {
      read_cycle_counts(&bytes, row, param, tally->n_cycles,
                        tally->read_pos_mean, tally->read_pos_var,
                        tally->mdfne, tally->cycle_bins);
      read_quality_counts(&bytes, row, tally->high_quality, tally->mean_quality,
                          param.high_base_quality);
    } else {
      tally->read_pos_mean[row] = R_NaN;
      tally->read_pos_var[row] = NA_REAL;
      tally->mdfne[row] = NA_REAL;
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
    tally->read_pos_mean_ref[r] = tally->read_pos_mean[ref_row];
    tally->read_pos_var_ref[r] = tally->read_pos_var[ref_row];
    tally->mdfne_ref[r] = tally->mdfne[ref_row];
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
    tally->read_pos_mean[ref_row] = NA_REAL;
    tally->read_pos_var[ref_row] = NA_REAL;
    tally->mdfne[ref_row] = NA_REAL;
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
  int n_rows = 0;
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
                          TallyParam param, TallyTable *tally, int row)
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
                           param, tally);
    if (deletion_offset - insertion_offset > 0)
      row += parse_indels(base + insertion_offset, row,
                          param, tally, true);
    if (allele_offset - deletion_offset > 0)
      row += parse_indels(base + deletion_offset, row,
                          param, tally, false);
    /* fill in position information */
    for (int r = ref_row; r < row; r++) {
      SET_STRING_ELT(tally->seqnames_R, r, divstring_R);
      tally->pos[r] = start + position;
    }      
  }
  UNPROTECT(1);
  return row;
}

static SEXP parse_all(IIT_T tally_iit, TallyParam param)
{
  int n_rows = 0;
  /* loop over the IIT, getting the total number of rows
     this is (num alts + 1) for every position. */
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    n_rows += count_rows_for_interval(tally_iit, index);
  }
  
  SEXP tally_R;
  PROTECT(tally_R = R_TallyTable_new(n_rows, param.n_cycle_bins));
  TallyTable *tally = TallyTable_new(tally_R);
  
  int row = 0;
  for (int index = 1; index <= IIT_total_nintervals(tally_iit); index++) {
    row = parse_interval(tally_iit, index, param, tally, row);
  }

  UNPROTECT(1);
  return tally_R;
}

static SEXP parse_some(IIT_T tally_iit, TallyParam param,
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
  PROTECT(tally_R = R_TallyTable_new(n_rows, param.n_cycle_bins));
  TallyTable *tally = TallyTable_new(tally_R);

  int row = 0;
  for (int i = 0; i < length(chr_R); i++) {
    for (int j = 0; j < nmatches[i]; j++) {
      row = parse_interval(tally_iit, indexes[i][j], param, tally, row);
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
                       SEXP high_base_quality_R, SEXP which_R,
                       SEXP read_length_R)
{
  IIT_T tally_iit = (IIT_T) R_ExternalPtrAddr(tally_iit_R);
  SEXP tally_R;
  TallyParam param;
  param.high_base_quality = asInteger(high_base_quality_R);
  param.cycle_breaks =
    cycle_breaks_R == R_NilValue ? NULL : INTEGER(cycle_breaks_R);
  param.n_cycle_bins =
    length(cycle_breaks_R) == 0 ? 0 : length(cycle_breaks_R) - 1;
  param.read_length = asInteger(read_length_R);
  if (param.read_length != NA_INTEGER) {
    param.mdfne_buf = (double *)R_alloc(sizeof(double), param.read_length);
  } else {
    param.mdfne_buf = NULL;
  }
  
  if (which_R == R_NilValue) {
    tally_R = parse_all(tally_iit, param);
  } else {
    SEXP chr_R = VECTOR_ELT(which_R, 0);
    int *start = INTEGER(VECTOR_ELT(which_R, 1));
    int *end = INTEGER(VECTOR_ELT(which_R, 2));
    tally_R = parse_some(tally_iit, param, chr_R, start, end);
  }
  
  return tally_R;
}
