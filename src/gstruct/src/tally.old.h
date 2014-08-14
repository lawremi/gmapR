#ifndef TALLY_INCLUDED
#define TALLY_INCLUDED
#include "iit-read.h"
#include "genomicpos.h"

extern long int
Tally_store_counts (long int *counts_tally, IIT_T tally_iit,
		    char *chr, unsigned int coordstart, unsigned int coordend);
extern long int
Tally_store_runlength (long int *counts_tally, IIT_T runlength_iit,
		       char *chr, unsigned int coordstart, unsigned int coordend);
extern void
Tally_add_runlength (long int *counts, IIT_T runlength_iit,
		     char *chr, unsigned int coordstart, unsigned int coordend);

extern void
Tally_compute_log (double *log_tally, int *counts_tally, Genomicpos_T chrstart, Genomicpos_T chrend);
extern long int
Tally_sum (long int *maxcount, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern double
Tally_mean (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern double
Tally_mean_double (double *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern long int
Tally_median (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern double
Tally_median_double (double *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern long int
Tally_quantile (long int *x, Genomicpos_T chrstart, Genomicpos_T chrend, double percentile);
extern void
Tally_stats (long int *minx, long int *maxx, double *mean, double *sdev,
	     long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);

extern void
Tally_range (long int *mincount, long int *maxcount, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern long int
Tally_maxcount (long int *x, Genomicpos_T chrlength, Genomicpos_T chrstart, Genomicpos_T chrend);

extern void
Tally_cumulate_int (long int *cum, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern void
Tally_cumulate_int_to_double (double *cum, long int *x, Genomicpos_T chrstart, Genomicpos_T chrend);
extern void
Tally_cumulate_double (double *cum, double *x, Genomicpos_T chrstart, Genomicpos_T chrend);

extern bool
Tally_introntest (long int *counts, Genomicpos_T low, Genomicpos_T high,
		  double percentile, int testregion);
extern bool
Tally_exontest_leftward (long int *counts, Genomicpos_T exonlow, Genomicpos_T exonhigh,
			 double percentile, int testregion);
extern bool
Tally_exontest_rightward (long int *counts, Genomicpos_T exonlow, Genomicpos_T exonhigh,
			  double percentile, int testregion);

extern Genomicpos_T
Tally_trace_leftward_exontest (long int *counts, Genomicpos_T exonhigh, int min_exonlength, int max_exonlength,
			       double percentile, int testregion);
extern Genomicpos_T
Tally_trace_rightward_exontest (long int *counts, Genomicpos_T exonlow, int min_exonlength, int max_exonlength,
				double percentile, int testregion);

extern Genomicpos_T
Tally_trace_leftward_pairing (long int *counts, Genomicpos_T low, Genomicpos_T high,
			      int min_exonlength, int max_exonlength,
			      double foldchange_down, double foldchange_up, int testregion);
extern Genomicpos_T
Tally_trace_rightward_pairing (long int *counts, Genomicpos_T low, Genomicpos_T high,
			       int min_exonlength, int max_exonlength,
			       double foldchange_down, double foldchange_up, int testregion);

extern bool
Tally_gaptest (long int *counts, Genomicpos_T low, Genomicpos_T high,
	       int testregion);

extern void
Tally_median_filter (long int *new_tally, long int *raw_tally, long int maxcount, Genomicpos_T median_halfwidth,
		     Genomicpos_T chrstart, Genomicpos_T chrend, Genomicpos_T chrlength);

extern void
Tally_setup_xcoef ();
extern Genomicpos_T
Tally_trace_leftward_linearfit (Genomicpos_T startpos, long int *counts_median, long int *cum_median,
				Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance);
extern Genomicpos_T
Tally_trace_leftward_loglinearfit (Genomicpos_T startpos, double *log_median, double *cumlog_median,
				   Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance);
extern Genomicpos_T
Tally_trace_rightward_linearfit (Genomicpos_T startpos, long int *counts_median, long int *cum_median,
				 Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance);
extern Genomicpos_T
Tally_trace_rightward_loglinearfit (Genomicpos_T startpos, double *log_median, double *cumlog_median,
				    Genomicpos_T chrstart, Genomicpos_T chrend, double tolerance);


extern Genomicpos_T
Tally_trace_leftward_loglik (Genomicpos_T downpos, long int *counts_tally, long int *cum_tally,
			     Genomicpos_T chrstart, Genomicpos_T chrend, 
			     int min_exonlength, int max_exonlength);
extern Genomicpos_T
Tally_trace_rightward_loglik (Genomicpos_T uppos, long int *counts_tally, long int *cum_tally,
			      Genomicpos_T chrstart, Genomicpos_T chrend,
			      int min_exonlength, int max_exonlength);

extern void
Tally_prior_init (int tau, int max_exonlength);

extern Genomicpos_T
Tally_trace_leftward_varmean (Genomicpos_T downpos, long int *counts_tally,
			      Genomicpos_T chrstart, Genomicpos_T chrend, 
			      int auto_exonlength, int max_exonlength, bool debugp);
extern Genomicpos_T
Tally_trace_rightward_varmean (Genomicpos_T uppos, long int *counts_tally,
			       Genomicpos_T chrstart, Genomicpos_T chrend, 
			       int auto_exonlength, int max_exonlength, bool debugp);

#endif

