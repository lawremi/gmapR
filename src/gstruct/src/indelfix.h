#ifndef INDELFIX_INCLUDED
#define INDELFIX_INCLUDED

#include "bool.h"
#include "bamread.h"
#include "genome.h"
#include "genomicpos.h"


extern void
Indelfix_run (Bamreader_T bamreader, int max_rlength, int max_glength, Genome_T genome, char *printchr,
	      Genomicpos_T chroffset, Genomicpos_T chrstart, Genomicpos_T chrend,
	      char *desired_read_group, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
	      bool need_concordant_p, bool need_unique_p, bool need_primary_p, bool ignore_duplicates_p,
	      bool ignore_lowend_p, bool ignore_highend_p, Genomicpos_T neighborhood_size,
	      bool allow_multiple_p, char *bamfile);

#endif

