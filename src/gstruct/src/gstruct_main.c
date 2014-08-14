static char rcsid[] = "$Id: gstruct_main.c 138430 2014-06-06 21:26:50Z twu $";
/* Note: Handles only paired-end data */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For qsort */

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "iit-write.h"		/* For IIT_create */
#include "interval.h"
#include "table.h"
#include "uinttable.h"
#include "maxent_hr.h"		/* For Maxent_hr_setup */

#ifdef BAM_INPUT
#include "bamread.h"
#include "multimap.h"
#endif

#include "samflags.h"
#include "samread.h"
#include "genome.h"
#include "complement.h"

#include "bamtally.h"
#include "splice.h"
#include "gstruct.h"
#include "splicegraph.h"
#include "cappaths.h"		/* For Cappaths_setup */

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"



/************************************************************************
 *   Global variables 
 ************************************************************************/

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_range = NULL;
static char *range_chromosome = NULL;	/* To process entire genome, use NULL */
static Genomicpos_T range_chrstart = 0U;
static Genomicpos_T range_chrend = -1U;


/* For tally */
static int quality_score_adj = 0;
static int alloclength = 400000;
static int blocksize = 1000;

static int min_depth = 1;
static int variant_strands = 0;
static bool genomic_diff_p = false;
static bool ignore_query_Ns_p = true;


/* For splicegraph */
static Gene_outputtype_T output_type = OUTPUT_GENES;

static bool altpaths_p = true;

/* According to Sakharkar, 2004, max exon length in human is 11,923
   bp.  However, that is probably an UTR */

static int min_exonlength = 5;	/* NM_012399 has an 8-nt exon */
static int max_exonlength = 10000; /* Can be large, since we use tally to find flat_bound */
static int auto_exonlength = 500;


static int min_intronlength = 30;
static unsigned int max_intronlength = 500000;

#if 0
static int short_form_support = 5;
static int long_form_support = 10;
static double slope_threshold = 0.1;
static int length_threshold = 50;
static int highlow_length = 70;
static int halfdata_genebounds = 20;
static int halfdata_cappaths = 10;
#endif


/* For sam_splices */

#if 0
static bool allow_translocations_p = false;
static bool trust_sam_p = true;
#endif

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = true;
static bool ignore_duplicates_p = false;

static bool need_canonical_p = false;
static Genomicpos_T max_pairlength = 400000;
static Genomicpos_T shortsplicedist = 400000;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;

static int mincount_end_alt = 2;
static int minsupport = 8;


static bool twopassp = false;




static struct option long_options[] = {
  /* Input options */
#if 0
  {"chr", required_argument, 0, 'c'}, /* chromosome */
#endif
  {"range", required_argument, 0, 'r'}, /* user_range, range_chromosome, range_chrstart, range_chrend */

  {"read-group", required_argument, 0, 0},   /* desired_read_group */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"pairmax", required_argument, 0, 0}, /* max_pairlength */

  /* Known genes */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  /* Options for splicegraph */
  {"format", required_argument, 0, 'A'}, /* output_type */

  {"minexonlength", required_argument, 0, 'e'}, /* min_exonlength */
  {"maxexonlength", required_argument, 0, 'E'}, /* min_exonlength */
  {"minintronlength", required_argument, 0, 'i'}, /* min_intronlength */

  {"main-only", no_argument, 0, 0},	/* altpaths_p */
  {"onepass", no_argument, 0, '1'},	/* twopassp */


  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_GSTRUCT\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: bam_gstruct [-d <genome>] [OPTIONS...] <bam file>\n\
\n\
Note: If SAM/BAM output contains XS: tags (generated by GSNAP, for example),\n\
      to indicate the genomic strands of splices, then the -d flag is not necessary\n\
\n\
Input options\n\
  -r, --range=STRING             Process only the given chromosome: or chromosome:start..end\n\
  --read-group=STRING            Require alignments to have this read group in the RG field of\n\
                                   the BAM line\n\
  -q, --mapq=INT                 Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --pairmax=INT                  Expected insert length (reserves memory for this amount, so\n\
                                     alignments longer than this value are discarded) [default=400000]\n\
  -9, --dump                     Dump graph\n\
\n\
  -e, --minexonlength=INT        Minimum exon length (default 30)\n\
  -E, --maxexonlength=INT        Maximum exon length (default 5000)\n\
  -i, --minintronlength=INT      Minimum intron length (default 30)\n\
  -I, --maxintronlength=INT      Maximum intron length (default 100000)\n\
\n\
Output options\n\
  -A, --format=STRING            Output type: splices, paths, genes, diffs\n\
\n\
");

  return;
}





static void
gstruct_solve_onepass (Genome_T genome, IIT_T chromosome_iit, IIT_T knowngenes_iit,
		       List_T bamfiles, Genomicpos_T max_pairlength, Genomicpos_T shortsplicedist,
		       Genomicpos_T max_exonlength, bool need_unique_p, bool need_primary_p,
		       bool need_canonical_p, int mincount_end_alt, int minsupport) {
  Gstruct_T gstruct, knowngenes_gstruct = NULL;
  IIT_T introns_iit = NULL, middle_exons_iit = NULL, end_exons_iit = NULL;
  char *bamfile;
  Bamreader_T bamreader;

  Splicegraph_T splicegraph;
  List_T splices, knowngenes_splices, obssplices, p;
  List_T obsgenes, knowngenes;
  /* List_T intervallist = NULL, labellist = NULL, datalist = NULL; */

  long int *fwd_extents, *rev_extents, *null_extents;
  long int *tally_matches_low, *tally_mismatches_low,
    *tally_matches_high, *tally_mismatches_high, *primary_extents, *crosshyb_extents;

  Genomicpos_T chrlength;
  char *chr;
  int quality_counts_match[256], quality_counts_mismatch[256], i;
  int divno;
  Genomicpos_T chroffset;
  Interval_T interval;

  int readlength, insertlength;


  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }

  if (knowngenes_iit != NULL) {
    knowngenes_gstruct = Gstruct_knowngenes_input(&introns_iit,&middle_exons_iit,&end_exons_iit,
						  knowngenes_iit,genome,chromosome_iit);
  }

  gstruct = Gstruct_bam_input(&readlength,&insertlength,bamfiles,range_chromosome,range_chrstart,range_chrend,
			      /*ngoodhits_low_table*/NULL,/*ngoodhits_high_table*/NULL,
			      /*genes_iit*/NULL,shortsplicedist,max_pairlength,genome,chromosome_iit,
			      desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_unique_p,need_primary_p,
			      ignore_duplicates_p,/*trust_sam_p*/true,need_canonical_p,
			      /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0);

  while ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				 &primary_extents,&crosshyb_extents,&chrlength,gstruct,
				 max_exonlength,mincount_end_alt,minsupport,need_canonical_p,
				 /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0)) != NULL) {

    FREE(null_extents);
    FREE(rev_extents);
    FREE(fwd_extents);
    fprintf(stderr,"Processing chromosome %s\n",chr);

    splicegraph = Splicegraph_new();

    knowngenes = (List_T) NULL;
    knowngenes_splices = (List_T) NULL;
    if (knowngenes_gstruct != NULL) {
      if ((divno = IIT_find_one(chromosome_iit,chr)) >= 0) {
	chrlength = Interval_length(IIT_interval(chromosome_iit,divno));
	knowngenes_splices = Gstruct_get_known(knowngenes_gstruct,chr);
	knowngenes = Splicegraph_solve_known(splicegraph,knowngenes_splices,
					     introns_iit,middle_exons_iit,end_exons_iit,
					     chr,chrlength);
	/* Splicegraph_print_genes(stdout,knowngenes,chr,bamreader); */
      }
    }


    if (output_type == OUTPUT_SPLICES) {
      for (p = splices; p != NULL; p = List_next(p)) {
	Splice_print((Splice_T) List_head(p),chr,/*map_iit*/NULL,/*show_invalid_p*/true);
      }

    } else if ((divno = IIT_find_one(chromosome_iit,chr)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file...skipping\n",chr);

    } else {
      /* Tally */
      interval = IIT_interval(chromosome_iit,divno);
      chroffset = Interval_low(interval);
      chrlength = Interval_length(interval);

      tally_matches_low = (long int *) CALLOC(chrlength+1,sizeof(long int));
      tally_mismatches_low = (long int *) CALLOC(chrlength+1,sizeof(long int));
      tally_matches_high = (long int *) CALLOC(chrlength+1,sizeof(long int));
      tally_mismatches_high = (long int *) CALLOC(chrlength+1,sizeof(long int));

      for (p = bamfiles; p != NULL; p = List_next(p)) {
	bamfile = (char *) List_head(p);
	bamreader = Bamread_new(bamfile);

	Bamread_limit_region(bamreader,chr,/*chrstart*/0U,/*chrend*/chrlength);
	Bamtally_run_lh(&tally_matches_low,&tally_mismatches_low,
			&tally_matches_high,&tally_mismatches_high,
			quality_counts_match,quality_counts_mismatch,
			bamreader,genome,chr,chroffset,/*chrstart*/0U,/*chrend*/chrlength,
			/*map_iit*/NULL,alloclength,
			desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
			need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
			blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
			ignore_query_Ns_p,/*verbosep*/false,/*readlevel_p*/false,
			/*max_softclip*/0,/*print_noncovered_p*/false);
	Bamread_unlimit_region(bamreader);
	Bamread_free(&bamreader);
      }

      obssplices = Splice_valid_list(splices);
      /* knownsplices = Splicegraph_splices(knowngenes); */
      obsgenes = Splicegraph_solve_obs(splicegraph,obssplices,knowngenes,chr,
				       tally_matches_low,tally_mismatches_low,
				       tally_matches_high,tally_mismatches_high,
				       primary_extents,crosshyb_extents,
				       genome,chroffset,chrlength,
				       insertlength,readlength,/*min_overhang*/15,
				       mincount_end_alt,minsupport,need_canonical_p,
				       min_exonlength,max_exonlength,auto_exonlength,
				       min_intronlength,max_intronlength,/*altpaths_p*/true);
      Splicegraph_print_genes(stdout,obsgenes,Gstruct_npairs(gstruct),knowngenes_iit,chr,
			      tally_matches_low,tally_mismatches_low,
			      tally_matches_high,tally_mismatches_high,
			      insertlength,readlength,/*min_overhang*/15);

      List_free(&obssplices);

      Splicegraph_genes_gc(&obsgenes);

      FREE(tally_mismatches_high);
      FREE(tally_matches_high);
      FREE(tally_mismatches_low);
      FREE(tally_matches_low);
    }

    if (knowngenes != NULL) {
      Splicegraph_genes_gc(&knowngenes);
    }

    if (knowngenes_splices != NULL) {
      /* Splice_gc(&knowngenes_splices); */
      List_free(&knowngenes_splices);
    }

    Splicegraph_free(&splicegraph);

    FREE(crosshyb_extents);
    FREE(primary_extents);
    List_free(&splices);
  }

  Gstruct_free(&gstruct);
  if (knowngenes_gstruct != NULL) {
    IIT_free(&introns_iit);
    IIT_free(&middle_exons_iit);
    IIT_free(&end_exons_iit);
    Gstruct_free(&knowngenes_gstruct);
  }

  return;
}


#if 0
static void
gstruct_solve_twopass (Genome_T genome, IIT_T chromosome_iit, IIT_T knowngenes_iit,
		       char *bamfile, Genomicpos_T max_pairlength, Genomicpos_T shortsplicedist,
		       Genomicpos_T max_exonlength, bool need_unique_p, bool need_primary_p,
		       bool need_canonical_p, int mincount_end_alt, int minsupport) {
  Gstruct_T gstruct, knowngenes_gstruct = NULL;
  IIT_T introns_iit = NULL, middle_exons_iit = NULL, end_exons_iit = NULL;
  Bamreader_T bamreader;

  IIT_T genes_iit;
  char Buffer[20];
  List_T genes, splices, valid_splices, p, q;
  List_T intervallist = NULL, labellist = NULL;

  char *divstring, *typestring, *label;
  List_T divlist = NULL, typelist = NULL;
  Table_T intervaltable, labeltable;
  Tableint_T ngoodhits_low_table, ngoodhits_high_table;

  long int *fwd_extents, *rev_extents, *null_extents;
  long int *tally_matches_low, *tally_mismatches_low,
    *tally_matches_high, *tally_mismatches_high, *primary_extents, *crosshyb_extents;
  Genomicpos_T fwd_chrlength, rev_chrlength, null_chrlength,
    primary_chrlength, crosshyb_chrlength, chrlength;
  char *chr, *copy;
  int quality_counts_match[256], quality_counts_mismatch[256], i;
  int divno, genei;
  Genomicpos_T chroffset;
  Interval_T interval;



  /* First pass */

  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }

  if (knowngenes_iit != NULL) {
    knowngenes_gstruct = Gstruct_knowngenes_input(&introns_iit,&middle_exons_iit,&end_exons_iit,
						  knowngenes_iit,genome,chromosome_iit);
  }

  bamreader = Bamread_new(bamfile);
  gstruct = Gstruct_bam_input(bamreader,/*ngoodhits_low_table*/NULL,/*ngoodhits_high_table*/NULL,
			      /*genes_iit*/NULL,shortsplicedist,max_pairlength,genome,chromosome_iit,
			      desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_unique_p,need_primary_p,
			      /*trust_sam_p*/true,need_canonical_p,
			      /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0);
  Bamread_free(&bamreader);


  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);

  bamreader = Bamread_new(bamfile);
  while ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				 &primary_extents,&crosshyb_extents,
				 &fwd_chrlength,&rev_chrlength,&null_chrlength,
				 &primary_chrlength,&crosshyb_chrlength,gstruct,
				 max_exonlength,mincount_end_alt,minsupport,need_canonical_p,
				 /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0)) != NULL) {
    FREE(null_extents);
    FREE(rev_extents);
    FREE(fwd_extents);


    copy = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(copy,chr);
    divlist = List_push(divlist,(void *) copy);

    if (output_type == OUTPUT_SPLICES) {
      for (p = splices; p != NULL; p = List_next(p)) {
	Splice_print((Splice_T) List_head(p),chr,/*map_iit*/NULL,/*show_invalid_p*/true);
      }
    } else {
      /* Tally */
      if ((divno = IIT_find_one(chromosome_iit,chr)) < 0) {
	fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",chr);
	exit(9);
      } else {
	interval = IIT_interval(chromosome_iit,divno);
	chroffset = Interval_low(interval);
	chrlength = Interval_length(interval);
      }

      Bamread_limit_region(bamreader,chr,/*chrstart*/0U,/*chrend*/chrlength);
      Bamtally_run_lh(&tally_matches_low,&tally_mismatches_low,
		      &tally_matches_high,&tally_mismatches_high,
		      quality_counts_match,quality_counts_mismatch,
		      bamreader,genome,chr,chroffset,/*chrstart*/0U,/*chrend*/chrlength,
		      /*map_iit*/NULL,alloclength,
		      minimum_mapq,good_unique_mapq,maximum_nhits,
		      need_concordant_p,need_unique_p,need_primary_p,
		      blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		      ignore_query_Ns_p,/*verbosep*/true,/*readlevel_p*/false,
		      /*max_softclip*/0);
      Bamread_unlimit_region(bamreader);

      valid_splices = Splice_valid_list(splices);
      genes = Splicegraph_solve_obs(valid_splices,/*knowngenes*/NULL,chr,
				    tally_matches_low,tally_mismatches_low,
				    tally_matches_high,tally_mismatches_high,
				    primary_extents,crosshyb_extents,genome,
				    chroffset,chrlength,bamreader,
				    mincount_end_alt,minsupport,need_canonical_p,
				    min_exonlength,max_exonlength,auto_exonlength,
				    min_intronlength,max_intronlength,/*altpaths_p*/false);
      List_free(&valid_splices);

      labellist = (List_T) NULL;
      for (p = genes, genei = 0; p != NULL; p = List_next(p), genei++) {
	sprintf(Buffer,"Gene%d",genei++);
	label = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
	strcpy(label,Buffer);
	labellist = List_push(labellist,(void *) label);
      }

      Table_put(intervaltable,(void *) copy,genes = List_reverse(genes));
      Table_put(labeltable,(void *) copy,List_reverse(labellist));

      FREE(tally_mismatches_high);
      FREE(tally_matches_high);
      FREE(tally_mismatches_low);
      FREE(tally_matches_low);
    }

    List_free(&splices);
  }

  Bamread_free(&bamreader);

  Gstruct_free(&gstruct);


  /* Make genes_iit */
  divlist = List_reverse(divlist);

  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(divlist,(void *) divstring);

  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(NULL,(void *) typestring);

  genes_iit = IIT_create(divlist,typelist,/*fieldlist*/NULL,intervaltable,
			 labeltable,/*datatable*/NULL,/*divsort*/NO_SORT,/*version*/IIT_LATEST_VERSION,
			 /*presortedp*/false);

  FREE(typestring);
  List_free(&typelist);

  for (q = divlist; q != NULL; q = List_next(q)) {
    divstring = (char *) List_head(q);

    labellist = (List_T) Table_get(labeltable,(void *) divstring);
    for (p = labellist; p != NULL; p = List_next(p)) {
      label = (char *) List_head(p);
      FREE(label);
    }
    List_free(&labellist);

    intervallist = (List_T) Table_get(intervaltable,(void *) divstring);
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      Interval_free(&interval);
    }
    List_free(&intervallist);

    FREE(divstring);
  }
  List_free(&divlist);

  Table_free(&labeltable);
  Table_free(&intervaltable);


  bamreader = Bamread_new(bamfile);
  Gstruct_recount_multimappers(&ngoodhits_low_table,&ngoodhits_high_table,bamreader,
			       genes_iit,maximum_nhits);
  Bamread_free(&bamreader);

  /* Second pass */
  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }

  bamreader = Bamread_new(bamfile);
  gstruct = Gstruct_bam_input(bamreader,ngoodhits_low_table,ngoodhits_high_table,
			      genes_iit,shortsplicedist,max_pairlength,genome,chromosome_iit,
			      desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
			      /*need_unique_p*/false,/*need_primary_p*/false,
			      /*trust_sam_p*/true,need_canonical_p,
			      /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0);
  Bamread_free(&bamreader);
  IIT_free(&genes_iit);


  bamreader = Bamread_new(bamfile);
  while ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				 &primary_extents,&crosshyb_extents,
				 &fwd_chrlength,&rev_chrlength,&null_chrlength,
				 &primary_chrlength,&crosshyb_chrlength,gstruct,
				 max_exonlength,mincount_end_alt,minsupport,need_canonical_p,
				 /*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0)) != NULL) {
    FREE(null_extents);
    FREE(rev_extents);
    FREE(fwd_extents);

    if (output_type == OUTPUT_SPLICES) {
      for (p = splices; p != NULL; p = List_next(p)) {
	Splice_print((Splice_T) List_head(p),chr,/*map_iit*/NULL,/*show_invalid_p*/true);
      }
    } else {
      /* Tally */
      if ((divno = IIT_find_one(chromosome_iit,chr)) < 0) {
	fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",chr);
	exit(9);
      } else {
	interval = IIT_interval(chromosome_iit,divno);
	chroffset = Interval_low(interval);
	chrlength = Interval_length(interval);
      }

      Bamread_limit_region(bamreader,chr,/*chrstart*/0U,/*chrend*/chrlength);
      Bamtally_run_lh(&tally_matches_low,&tally_mismatches_low,
		      &tally_matches_high,&tally_mismatches_high,
		      quality_counts_match,quality_counts_mismatch,
		      bamreader,genome,chr,chroffset,/*chrstart*/0U,/*chrend*/chrlength,
		      /*map_iit*/NULL,alloclength,
		      minimum_mapq,good_unique_mapq,maximum_nhits,
		      need_concordant_p,need_unique_p,need_primary_p,
		      blocksize,quality_score_adj,min_depth,variant_strands,genomic_diff_p,
		      ignore_query_Ns_p,/*verbosep*/true,/*readlevel_p*/false,
		      /*max_softclip*/0);
      Bamread_unlimit_region(bamreader);

      valid_splices = Splice_valid_list(splices);
      Splicegraph_solve_obs(valid_splices,/*knowngenes*/NULL,chr,
			    tally_matches_low,tally_mismatches_low,
			    tally_matches_high,tally_mismatches_high,
			    primary_extents,crosshyb_extents,genome,
			    chroffset,chrlength,bamreader,
			    mincount_end_alt,minsupport,need_canonical_p,
			    min_exonlength,max_exonlength,auto_exonlength,
			    min_intronlength,max_intronlength,/*altpaths_p*/true);
      List_free(&valid_splices);

      FREE(tally_mismatches_high);
      FREE(tally_matches_high);
      FREE(tally_mismatches_low);
      FREE(tally_matches_low);
    }

    /* Splice_gc(&splices); */
    List_free(&splices);
  }

  Bamread_free(&bamreader);

  Gstruct_free(&gstruct);


#if 0
  if ((n = Tableint_length(ngoodhits_low_table)) > 0) {
    keys = (char **) Tableint_keys(ngoodhits_low_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
  }
#endif
  Tableint_free(&ngoodhits_low_table);


#if 0
  if ((n = Tableint_length(ngoodhits_high_table)) > 0) {
    keys = (char **) Tableint_keys(ngoodhits_high_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
  }
#endif
  Tableint_free(&ngoodhits_high_table);

  return;
}

#endif


int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL;
  char *mapdir;
  IIT_T knowngenes_iit = NULL;

  char *iitfile;
  IIT_T chromosome_iit = NULL;
  Genome_T genome = NULL;

  /* For results of parsing -r or --range */
  Genomicpos_T genomicstart, genomiclength;
  Genomicpos_T chroffset = 0U, chrlength;
  bool revcomp;

#ifdef BAM_INPUT
  List_T bamfiles;
#endif

  int i;
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:r:M:m:q:n:C:U:P:A:Me:E:i:1V?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"main-only")) {
	altpaths_p = false;
      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;

      } else if (!strcmp(long_name,"read-group")) {
	desired_read_group = optarg;

      } else if (!strcmp(long_name,"pairmax")) {
	max_pairlength = atoi(optarg);
	shortsplicedist = atoi(optarg);
	alloclength = atoi(optarg);
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'r': user_range = optarg; break;

    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;

    case 'q': minimum_mapq = atoi(optarg); break;
    case 'n': maximum_nhits = atoi(optarg); break;

    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'U':
      switch (atoi(optarg)) {
      case 0: need_unique_p = false; break;
      case 1: need_unique_p = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'P':
      switch (atoi(optarg)) {
      case 0: need_primary_p = false; break;
      case 1: need_primary_p = true; break;
      default: fprintf(stderr,"Primary mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'A':
      if (!strcmp(optarg,"splices")) {
	output_type = OUTPUT_SPLICES;
      } else if (!strcmp(optarg,"paths")) {
	output_type = OUTPUT_PATHS;
      } else if (!strcmp(optarg,"genes")) {
	output_type = OUTPUT_GENES;
      } else if (!strcmp(optarg,"diffs")) {
	output_type = OUTPUT_DIFFS;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
      }
      break;

    case 'e': min_exonlength = atoi(optarg); break;
    case 'E': max_exonlength = atoi(optarg); break;
    case 'i': min_intronlength = atoi(optarg); break;

    case '1': twopassp = false; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;


  if (dbroot != NULL) {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
			/*access*/USE_MMAP_ONLY);
    Maxent_hr_setup(Genome_blocks(genome));

    if (map_iitfile != NULL) {
      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      if ((knowngenes_iit = IIT_read(map_iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
				     /*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Read known genes from %s\n",map_iitfile);
      } else {
	iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+strlen(map_iitfile)+1,sizeof(char));
	sprintf(iitfile,"%s/%s",mapdir,map_iitfile);
	if ((knowngenes_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
				       /*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	  fprintf(stderr,"Read known genes from %s\n",iitfile);
	} else {
	  fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	}
	FREE(iitfile);
      }
      FREE(mapdir);
    }

    if (user_range != NULL) {
      Parserange_universal(&range_chromosome,&revcomp,&genomicstart,&genomiclength,&range_chrstart,&range_chrend,
			   &chroffset,&chrlength,user_range,genomesubdir,fileroot);
    }
  }


  Cappaths_setup();

#if 0
  Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
		       &chroffset,&chrlength,argv[1],genomesubdir,fileroot);
  /* Bamread_limit_region(bamreader,chromosome,chrstart,chrend); */
  gstruct_solve_twopass(genome,chromosome_iit,knowngenes_iit,bamfile,
			max_pairlength,shortsplicedist,max_exonlength,
			need_unique_p,need_primary_p,need_canonical_p,mincount_end_alt,minsupport);
  /* Bamread_unlimit_region(bamreader); */
#endif
#if 0
  /* Handles all chromosomes */
  gstruct_solve_twopass(genome,chromosome_iit,knowngenes_iit,bamfile,
			max_pairlength,shortsplicedist,max_exonlength,
			need_unique_p,need_primary_p,need_canonical_p,mincount_end_alt,minsupport);
#endif

  bamfiles = (List_T) NULL;
  for (i = 0; i < argc; i++) {
    fprintf(stderr,"BAM file: %s\n",argv[i]);
    bamfiles = List_push(bamfiles,(void *) argv[i]);
  }
  bamfiles = List_reverse(bamfiles);

  gstruct_solve_onepass(genome,chromosome_iit,knowngenes_iit,bamfiles,
			max_pairlength,shortsplicedist,max_exonlength,
			need_unique_p,need_primary_p,need_canonical_p,mincount_end_alt,minsupport);

  List_free(&bamfiles);

  if (knowngenes_iit != NULL) {
    IIT_free(&knowngenes_iit);
  }

  Genome_free(&genome);
  IIT_free(&chromosome_iit);

  if (genomesubdir != NULL) {
    FREE(fileroot);
    FREE(dbversion);
    FREE(genomesubdir);
  }

  if (dbroot != NULL) {
    FREE(dbroot);
  }

  fprintf(stderr,"Completed successfully\n");

  return 0;
}
