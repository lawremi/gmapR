static char rcsid[] = "$Id: exon-exon-junctions.c 143393 2014-08-05 16:01:25Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <math.h>		/* For qsort */
#include <ctype.h>		/* For isspace */


#include "bool.h"
#include "mem.h"
#include "genomicpos.h"
#include "genome.h"
#include "iit-read.h"
#include "list.h"
#include "oligoindex.h"
#include "table.h"
#include "tableuint.h"
#include "datadir.h"
#include "complement.h"
#include "translation.h"
#include "bamread.h"
#include "bamtally.h"
#include "getopt.h"

#define MAX_GENENAMES 5
#define MAX_EXONNAMES 5
#define SIMILARITY_WIDTH 25

static char *ignore_filename = NULL;
static Table_T ignore_table = NULL;
static char *intertranscript_outer_iitfile = NULL;
static char *intertranscript_inner_iitfile = NULL;
static bool rename_types_p = true;

static int halflength5 = 75;
static int halflength3 = 75;
static int extra_exon = 1;	/* To capture alternate splice of gene fusion */

static bool all_exons_p = false;
static bool print_proteins_p = false;
static bool inframe_only_p = false;

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;

/* For bamreader */
static int alloclength = 200000;
static char *bamfile = NULL;
static Bamreader_T bamreader = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_concordant_p = false;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;
static bool blockp = true;
static int blocksize = 1000;
static int quality_score_adj = 0;
static int quality_counts_match[256];
static int quality_counts_mismatch[256];


static char *countsfile_name = NULL;


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Extending to other exons */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Intragenic test */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static bool exists = true;	/* For setting a value for acc_used_tables */
static char *bam_lacks_chr = NULL;
static int bam_lacks_chr_length = 0;


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"genome", required_argument, 0, 'd'}, /* dbroot */

  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  {"bamfile", required_argument, 0, 'B'}, /* bamfile */
  {"bam-lacks-chr", required_argument, 0, 'P'}, /* bam_lacks_chr */

  {"ignore", required_argument, 0, 'I'},  /* ignore_filename */
  {"outer", required_argument, 0, 'o'}, /* intertranscript_outer_iitfile */
  {"inner", required_argument, 0, 'i'}, /* intertranscript_inner_iitfile */
  {"dont-rename-types", no_argument, 0, 0},	    /* rename_types_p */

  {"extend5", required_argument, 0, '5'}, /* halflength5 */
  {"extend3", required_argument, 0, '3'}, /* halflength3 */
  {"insertlength", required_argument, 0, 'l'}, /* max_insertlength */
  {"all-exons", no_argument, 0, 'A'},	       /* all_exons_p */
  {"countsfile", required_argument, 0, 'c'},       /* countsfile_name */
  {"all-proteins", no_argument, 0, 0},		   /* print_proteins_p, inframe_only_p */
  {"inframe-proteins", no_argument, 0, 0}, /* print_proteins_p, inframe_only_p */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"exon-exon-junctions\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


typedef struct AA_binding_T *AA_binding_T;
struct AA_binding_T {
  char *aafragment;
  List_T exonnames;
};

static AA_binding_T
AA_binding_new (char *aafragment) {
  AA_binding_T new = (AA_binding_T) MALLOC(sizeof(*new));

  new->aafragment = aafragment;
  new->exonnames = (List_T) NULL;
  return new;
}

static void
AA_binding_free (AA_binding_T *old) {
  FREE((*old)->aafragment);
  List_free(&((*old)->exonnames));
  /* Don't free actual exonnames, which is stored primarily in Site_T */
  FREE(*old);
  return;
}


typedef struct Site_T *Site_T;

struct Site_T {
  char *ntfragment;
  List_T aabindings;
  Genomicpos_T chrpos;
  int sign;
  double prob;
  long int tally;
  List_T accessions;
  List_T exonnames;
  bool origp;			/* if matches the given chrpos, when
				   exact_junction_p is true
				   (orig_source == 'X' or 'H').
				   Otherwise, it is an alternate
				   site.  */
  bool novelp;
  bool lastp;			/* last donor (exon <n-1>/n) */
  bool firstp;			/* first acceptor (exon 2) */
};

static Site_T
Site_new (char *ntfragment, Genomicpos_T chrpos, int sign, double prob, long int tally, bool novelp) {
  Site_T new = (Site_T) MALLOC(sizeof(*new));

  new->ntfragment = ntfragment;
  new->aabindings = (List_T) NULL;
  new->chrpos = chrpos;
  new->sign = sign;
  new->prob = prob;
  new->tally = tally;
  new->accessions = (List_T) NULL;
  new->exonnames = (List_T) NULL;
  new->origp = false;
  new->novelp = novelp;
  new->lastp = false;
  new->firstp = false;
  return new;
}

static void
Site_free (Site_T *old) {
  List_T p;
  AA_binding_T aabinding;
  char *name;

  for (p = (*old)->aabindings; p != NULL; p = List_next(p)) {
    aabinding = (AA_binding_T) List_head(p);
    AA_binding_free(&aabinding);
  }
  List_free(&((*old)->aabindings));

  for (p = (*old)->accessions; p != NULL; p = List_next(p)) {
    name = (char *) List_head(p);
    FREE(name);
  }
  List_free(&((*old)->accessions));

  for (p = (*old)->exonnames; p != NULL; p = List_next(p)) {
    name = (char *) List_head(p);
    FREE(name);
  }
  List_free(&((*old)->exonnames));

  FREE(*old);
  return;
}

static void
Site_aabinding_add (Site_T this, char *aafragment, char *exonname) {
  List_T p;
  AA_binding_T aabinding;

  for (p = this->aabindings; p != NULL; p = List_next(p)) {
    aabinding = (AA_binding_T) List_head(p);
    if (!strcmp(aabinding->aafragment,aafragment)) {
      FREE(aafragment);
      aabinding->exonnames = List_push(aabinding->exonnames,(void *) exonname);
      return;
    }
  }

  aabinding = AA_binding_new(aafragment);
  aabinding->exonnames = List_push(aabinding->exonnames,(void *) exonname);
  this->aabindings = List_push(this->aabindings,(void *) aabinding);
  return;
}

static int
Site_cmp (const void *a, const void *b) {
  Site_T x = * (Site_T *) a;
  Site_T y = * (Site_T *) b;

  if (x->chrpos < y->chrpos) {
    return -1;
  } else if (y->chrpos < x->chrpos) {
    return +1;
  } else {
    return 0;
  }
}


/* Two accessions are in the same gene, even if they don't have the
   same accession name, if they share a donor or acceptor site */
static bool
same_gene_p (char *accession1, int acclength1, char *accession2, int acclength2,
	     Site_T *sites, int nsites) {
  List_T p, q;
  char *exonname1, *exonname2;
  int i;

  for (i = 0; i < nsites; i++) {
    for (p = sites[i]->exonnames; p != NULL; p = List_next(p)) {
      exonname1 = (char *) List_head(p);
      if (!strncmp(exonname1,accession1,acclength1)) {
	for (q = sites[i]->exonnames; q != NULL; q = List_next(q)) {
	  if (p == q) {
	    /* Skip */
	  } else {
	    exonname2 = (char *) List_head(q);
	    if (!strncmp(exonname2,accession2,acclength2)) {
	      debug2(fprintf(stderr,"accessions %s and %s share a site at %u\n",exonname1,exonname2,sites[i]->chrpos));
	      return true;
	    }
	  }
	}
      }
    }
  }

  return false;
}


static void
build_accession_pairs (Table_T accession_pair_table, Site_T *sites, int nsites) {
  List_T p, q;
  char *exonname1, *exonname2;
  char *accession_pair, *a, *b;
  int acclength1, acclength2;
  int i;

  for (i = 0; i < nsites; i++) {
    for (p = sites[i]->exonnames; p != NULL; p = List_next(p)) {
      exonname1 = (char *) List_head(p);
      a = strstr(exonname1,"_exon");
      acclength1 = a - exonname1;

      for (q = sites[i]->exonnames; q != NULL; q = List_next(q)) {
	if (p == q) {
	  /* Skip */
	} else {
	  exonname2 = (char *) List_head(q);
	  b = strstr(exonname2,"_exon");
	  acclength2 = b - exonname2;

	  accession_pair = (char *) CALLOC(acclength1 + strlen("_") + acclength2 + 1,sizeof(char));
	  sprintf(accession_pair,"%.*s_%.*s",acclength1,exonname1,acclength2,exonname2);
	  if (Table_get(accession_pair_table,accession_pair) != NULL) {
	    FREE(accession_pair);
	  } else {
	    Table_put(accession_pair_table,accession_pair,(void *) true);
	    /* printf("%s\n",accession_pair); */
	  }

	  accession_pair = (char *) CALLOC(acclength2 + strlen("_") + acclength1 + 1,sizeof(char));
	  sprintf(accession_pair,"%.*s_%.*s",acclength2,exonname2,acclength1,exonname1);
	  if (Table_get(accession_pair_table,accession_pair) != NULL) {
	    FREE(accession_pair);
	  } else {
	    Table_put(accession_pair_table,accession_pair,(void *) true);
	    /* printf("%s\n",accession_pair); */
	  }

	}
      }
    }
  }

  return;
}


static bool
accession_pair_p (Table_T accession_pair_table, char *donor_exonname, int donor_acclength,
		  char *acceptor_exonname, int acceptor_acclength) {
  void *result;
  char *accession_pair;

  accession_pair = (char *) CALLOC(donor_acclength + strlen("_") + acceptor_acclength + 1,sizeof(char));
  sprintf(accession_pair,"%.*s_%.*s",donor_acclength,donor_exonname,acceptor_acclength,acceptor_exonname);

  result = Table_get(accession_pair_table,accession_pair);
  FREE(accession_pair);
  if (result == NULL) {
    return false;
  } else {
    return true;
  }
}


/* Used to look at accession_pair_table */
static bool
Site_intragenic_pair_p (Site_T donor, Site_T acceptor) {
  List_T p, q;
  char *donor_exonname, *acceptor_exonname, *a, *b;
  int donor_acclength, acceptor_acclength;
  
  for (p = donor->exonnames; p != NULL; p = List_next(p)) {
    donor_exonname = (char *) List_head(p);
    if ((a = strstr(donor_exonname,"_exon")) == NULL) {
      debug2(fprintf(stderr,"Could not find _exon in donor_exonname %s\n",donor_exonname));
    } else {
      for (q = acceptor->exonnames; q != NULL; q = List_next(q)) {
	acceptor_exonname = (char *) List_head(q);
	if ((b = strstr(acceptor_exonname,"_exon")) == NULL) {
	  debug2(fprintf(stderr,"Could not find _exon in acceptor_exonname %s\n",acceptor_exonname));
	} else {
	  donor_acclength = a - donor_exonname;
	  acceptor_acclength = b - acceptor_exonname;
	  debug2(fprintf(stderr,"Comparing %s and %s\n",donor_exonname,acceptor_exonname));
	  if (donor_acclength == acceptor_acclength && !strncmp(donor_exonname,acceptor_exonname,donor_acclength)) {
	    return true;
#if 0
	  } else if (accession_pair_p(accession_pair_table,donor_exonname,donor_acclength,acceptor_exonname,acceptor_acclength) == true) {
	    /*
	    fprintf(stderr,"For a deletion, donor %.*s and acceptor %.*s are the same gene\n",
		    donor_acclength,donor_exonname,acceptor_acclength,acceptor_exonname);
	    */
	    return true;
#endif
	  }
	}
      }
    }
  }

  return false;
}

static void
Site_store_acc (Site_T site, char *desired_acc) {
  List_T p;
  char *copy, *acc;

  for (p = site->accessions; p != NULL; p = List_next(p)) {
    acc = (char *) List_head(p);
    if (!strcmp(acc,desired_acc)) {
      return;
    }
  }
  copy = (char *) CALLOC(strlen(desired_acc)+1,sizeof(char));
  strcpy(copy,desired_acc);
  site->accessions = List_push(site->accessions,(void *) copy);
  return;
}


static void
Site_store_accs (Site_T site, Table_T used_table) {
  List_T p;
  char *acc, *copy;
  
  /* fprintf(stderr,"Storing accs from site at %u\n",site->chrpos); */
  for (p = site->accessions; p != NULL; p = List_next(p)) {
    acc = (char *) List_head(p);
    if (Table_get(used_table,acc) == NULL) {
      copy = (char *) CALLOC(strlen(acc)+1,sizeof(char));
      strcpy(copy,acc);
      Table_put(used_table,copy,(void *) &exists);
    }
  }

  return;
}


#if 0
static void
Site_store_accs_old (Site_T site, Table_T used_table) {
  List_T p;
  char *exonname, *a;
  int acclength;
  char *acc;
  
  for (p = site->exonnames; p != NULL; p = List_next(p)) {
    exonname = (char *) List_head(p);
    if (exonname[0] == '*') {
      /* Skip novel exon, which may not correspond to the given acc */
    } else if ((a = strstr(exonname,"_exon")) == NULL) {
      debug2(fprintf(stderr,"Could not find _exon in donor_exonname %s\n",exonname));
    } else {
      acclength = a - exonname;
      acc = (char *) CALLOC(acclength+1,sizeof(char));
      strncpy(acc,exonname,acclength);
      if (Table_get(used_table,acc) != NULL) {
	FREE(acc);
      } else {
	Table_put(used_table,acc,(void *) &exists);
      }
    }
  }

  return;
}
#endif


typedef struct Exon_T *Exon_T;

struct Exon_T {
  Genomicpos_T exonlow;
  Genomicpos_T exonhigh;
  Genomicpos_T exonstart;
  Genomicpos_T exonend;
  Genomicpos_T length;
  int sign;
  char *ntstring;
  char *aastring;
  long int *tally_matches;
  long int *tally_mismatches;
};


static char complCode[128] = COMPLEMENT_LC;


static void
make_complement_inplace (char *sequence, Genomicpos_T length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


static Exon_T
Exon_new (Genomicpos_T exonlow, Genomicpos_T exonhigh, char *chr, Genomicpos_T chroffset,
	  int sign, Genome_T genome) {
  Exon_T new = (Exon_T) MALLOC(sizeof(*new));
  List_T intervallist = NULL, labellist = NULL, datalist = NULL;
  char *chrptr;


  new->exonlow = exonlow;
  new->exonhigh = exonhigh;
  new->length = exonhigh - exonlow + 1;
  new->ntstring = (char *) CALLOC(new->length + 1,sizeof(char));
  new->aastring = (char *) CALLOC(new->length + 1,sizeof(char));
  new->sign = sign;

  if (sign > 0) {
    new->exonstart = exonlow;
    new->exonend = exonhigh;
    Genome_fill_buffer_simple(genome,chroffset+exonlow-1U,new->length,new->ntstring);
  } else {
    new->exonstart = exonhigh;
    new->exonend = exonlow;
    Genome_fill_buffer_simple(genome,chroffset+exonlow-1U,new->length,new->ntstring);
    make_complement_inplace(new->ntstring,new->length);
  }

  if (bamreader == NULL) {
    new->tally_matches = (long int *) NULL;
    new->tally_mismatches = (long int *) NULL;
  } else {
    if (bam_lacks_chr == NULL) {
      chrptr = chr;
    } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
      chrptr = &(chr[bam_lacks_chr_length]);
    } else {
      chrptr = chr;
    }

    Bamread_limit_region(bamreader,chrptr,/*chrstart*/exonlow,/*chrend*/exonhigh);
    Bamtally_run(&new->tally_matches,&new->tally_mismatches,
		 &intervallist,&labellist,&datalist,
		 quality_counts_match,quality_counts_mismatch,
		 bamreader,genome,/*printchr*/chr,chroffset,/*chrstart*/exonlow,/*chrend*/exonhigh,
		 /*map_iit*/NULL,alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		 /*desired_read_group*/NULL,minimum_mapq,good_unique_mapq,maximum_nhits,
		 need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		 /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		 /*output_type*/OUTPUT_TALLY,blockp,blocksize,
		 quality_score_adj,/*min_depth*/1,/*variant_strands*/0,
		 /*genomic_diff_p*/false,/*signed_counts_p*/false,/*ignore_query_Ns_p*/true,
		 /*print_indels_p*/false,/*print_totals_p*/false,
		 /*print_cycles_p*/false,/*print_quality_scores_p*/false,
		 /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		 /*verbosep*/false,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		 /*bamfile*/NULL);
    Bamread_unlimit_region(bamreader);
  }

  debug(printf("Exon %u..%u (%d): %s\n",exonlow,exonhigh,new->length,new->ntstring));

  return new;
}

static void
Exon_free (Exon_T *old) {
  if ((*old)->tally_matches != NULL) {
    FREE((*old)->tally_matches);
    FREE((*old)->tally_mismatches);
  }
  FREE((*old)->aastring);
  FREE((*old)->ntstring);
  FREE(*old);
  return;
}


static int
string_cmp (const void *a, const void *b) {
  char *x = * (char **) a;
  char *y = * (char **) b;

  return strcmp(x,y);
}


static List_T
list_sort (List_T names) {
  List_T sorted = NULL;
  char **array;
  int n, i;

  n = List_length(names);
  array = (char **) List_to_array(names,NULL);
  qsort(array,n,sizeof(char *),string_cmp);
  for (i = n-1; i >= 0; i--) {
    sorted = List_push(sorted,(void *) array[i]);
  }
  FREE(array);
  List_free(&names);
  return sorted;
}

static List_T
list_sort_uniq (List_T names) {
  List_T sorted = NULL;
  char **array;
  int n, i, lasti;

  n = List_length(names);
  array = (char **) List_to_array(names,NULL);
  qsort(array,n,sizeof(char *),string_cmp);
  sorted = List_push(sorted,(void *) array[n-1]);
  lasti = n-1;
  for (i = n-2; i >= 0; i--) {
    if (strcmp(array[i],array[lasti]) != 0) {
      sorted = List_push(sorted,(void *) array[i]);
      lasti = i;
    } else {
      FREE(array[i]);
    }
  }
  FREE(array);
  List_free(&names);
  return sorted;
}


static char *
concatenate_exons (Genomicpos_T *length, Exon_T *exons, int nexons) {
  char *gene_sequence;
  int i;

  *length = 0;
  for (i = 0; i < nexons; i++) {
    *length += exons[i]->length;
  }

  gene_sequence = MALLOC((*length + 1) * sizeof(char));
  *length = 0;
  for (i = 0; i < nexons; i++) {
    strcpy(&(gene_sequence[*length]),exons[i]->ntstring);
    *length += exons[i]->length;
  }
  gene_sequence[*length] = '\0';
  return gene_sequence;
}

static void
add_aa_to_exons (char *aa_sequence, Exon_T *exons, int nexons) {
  Genomicpos_T length;
  int i;

  length = 0;
  for (i = 0; i < nexons; i++) {
    strncpy(exons[i]->aastring,&(aa_sequence[length]),exons[i]->length);
    length += exons[i]->length;
  }
  return;
}


static void
print_gene (char *gene_sequence, char *acc) {
  printf(">Gene-%s\n",acc);
  printf("%s\n",gene_sequence);
  
  return;
}


/* Assumes gene name is first word in comment */
static char *
get_genename (char *chr, Genomicpos_T chrpos, int sign, IIT_T map_iit) {
  List_T names = NULL, l;
  char *genename, *annot, *restofheader, *acc, *p, *q;
  int *matches, nmatches, index, gene_namelength, i, k;
  bool allocp;
  

  matches = IIT_get_signed(&nmatches,map_iit,chr,chrpos,chrpos,sign,/*sortp*/false);
  for (i = 0; i < nmatches && i < MAX_GENENAMES; i++) {
    index = matches[i];
    acc = IIT_label(map_iit,index,&allocp);
    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
      /* Skip readthrough acc */
      if (allocp) FREE(acc);
    } else {
      if (allocp) FREE(acc);
      annot = IIT_annotation(&restofheader,map_iit,index,&allocp);

      /* Skip header */
      p = annot;
      while (*p != '\0' && *p != '\n' && !isspace(*p)) {
	p++;
      }

      gene_namelength = p - annot;
      genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
      strncpy(genename,annot,gene_namelength);
      for (p = genename, k = 0; k < gene_namelength; p++, k++) {
	if (*p == '-') {
	  *p = '.';
	}
      }

      names = List_push(names,(void *) genename);
    }
  }
  FREE(matches);

  if (names == NULL) {
    genename = (char *) CALLOC(strlen("Novel")+1,sizeof(char));
    strcpy(genename,"Novel");

  } else {
    names = list_sort_uniq(names);
    gene_namelength = 0;
    for (l = names; l != NULL; l = List_next(l)) {
      q = (char *) List_head(l);
      gene_namelength += strlen(q) + 1;
    }

    p = genename = (char *) CALLOC(gene_namelength,sizeof(char));
    for (l = names; l != NULL; l = List_next(l)) {
      q = (char *) List_head(l);
      strcpy(p,q);
      p += strlen(q);
      *p++ = ',';
      FREE(q);
    }
    p--;
    *p = '\0';

    List_free(&names);
  }

  return genename;
}


static Exon_T *
get_exons (int *nexons, int *margin5, int *margin3, char *annot, char *chr,
	   Genomicpos_T chroffset, Genome_T genome, int halflength) {
  Exon_T *array;
  List_T exons;
  Genomicpos_T exonlow, exonhigh, exonstart, exonend;
  int exonlength;
  int sign;
  char *p;
  bool firstp = true, lastp = false;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;

  exons = (List_T) NULL;
  *margin5 = *margin3 = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {

      /* Advance to next exon to see if this is the last one */
      while (*p != '\0' && *p != '\n') p++;
      if (*p == '\n') p++;
      if (*p == '\0') {
	debug(printf("This is the last exon\n"));
	lastp = true;
      }

      if (exonstart <= exonend) {
	sign = +1;
	exonlow = exonstart;
	exonhigh = exonend;
	exonlength = (int) (exonhigh - exonlow + 1);
	if (exonlength >= halflength) {
	  /* Do nothing */
	} else if (firstp == true) {
	  exonlow = exonend - (halflength - 1);
	  *margin5 = exonstart - exonlow;
	} else if (lastp == true) {
	  exonhigh = exonstart + (halflength - 1);
	  *margin3 = exonhigh - exonend;
	}

      } else {
	sign = -1;
	exonlow = exonend;
	exonhigh = exonstart;
	exonlength = (int) (exonhigh - exonlow + 1);
	if (exonlength >= halflength) {
	  /* Do nothing */
	} else if (firstp == true) {
	  exonhigh = exonend + (halflength - 1);
	  *margin5 = exonhigh - exonstart;
	} else if (lastp == true) {
	  exonlow = exonstart - (halflength - 1);
	  *margin3 = exonend - exonlow;
	}
      }

      firstp = false;
      exons = List_push(exons,(void *) Exon_new(exonlow,exonhigh,chr,chroffset,sign,genome));
    }

  }

  exons = List_reverse(exons);
  *nexons = List_length(exons);
  array = (Exon_T *) List_to_array(exons,NULL);
  List_free(&exons);

  return array;
}


static char *
get_fragment_left (char *ntptr, char *aaptr, long int *tally, int exoni, Exon_T *exons, int remainder) {
  int length, i;

  if (exoni < 0) {
    return (char *) NULL;
  } else if ((length = exons[exoni]->length) >= remainder) {
    ntptr -= remainder;
    aaptr -= remainder;
    strncpy(ntptr,&(exons[exoni]->ntstring[length - remainder]),remainder);
    strncpy(aaptr,&(exons[exoni]->aastring[length - remainder]),remainder);

    if (exons[exoni]->tally_matches != NULL) {
      if (exons[exoni]->sign > 0) {
	for (i = length - remainder; i < length; i++) {
	  *tally += exons[exoni]->tally_matches[i];
	  *tally += exons[exoni]->tally_mismatches[i];
	}
      } else {
	for (i = 0; i < remainder; i++) {
	  *tally += exons[exoni]->tally_matches[i];
	  *tally += exons[exoni]->tally_mismatches[i];
	}
      }
    }
    return ntptr;

  } else {
    ntptr -= length;
    aaptr -= length;
    strncpy(ntptr,exons[exoni]->ntstring,length);
    strncpy(aaptr,exons[exoni]->aastring,length);

    if (exons[exoni]->tally_matches != NULL) {
      for (i = 0; i < length; i++) {
	*tally += exons[exoni]->tally_matches[i];
	*tally += exons[exoni]->tally_mismatches[i];
      }
    }
    return get_fragment_left(ntptr,aaptr,&(*tally),exoni-1,exons,remainder - length);
  }
}


static bool
get_donor_bounds_plus (int *firstexoni, int *lastexoni, Exon_T *exons, int nexons,
		       Genomicpos_T chrpos, int max_insertlength, Genomicpos_T overhang,
		       bool exact_junction_p) {
  int length;
  int exoni;

  debug1(printf("get_donor_bounds_plus called with %d exons and chrpos %u",
		nexons,chrpos));
  if (chrpos < overhang) {
    chrpos = 0;
  } else {
    chrpos -= overhang;
  }
  debug1(printf(" (%u after overhang)\n",chrpos));


  exoni = 0;			/* was -1 */
  while (exoni < nexons - 1 && exons[exoni]->exonend < chrpos) {
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    exoni++;
  }
  if (exons[exoni]->exonend < chrpos) {
    debug1(printf("  went to last exon\n"));
    return false;
  } else if (exact_junction_p == true && exons[exoni]->exonend != chrpos) {
    debug1(printf("  not exact junction\n"));
    return false;
  } else {
    *firstexoni = exoni;
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    debug1(printf("  first exon ending right of chrpos is %d\n\n",*firstexoni));
  }

  length = 0;
  while (exoni + 1 < nexons - 1 && length < max_insertlength) {
    exoni++;
    length += exons[exoni]->length;
    debug1(printf("  exon %d, %u..%u, length %d\n",
		  exoni,exons[exoni]->exonstart,exons[exoni]->exonend,length));
  }
  *lastexoni = exoni;
  debug1(printf("  last exon within insert length of %d is %d\n",
		max_insertlength,*lastexoni));

  if (*firstexoni - extra_exon >= 0) {
    debug1(printf(" extra_exon of %d allows first exon to change from %d to %d\n",
		  extra_exon,*firstexoni,*firstexoni - extra_exon));
    *firstexoni -= extra_exon;
  }
  if (*lastexoni + extra_exon < nexons - 1) {
    debug1(printf(" extra_exon of %d allows last exon to change from %d to %d\n",
		  extra_exon,*lastexoni,*lastexoni + extra_exon));
    *lastexoni += extra_exon;
  }

  return true;
}

static bool
get_donor_bounds_minus (int *firstexoni, int *lastexoni, Exon_T *exons, int nexons,
			Genomicpos_T chrpos, int max_insertlength, Genomicpos_T overhang,
			bool exact_junction_p) {
  int length;
  int exoni;

  debug1(printf("get_donor_bounds_minus called with %d exons and chrpos %u",
		nexons,chrpos));
  chrpos += overhang;
  debug1(printf(" (%u after overhang)\n",chrpos));


  exoni = 0;			/* was -1 */
  while (exoni < nexons - 1 && exons[exoni]->exonend > chrpos) {
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    exoni++;
  }
  if (exons[exoni]->exonend > chrpos) {
    debug1(printf("  went to last exon\n"));
    return false;
  } else if (exact_junction_p == true && exons[exoni]->exonend != chrpos) {
    debug1(printf("  not exact junction\n"));
    return false;
  } else {
    *firstexoni = exoni;
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    debug1(printf("  first exon ending right of chrpos is %d\n\n",*firstexoni));
  }

  length = 0;
  while (exoni + 1 < nexons - 1 && length < max_insertlength) {
    exoni++;
    length += exons[exoni]->length;
    debug1(printf("  exon %d, %u..%u, length %d\n",
		  exoni,exons[exoni]->exonstart,exons[exoni]->exonend,length));
  }
  *lastexoni = exoni;
  debug1(printf("  last exon within insert length of %d is %d\n",
		max_insertlength,*lastexoni));

  if (*firstexoni - extra_exon >= 0) {
    debug1(printf(" extra_exon of %d allows first exon to change from %d to %d\n",
		  extra_exon,*firstexoni,*firstexoni - extra_exon));
    *firstexoni -= extra_exon;
  }
  if (*lastexoni + extra_exon < nexons - 1) {
    debug1(printf(" extra_exon of %d allows last exon to change from %d to %d\n",
		  extra_exon,*lastexoni,*lastexoni + extra_exon));
    *lastexoni += extra_exon;
  }

  return true;
}


static char *
get_fragment_right (char *ntptr, char *aaptr, long int *tally, int exoni, int nexons, Exon_T *exons, int remainder) {
  int length, i;

  if (exoni >= nexons) {
    return (char *) NULL;
  } else if ((length = exons[exoni]->length) >= remainder) {
    strncpy(ntptr,exons[exoni]->ntstring,remainder);
    strncpy(aaptr,exons[exoni]->aastring,remainder);

    if (exons[exoni]->tally_matches != NULL) {
      if (exons[exoni]->sign > 0) {
	for (i = 0; i < remainder; i++) {
	  *tally += exons[exoni]->tally_matches[i];
	  *tally += exons[exoni]->tally_mismatches[i];
	}
      } else {
	for (i = length - remainder; i < length; i++) {
	  *tally += exons[exoni]->tally_matches[i];
	  *tally += exons[exoni]->tally_mismatches[i];
	}
      }
    }
    return ntptr + remainder;

  } else {
    strncpy(ntptr,exons[exoni]->ntstring,length);
    strncpy(aaptr,exons[exoni]->aastring,length);

    if (exons[exoni]->tally_matches != NULL) {
      for (i = 0; i < length; i++) {
	*tally += exons[exoni]->tally_matches[i];
	*tally += exons[exoni]->tally_mismatches[i];
      }
    }

    return get_fragment_right(ntptr + length,aaptr + length,&(*tally),exoni+1,nexons,exons,remainder - length);
  }
}

static bool
get_acceptor_bounds_plus (int *firstexoni, int *lastexoni, Exon_T *exons, int nexons,
			  Genomicpos_T chrpos, int max_insertlength, Genomicpos_T overhang,
			  bool exact_junction_p) {
  int length;
  int exoni;

  debug1(printf("get_acceptor_bounds_plus called with %d exons and chrpos %u",
		nexons,chrpos));
  chrpos += overhang;
  debug1(printf(" (%u after overhang)\n",chrpos));


  exoni = nexons - 1;		/* was nexons */
  while (exoni >= 1 && exons[exoni]->exonstart > chrpos) {
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    exoni--;
  }
  if (exons[exoni]->exonstart > chrpos) {
    debug1(printf("  went to first exon\n"));
    return false;
  } else if (exact_junction_p == true && exons[exoni]->exonstart != chrpos) {
    debug1(printf("  not exact junction\n"));
    return false;
  } else {
    *firstexoni = exoni;
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    debug1(printf("  first exon ending left of chrpos is %d\n\n",*firstexoni));
  }

  length = 0;
  while (exoni - 1 >= 1 && length < max_insertlength) {
    exoni--;
    length += exons[exoni]->length;
    debug1(printf("  exon %d, %u..%u, length %d\n",
		  exoni,exons[exoni]->exonstart,exons[exoni]->exonend,length));
  }
  *lastexoni = exoni;
  debug1(printf("  last exon within insert length of %d is %d\n",
		max_insertlength,*lastexoni));

  if (*firstexoni + extra_exon <= nexons - 1) {
    debug1(printf(" extra_exon of %d allows first exon to change from %d to %d\n",
		  extra_exon,*firstexoni,*firstexoni + extra_exon));
    *firstexoni += extra_exon;
  }
  if (*lastexoni - extra_exon > 0) {
    debug1(printf(" extra_exon of %d allows last exon to change from %d to %d\n",
		  extra_exon,*lastexoni,*lastexoni - extra_exon));
    *lastexoni -= extra_exon;
  }

  return true;
}

static bool
get_acceptor_bounds_minus (int *firstexoni, int *lastexoni, Exon_T *exons, int nexons,
			   Genomicpos_T chrpos, int max_insertlength, Genomicpos_T overhang,
			   bool exact_junction_p) {

  int length;
  int exoni;

  debug1(printf("get_acceptor_bounds_minus called with %d exons and chrpos %u",
		nexons,chrpos));
  if (chrpos < overhang) {
    chrpos = 0;
  } else {
    chrpos -= overhang;
  }
  debug1(printf(" (%u after overhang)\n",chrpos));


  exoni = nexons - 1;		/* was nexons */
  while (exoni >= 1 && exons[exoni]->exonstart < chrpos) {
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    exoni--;
  }
  if (exons[exoni]->exonstart < chrpos) {
    debug1(printf("  went to first exon\n"));
    return false;
  } else if (exact_junction_p == true && exons[exoni]->exonstart != chrpos) {
    debug1(printf("  not exact junction\n"));
    return false;
  } else {
    *firstexoni = exoni;
    debug1(printf("  exon %d, %u..%u\n",exoni,exons[exoni]->exonstart,exons[exoni]->exonend));
    debug1(printf("  first exon starting right of chrpos is %d\n\n",*firstexoni));
  }

  length = 0;
  while (exoni - 1 >= 1 && length < max_insertlength) {
    exoni--;
    length += exons[exoni]->length;
    debug1(printf("  exon %d, %u..%u, length %d\n",
		  exoni,exons[exoni]->exonstart,exons[exoni]->exonend,length));
  }
  *lastexoni = exoni;
  debug1(printf("  last exon within insert length of %d is %d\n",
		max_insertlength,*lastexoni));

  if (*firstexoni + extra_exon <= nexons - 1) {
    debug1(printf(" extra_exon of %d allows first exon to change from %d to %d\n",
		  extra_exon,*firstexoni,*firstexoni + extra_exon));
    *firstexoni += extra_exon;
  }
  if (*lastexoni - extra_exon > 0) {
    debug1(printf(" extra_exon of %d allows last exon to change from %d to %d\n",
		  extra_exon,*lastexoni,*lastexoni - extra_exon));
    *lastexoni -= extra_exon;
  }

  return true;
}



static int
numlen (int number) {
  if (number < 10) {
    return 1;
  } else if (number < 100) {
    return 2;
  } else if (number < 1000) {
    return 3;
  } else if (number < 10000) {
    return 4;
  } else if (number < 100000) {
    return 5;
  } else {
    abort();
  }
}


/* exonnames taken from input file */
static void
add_novel_donor (Table_T donor_fragment_table, char sign, char *chr, Genomicpos_T chrpos, double prob,
		 IIT_T map_iit, IIT_T chromosome_iit, Genome_T genome, char *exons) {
  char *ntfragment, *exonname, *acc, *copy, *p, *q, *a, *b;
  int acclength;
  char *chrptr;
  Site_T donor;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;

  int i;
  long int tally, *tally_matches, *tally_mismatches;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL;

  debug(fprintf(stderr,"Adding novel donor\n"));
  chrnum = IIT_find_one(chromosome_iit,chr);
  chroffset = IIT_interval_low(chromosome_iit,chrnum);

  ntfragment = (char *) CALLOC(halflength5+1,sizeof(char));
  ntfragment[halflength5] = '\0';
  if (sign == '+') {
    Genome_fill_buffer_simple(genome,chroffset+chrpos-halflength5,/*length*/halflength5,ntfragment);
    if (bamreader == NULL) {
      tally = 0;
    } else {
      if (bam_lacks_chr == NULL) {
	chrptr = chr;
      } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chr[bam_lacks_chr_length]);
      } else {
	chrptr = chr;
      }

      Bamread_limit_region(bamreader,chrptr,/*chrstart*/chrpos-halflength5,/*chrend*/chrpos);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,/*printchr*/chr,chroffset,/*chrstart*/chrpos-halflength5,/*chrend*/chrpos,
		   /*map_iit*/NULL,alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   /*desired_read_group*/NULL,minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_TALLY,blockp,blocksize,
		   quality_score_adj,/*min_depth*/1,/*variant_strands*/0,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,/*ignore_query_Ns_p*/true,
		   /*print_indels_p*/false,/*print_totals_p*/false,
		   /*print_cycles_p*/false,/*print_quality_scores_p*/false,
		   /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		   /*verbosep*/false,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		   /*bamfile*/NULL);
      tally = 0;
      for (i = 0; i < halflength5; i++) {
	tally += tally_matches[i] + tally_mismatches[i];
      }
      FREE(tally_mismatches);
      FREE(tally_matches);
      Bamread_unlimit_region(bamreader);
    }

    donor = Site_new(ntfragment,/*splicesite*/chrpos,/*sign*/+1,prob,tally,/*novelp*/true);
    /* No aafragment */

  } else if (sign == '-') {
    Genome_fill_buffer_simple(genome,chroffset+chrpos-1,/*length*/halflength5,ntfragment);
    make_complement_inplace(ntfragment,halflength5);
    if (bamreader == NULL) {
      tally = 0;
    } else {
      if (bam_lacks_chr == NULL) {
	chrptr = chr;
      } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chr[bam_lacks_chr_length]);
      } else {
	chrptr = chr;
      }

      Bamread_limit_region(bamreader,chrptr,/*chrstart*/chrpos,/*chrend*/chrpos+halflength5);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,/*printchr*/chr,chroffset,/*chrstart*/chrpos-halflength5,/*chrend*/chrpos,
		   /*map_iit*/NULL,alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   /*desired_read_group*/NULL,minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_TALLY,blockp,blocksize,
		   quality_score_adj,/*min_depth*/1,/*variant_strands*/0,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,/*ignore_query_Ns_p*/true,
		   /*print_indels_p*/false,/*print_totals_p*/false,
		   /*print_cycles_p*/false,/*print_quality_scores_p*/false,
		   /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		   /*verbosep*/false,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		   /*bamfile*/NULL);
      tally = 0;
      for (i = 0; i < halflength5; i++) {
	tally += tally_matches[i] + tally_mismatches[i];
      }
      FREE(tally_mismatches);
      FREE(tally_matches);
      Bamread_unlimit_region(bamreader);
    }

    donor = Site_new(ntfragment,/*splicesite*/chrpos,/*sign*/-1,prob,tally,/*novelp*/true);
    /* No aafragment */

  } else {
    fprintf(stderr,"Sign %c not recognized\n",sign);
    abort();
  }

  donor->origp = true;		/* Novel site was found originally by alignment */
  Table_put(donor_fragment_table,ntfragment,donor);

  if ((p = exons) == NULL) {
    exonname = (char *) CALLOC(strlen("NovelExon") + 1,sizeof(char));
    sprintf(exonname,"NovelExon");
    donor->exonnames = List_push(donor->exonnames,(void *) exonname);
  } else {
    acc = p;
    while (*p != '\0') {
      if (*p == ',') {
	exonname = (char *) CALLOC(p - acc + 1,sizeof(char));
	strncpy(exonname,acc,p-acc);
	q = exonname;

	a = exonname;
	while (*a == '*') {
	  a++;
	}
	b = strstr(exonname,"_exon");
	acclength = b - a;
	copy = (char *) CALLOC(acclength+1,sizeof(char));
	strncpy(copy,a,acclength);
	Site_store_acc(donor,copy);
	FREE(copy);
      
	while (*q != '\0') {
	  if (*q == '-') {
	    *q = '_';
	  }
	  q++;
	}
	donor->exonnames = List_push(donor->exonnames,(void *) exonname);
	p++;
	acc = p;
      } else {
	p++;
      }
    }

    exonname = (char *) CALLOC(p - acc + 1,sizeof(char));
    strncpy(exonname,acc,p-acc);
    q = exonname;

    a = exonname;
    while (*a == '*') {
      a++;
    }
    b = strstr(exonname,"_exon");
    acclength = b - a;
    copy = (char *) CALLOC(acclength+1,sizeof(char));
    strncpy(copy,a,acclength);
    Site_store_acc(donor,copy);
    FREE(copy);
    
    while (*q != '\0') {
      if (*q == '-') {
	*q = '_';
      }
      q++;
    }
    donor->exonnames = List_push(donor->exonnames,(void *) exonname);
  }

  return;
}


static char *
process_donors (bool *found_exact_p, Genomicpos_T *minlow, Genomicpos_T *maxhigh,
		Table_T donor_fragment_table, Table_T donor_continuation_table,
		char *acc, char *chr, Genomicpos_T chrpos, double prob,
		IIT_T map_iit, IIT_T chromosome_iit, Genome_T genome,
		int max_insertlength, Genomicpos_T overhang,
		bool exact_junction_p) {
  char *gene_sequence = NULL, *aa_sequence;
  int margin5, margin3;
  int translation_leftpos, translation_rightpos, translation_length;
  int genestrand;
  char *annot, *restofheader, *divstring;
  bool allocp;
  int index;

  Chrnum_T chrnum;
  Genomicpos_T chroffset, low, high;
  Genomicpos_T genelength;

  Exon_T *exons;
  int nexons, firstexoni, lastexoni, exoni;
  char *ntfragment, *aafragment;
  long int tally;

  Site_T donor, donor_cont;
  Genomicpos_T splicesite;
  char *exonname;


  if (map_iit == NULL) {
    index = -1;
  } else {
    index = IIT_find_one(map_iit,acc);
  }
  if (index < 0) {
    fprintf(stderr,"Cannot find gene for %s\n",acc);
  } else {
    genestrand = IIT_interval_sign(map_iit,index);
    divstring = IIT_divstring_from_index(map_iit,index);
    annot = IIT_annotation(&restofheader,map_iit,index,&allocp);
    chrnum = IIT_find_one(chromosome_iit,chr);
    chroffset = IIT_interval_low(chromosome_iit,chrnum);
    exons = get_exons(&nexons,&margin5,&margin3,annot,/*chr*/divstring,chroffset,genome,halflength5);
    if (allocp) {
      FREE(restofheader);
    }

    gene_sequence = concatenate_exons(&genelength,exons,nexons);
    aa_sequence = Translation_via_genomic(&translation_leftpos,&translation_rightpos,&translation_length,
					  gene_sequence,/*startpos*/margin5,/*endpos*/genelength-margin3,genelength,
					  /*backwardp*/false,/*revcompp*/false,/*fulllengthp*/true,/*cds_startpos*/0);
    add_aa_to_exons(aa_sequence,exons,nexons);
    FREE(aa_sequence);

    /* donor firstexoni goes upward to lastexoni */
    if (nexons == 0) {
      fprintf(stderr,"Gene %s has no exons\n",acc);
      /* Skip */

    } else if (genestrand > 0) {
      if (get_donor_bounds_plus(&firstexoni,&lastexoni,exons,nexons,chrpos,
				max_insertlength,overhang,exact_junction_p) == true) {
	if ((low = exons[firstexoni]->exonlow) < *minlow) {
	  *minlow = low;
	}
	if ((high = exons[lastexoni]->exonhigh) > *maxhigh) {
	  *maxhigh = high;
	}
      }

    } else if (genestrand < 0) {
      if (get_donor_bounds_minus(&firstexoni,&lastexoni,exons,nexons,chrpos,
				 max_insertlength,overhang,exact_junction_p) == true) {
	if ((high = exons[firstexoni]->exonhigh) > *maxhigh) {
	  *maxhigh = high;
	}
	if ((low = exons[lastexoni]->exonlow) < *minlow) {
	  *minlow = low;
	}
      }

    } else {
      fprintf(stderr,"genestrand is 0\n");
      abort();
    }

    debug(printf("For donor acc %s (strand %d), exons %d..%d\n",acc,genestrand,firstexoni,lastexoni));
    /* Was (exoni = firstexoni; exoni <= lastexoni && exoni < nexons - 1; exoni++) */
    for (exoni = 0; exoni < nexons - 1; exoni++) {
      ntfragment = (char *) CALLOC(halflength5+1,sizeof(char));
      aafragment = (char *) CALLOC(halflength5+1,sizeof(char));
      ntfragment[halflength5] = aafragment[halflength5] = '\0';
      tally = 0;
      if (get_fragment_left(&(ntfragment[halflength5]),&(aafragment[halflength5]),&tally,exoni,exons,halflength5) != NULL) {
	debug(printf("left fragment %d/%d: %s\n",exoni,nexons,ntfragment));
	if ((donor = Table_get(donor_fragment_table,ntfragment)) != NULL) {
	  FREE(ntfragment);
	} else {
	  splicesite = exons[exoni]->exonend;
	  donor = Site_new(ntfragment,splicesite,exons[exoni]->sign,prob,tally,/*novelp*/false);
	  Table_put(donor_fragment_table,ntfragment,donor);
	}
	
	Site_store_acc(donor,acc);
	exonname = (char *) CALLOC(strlen(acc) + strlen("_exon") + 
				   numlen(exoni+1) + strlen("/") + numlen(nexons) + 1,sizeof(char));
	sprintf(exonname,"%s_exon%d/%d",acc,exoni+1,nexons);
	donor->exonnames = List_push(donor->exonnames,(void *) exonname);
	if (exoni+1 == nexons - 1) {
	  donor->lastp = true;
	}

	Site_aabinding_add(donor,aafragment,exonname);
      }
      if (exons[exoni]->exonend == chrpos) {
	donor->origp = true;
	*found_exact_p = true;
      }
    }

    /* Do continuation */
    /* Was (exoni = firstexoni + 1; exoni <= lastexoni + 1 && exoni < nexons; exoni++) */
    for (exoni = 1; exoni < nexons; exoni++) {
      ntfragment = (char *) CALLOC(halflength5+1,sizeof(char));
      aafragment = (char *) CALLOC(halflength5+1,sizeof(char));
      ntfragment[halflength5] = aafragment[halflength5] = '\0';
      tally = 0;
      if (get_fragment_right(&(ntfragment[0]),&(aafragment[0]),&tally,exoni,nexons,exons,halflength5) != NULL) {
	debug(printf("right continuation fragment %d/%d: %s\n",exoni,nexons,ntfragment));
	if ((donor_cont = Table_get(donor_continuation_table,ntfragment)) != NULL) {
	  FREE(ntfragment);
	} else {
	  splicesite = exons[exoni]->exonstart;
	  donor_cont = Site_new(ntfragment,splicesite,exons[exoni]->sign,prob,tally,/*novelp*/false);
	  Table_put(donor_continuation_table,ntfragment,donor_cont);
	}
	
	/* Site_store_acc(donor_cont,acc); */
	exonname = (char *) CALLOC(strlen(acc) + strlen("_exon") + 
				   numlen(exoni+1) + strlen("/") + numlen(nexons) + 1,sizeof(char));
	sprintf(exonname,"%s_exon%d/%d",acc,exoni+1,nexons);
	donor_cont->exonnames = List_push(donor_cont->exonnames,(void *) exonname);
	if (exoni+1 == nexons - 1) {
	  donor_cont->lastp = true;
	}

	Site_aabinding_add(donor_cont,aafragment,exonname);
      }
    }

    for (exoni = 0; exoni < nexons; exoni++) {
      Exon_free(&(exons[exoni]));
    }
    FREE(exons);
  }

  return gene_sequence;
}


/* exonnames taken from input file */
static void
add_novel_acceptor (Table_T acceptor_fragment_table, char sign, char *chr, Genomicpos_T chrpos, double prob,
		    IIT_T map_iit, IIT_T chromosome_iit, Genome_T genome, char *exons) {
  char *ntfragment, *exonname, *acc, *copy, *p, *q, *a, *b;
  int acclength;
  char *chrptr;
  Site_T acceptor;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;

  int i;
  long int tally, *tally_matches, *tally_mismatches;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL;

  debug(fprintf(stderr,"Adding novel acceptor\n"));
  chrnum = IIT_find_one(chromosome_iit,chr);
  chroffset = IIT_interval_low(chromosome_iit,chrnum);

  ntfragment = (char *) CALLOC(halflength3+1,sizeof(char));
  ntfragment[halflength3] = '\0';
  if (sign == '+') {
    Genome_fill_buffer_simple(genome,chroffset+chrpos-1,/*length*/halflength3,ntfragment);
    if (bamreader == NULL) {
      tally = 0;
    } else {
      if (bam_lacks_chr == NULL) {
	chrptr = chr;
      } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chr[bam_lacks_chr_length]);
      } else {
	chrptr = chr;
      }

      Bamread_limit_region(bamreader,chrptr,/*chrstart*/chrpos,/*chrend*/chrpos+halflength3);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,/*printchr*/chr,chroffset,/*chrstart*/chrpos-halflength3,/*chrend*/chrpos,
		   /*map_iit*/NULL,alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   /*desired_read_group*/NULL,minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_TALLY,blockp,blocksize,
		   quality_score_adj,/*min_depth*/1,/*variant_strands*/0,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,/*ignore_query_Ns_p*/true,
		   /*print_indels_p*/false,/*print_totals_p*/false,
		   /*print_cycles_p*/false,/*print_quality_scores_p*/false,
		   /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		   /*verbosep*/false,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		   /*bamfile*/NULL);
      tally = 0;
      for (i = 0; i < halflength3; i++) {
	tally += tally_matches[i] + tally_mismatches[i];
      }
      FREE(tally_mismatches);
      FREE(tally_matches);
      Bamread_unlimit_region(bamreader);
    }

    acceptor = Site_new(ntfragment,/*splicesite*/chrpos,/*sign*/+1,prob,tally,/*novelp*/true);
    /* No aafragment */

  } else if (sign == '-') {
    Genome_fill_buffer_simple(genome,chroffset+chrpos-halflength3,/*length*/halflength3,ntfragment);
    make_complement_inplace(ntfragment,halflength3);
    if (bamreader == NULL) {
      tally = 0;
    } else {
      if (bam_lacks_chr == NULL) {
	chrptr = chr;
      } else if (!strncmp(chr,bam_lacks_chr,bam_lacks_chr_length)) {
	chrptr = &(chr[bam_lacks_chr_length]);
      } else {
	chrptr = chr;
      }

      Bamread_limit_region(bamreader,chrptr,/*chrstart*/chrpos-halflength3,/*chrend*/chrpos);
      Bamtally_run(&tally_matches,&tally_mismatches,
		   &intervallist,&labellist,&datalist,
		   quality_counts_match,quality_counts_mismatch,
		   bamreader,genome,/*printchr*/chr,chroffset,/*chrstart*/chrpos-halflength3,/*chrend*/chrpos,
		   /*map_iit*/NULL,alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		   /*desired_read_group*/NULL,minimum_mapq,good_unique_mapq,maximum_nhits,
		   need_concordant_p,need_unique_p,need_primary_p,ignore_duplicates_p,
		   /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		   /*output_type*/OUTPUT_TALLY,blockp,blocksize,
		   quality_score_adj,/*min_depth*/1,/*variant_strands*/0,
		   /*genomic_diff_p*/false,/*signed_counts_p*/false,/*ignore_query_Ns_p*/true,
		   /*print_indels_p*/false,/*print_totals_p*/false,
		   /*print_cycles_p*/false,/*print_quality_scores_p*/false,
		   /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		   /*verbosep*/false,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		   /*bamfile*/NULL);
      tally = 0;
      for (i = 0; i < halflength3; i++) {
	tally += tally_matches[i] + tally_mismatches[i];
      }
      FREE(tally_mismatches);
      FREE(tally_matches);
      Bamread_unlimit_region(bamreader);
    }

    acceptor = Site_new(ntfragment,/*splicesite*/chrpos,/*sign*/-1,prob,tally,/*novelp*/true);
    /* No aafragment */

  } else {
    fprintf(stderr,"Sign %c not recognized\n",sign);
    abort();
  }

  acceptor->origp = true; /* Novel site was found originally by alignment */
  Table_put(acceptor_fragment_table,ntfragment,acceptor);

  if ((p = exons) == NULL) {
    exonname = (char *) CALLOC(strlen("NovelExon") + 1,sizeof(char));
    sprintf(exonname,"NovelExon");
    acceptor->exonnames = List_push(acceptor->exonnames,(void *) exonname);
  } else {
    acc = p;
    while (*p != '\0') {
      if (*p == ',') {
	exonname = (char *) CALLOC(p - acc + 1,sizeof(char));
	strncpy(exonname,acc,p-acc);
	q = exonname;

	a = exonname;
	while (*a == '*') {
	  a++;
	}
	b = strstr(exonname,"_exon");
	acclength = b - a;
	copy = (char *) CALLOC(acclength+1,sizeof(char));
	strncpy(copy,a,acclength);
	Site_store_acc(acceptor,copy);
	FREE(copy);

	while (*q != '\0') {
	  if (*q == '-') {
	    *q = '_';
	  }
	  q++;
	}
	/* Site_store_acc(acceptor,acc) */
	acceptor->exonnames = List_push(acceptor->exonnames,(void *) exonname);
	p++;
	acc = p;
      } else {
	p++;
      }
    }

    exonname = (char *) CALLOC(p - acc + 1,sizeof(char));
    strncpy(exonname,acc,p-acc);
    q = exonname;

    a = exonname;
    while (*a == '*') {
      a++;
    }
    b = strstr(exonname,"_exon");
    acclength = b - a;
    copy = (char *) CALLOC(acclength+1,sizeof(char));
    strncpy(copy,a,acclength);
    Site_store_acc(acceptor,copy);
    FREE(copy);

    while (*q != '\0') {
      if (*q == '-') {
	*q = '_';
      }
      q++;
    }
    acceptor->exonnames = List_push(acceptor->exonnames,(void *) exonname);
  }

  return;
}



static char *
process_acceptors (bool *found_exact_p, Genomicpos_T *minlow, Genomicpos_T *maxhigh,
		   Table_T acceptor_fragment_table, Table_T acceptor_continuation_table,
		   char *acc, char *chr, Genomicpos_T chrpos, double prob,
		   IIT_T map_iit, IIT_T chromosome_iit, Genome_T genome,
		   int max_insertlength, Genomicpos_T overhang,
		   bool exact_junction_p) {
  char *gene_sequence = NULL, *aa_sequence;
  int translation_leftpos, translation_rightpos, translation_length;
  int margin5, margin3;
  int genestrand;
  char *annot, *restofheader, *divstring;
  bool allocp;
  int index;

  Chrnum_T chrnum;
  Genomicpos_T chroffset, low, high;
  Genomicpos_T genelength;

  Exon_T *exons;
  int nexons, firstexoni, lastexoni, exoni;
  char *ntfragment, *aafragment;
  long int tally;

  Site_T acceptor, acceptor_cont;
  Genomicpos_T splicesite;
  char *exonname;


  if (map_iit == NULL) {
    index = -1;
  } else {
    index = IIT_find_one(map_iit,acc);
  }
  if (index < 0) {
    fprintf(stderr,"Cannot find gene for %s\n",acc);
  } else {
    genestrand = IIT_interval_sign(map_iit,index);
    divstring = IIT_divstring_from_index(map_iit,index);
    annot = IIT_annotation(&restofheader,map_iit,index,&allocp);
    chrnum = IIT_find_one(chromosome_iit,chr);
    chroffset = IIT_interval_low(chromosome_iit,chrnum);
    exons = get_exons(&nexons,&margin5,&margin3,annot,/*chr*/divstring,chroffset,genome,halflength3);
    if (allocp) {
      FREE(restofheader);
    }

    gene_sequence = concatenate_exons(&genelength,exons,nexons);
    aa_sequence = Translation_via_genomic(&translation_leftpos,&translation_rightpos,&translation_length,
					  gene_sequence,/*startpos*/margin5,/*endpos*/genelength-margin3,genelength,
					  /*backwardp*/false,/*revcompp*/false,/*fulllengthp*/true,/*cds_startpos*/0);
    add_aa_to_exons(aa_sequence,exons,nexons);
    FREE(aa_sequence);

    /* acceptor firstexoni goes downward to lastexoni */
    if (nexons == 0) {
      /* Skip */
      fprintf(stderr,"Gene %s has no exons\n",acc);

    } else if (genestrand > 0) {
      if (get_acceptor_bounds_plus(&firstexoni,&lastexoni,exons,nexons,chrpos,
				   max_insertlength,overhang,exact_junction_p) == true) {
	if ((low = exons[lastexoni]->exonlow) < *minlow) {
	  *minlow = low;
	}
	if ((high = exons[firstexoni]->exonhigh) > *maxhigh) {
	  *maxhigh = high;
	}
      }

    } else if (genestrand < 0) {
      if (get_acceptor_bounds_minus(&firstexoni,&lastexoni,exons,nexons,chrpos,
				    max_insertlength,overhang,exact_junction_p) == true) {
	if ((high = exons[lastexoni]->exonhigh) > *maxhigh) {
	  *maxhigh = high;
	}
	if ((low = exons[firstexoni]->exonlow) < *minlow) {
	  *minlow = low;
	}
      }

    } else {
      fprintf(stderr,"genestrand is 0\n");
      abort();
    }

    debug(printf("For acceptor acc %s (strand %d), exons %d..%d\n",acc,genestrand,firstexoni,lastexoni));
    /* Was (exoni = firstexoni; exoni >= lastexoni && exoni >= 1; exoni--) */
    for (exoni = nexons - 1; exoni >= 1; exoni--) {
      ntfragment = (char *) CALLOC(halflength3+1,sizeof(char));
      aafragment = (char *) CALLOC(halflength3+1,sizeof(char));
      ntfragment[halflength3] = aafragment[halflength3] = '\0';
      tally = 0;
      if (get_fragment_right(&(ntfragment[0]),&(aafragment[0]),&tally,exoni,nexons,exons,halflength3) != NULL) {
	debug(printf("right fragment %d/%d: %s\n",exoni,nexons,ntfragment));
	if ((acceptor = Table_get(acceptor_fragment_table,ntfragment)) != NULL) {
	  FREE(ntfragment);
	} else {
	  splicesite = exons[exoni]->exonstart;
	  acceptor = Site_new(ntfragment,splicesite,exons[exoni]->sign,prob,tally,/*novelp*/false);
	  Table_put(acceptor_fragment_table,ntfragment,acceptor);
	}

	Site_store_acc(acceptor,acc);
	exonname = (char *) CALLOC(strlen(acc) + strlen("_exon") + 
				   numlen(exoni+1) + strlen("/") + numlen(nexons) + 1,sizeof(char));
	sprintf(exonname,"%s_exon%d/%d",acc,exoni+1,nexons);
	acceptor->exonnames = List_push(acceptor->exonnames,(void *) exonname);
	if (exoni+1 == 2) {
	  acceptor->firstp = true;
	}

	Site_aabinding_add(acceptor,aafragment,exonname);
      }
      if (exons[exoni]->exonstart == chrpos) {
	acceptor->origp = true;
	*found_exact_p = true;
      }
    }


    /* Do continuation */
    /* Was (exoni = firstexoni - 1; exoni >= lastexoni - 1 && exoni >= 0; exoni--) */
    for (exoni = nexons - 2; exoni >= 0; exoni--) {
      ntfragment = (char *) CALLOC(halflength3+1,sizeof(char));
      aafragment = (char *) CALLOC(halflength3+1,sizeof(char));
      ntfragment[halflength3] = aafragment[halflength3] = '\0';
      tally = 0;
      if (get_fragment_left(&(ntfragment[halflength3]),&(aafragment[halflength3]),&tally,exoni,exons,halflength3) != NULL) {
	debug(printf("left continuation fragment %d/%d: %s\n",exoni,nexons,ntfragment));
	if ((acceptor_cont = Table_get(acceptor_continuation_table,ntfragment)) != NULL) {
	  FREE(ntfragment);
	} else {
	  splicesite = exons[exoni]->exonend;
	  acceptor_cont = Site_new(ntfragment,splicesite,exons[exoni]->sign,prob,tally,/*novelp*/false);
	  Table_put(acceptor_continuation_table,ntfragment,acceptor_cont);
	}

	/* Site_store_acc(acceptor_cont,acc); */
	exonname = (char *) CALLOC(strlen(acc) + strlen("_exon") + 
				   numlen(exoni+1) + strlen("/") + numlen(nexons) + 1,sizeof(char));
	sprintf(exonname,"%s_exon%d/%d",acc,exoni+1,nexons);
	acceptor_cont->exonnames = List_push(acceptor_cont->exonnames,(void *) exonname);
	if (exoni+1 == 2) {
	  acceptor_cont->firstp = true;
	}

	Site_aabinding_add(acceptor_cont,aafragment,exonname);
      }
    }

    for (exoni = 0; exoni < nexons; exoni++) {
      Exon_free(&(exons[exoni]));
    }
    FREE(exons);
  }

  return gene_sequence;
}


static bool
junction_okay_p (char *type, bool intragenicp,
		 bool donor_novelp, int donor_sign, char *donor_chr, Genomicpos_T donor_chrpos,
		 bool acceptor_novelp, int acceptor_sign, char *acceptor_chr, Genomicpos_T acceptor_chrpos,
		 Genomicpos_T donor_minlow, Genomicpos_T donor_maxhigh,
		 Genomicpos_T acceptor_minlow, Genomicpos_T acceptor_maxhigh) {

  if (!strcmp(type,"Unknown")) {
    return false;
  }

  if (intragenicp == true && !strcmp(type,"Deletion")) {
    /* Intragenic deletions are alternate splices, not fusion events */
    return false;
  }

#if 0
  if (donor_novelp == true && acceptor_novelp == true && intragenicp == true) {
    return false;
  }
#endif

  if (donor_novelp == true) {
    /* Skip this test */
  } else if (donor_chrpos < donor_minlow) {
    return false;
  } else if (donor_chrpos > donor_maxhigh) {
    return false;
  }

  if (acceptor_novelp == true) {
    /* Skip this test */
  } else if (acceptor_chrpos < acceptor_minlow) {
    return false;
  } else if (acceptor_chrpos > acceptor_maxhigh) {
    return false;
  }

  if (!strcmp(type,"Deletion")) {
    if (donor_sign == +1 && acceptor_sign == +1) {
      if (donor_chrpos < acceptor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else if (donor_sign == -1 && acceptor_sign == -1) {
      if (donor_chrpos > acceptor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else {
      fprintf(stderr,"Deletion at %s:%u..%u involves different signs\n",
	      donor_chr,donor_chrpos,acceptor_chrpos);
      return false;
    }

  } else if (!strcmp(type,"Scramble")) {
    if (donor_sign == +1 && acceptor_sign == +1) {
      if (acceptor_chrpos < donor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else if (donor_sign == -1 && acceptor_sign == -1) {
      if (donor_chrpos < acceptor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else {
      fprintf(stderr,"Scramble at %s:%u..%u involves different signs\n",
	      donor_chr,donor_chrpos,acceptor_chrpos);
      return false;
    }


  } else if (!strcmp(type,"Intragenic")) {
    if (donor_sign == +1 && acceptor_sign == +1) {
      if (acceptor_chrpos > donor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else if (donor_sign == -1 && acceptor_sign == -1) {
      if (donor_chrpos > acceptor_chrpos) {
	return true;
      } else {
	return false;
      }
    } else {
      fprintf(stderr,"Intragenic at %s:%u..%u involves different signs\n",
	      donor_chr,donor_chrpos,acceptor_chrpos);
      return false;
    }

  } else {
    return true;
  }
}



static bool
donor_novel5_p (Site_T donor, Site_T acceptor) {
  List_T q, r;
  char *p, *accstart;
  int exoni, nexons;
  int acclength;
  char dir;
  Genomicpos_T distance;

  if (donor->novelp == false) {
    return false;
  } else if (acceptor->novelp == true) {
    return false;
  } else {
    for (q = donor->exonnames; q != NULL; q = List_next(q)) {
      p = (char *) List_head(q);
      if (*p == '*') {
	p++;
	accstart = p;
	if ((p = strstr(p,"_exon")) == NULL) {
	  fprintf(stderr,"Could not find the word exon after * in %s\n",(char *) List_head(q));
	} else {
	  acclength = p - accstart;
	  p += strlen("_exon");
	  if (sscanf(p,"%d/%d%c%u",&exoni,&nexons,&dir,&distance) != 4) {
	    fprintf(stderr,"Could not find exon info after the word exon in %s\n",(char *) List_head(q));
	  } else {
	    if (exoni == 1 && dir == '_') {
	      for (r = acceptor->exonnames; r != NULL; r = List_next(r)) {
		if (!strncmp((char *) List_head(r),accstart,acclength)) {
		  return true;
		}
	      }
	    }
	  }
	}
      }
    }
    return false;
  }
}

static bool
acceptor_novel3_p (Site_T donor, Site_T acceptor) {
  List_T q, r;
  char *p, *accstart;
  int exoni, nexons;
  int acclength;
  char dir;
  Genomicpos_T distance;

  if (acceptor->novelp == false) {
    return false;
  } else if (donor->novelp == true) {
    return false;
  } else {
    for (q = acceptor->exonnames; q != NULL; q = List_next(q)) {
      p = (char *) List_head(q);
      if (*p == '*') {
	p++;
	accstart = p;
	if ((p = strstr(p,"_exon")) == NULL) {
	  fprintf(stderr,"Could not find the word exon after * in %s\n",(char *) List_head(q));
	} else {
	  acclength = p - accstart;
	  p += strlen("_exon");
	  if (sscanf(p,"%d/%d%c%u",&exoni,&nexons,&dir,&distance) != 4) {
	    fprintf(stderr,"Could not find exon info after the word exon in %s\n",(char *) List_head(q));
	  } else {
	    if (exoni == nexons && dir == '+') {
	      for (r = donor->exonnames; r != NULL; r = List_next(r)) {
		if (!strncmp((char *) List_head(r),accstart,acclength)) {
		  return true;
		}
	      }
	    }
	  }
	}
      }
    }
    return false;
  }
}


static char *
deletion_type (char *chr, Site_T donor, Site_T acceptor, bool intragenicp,
	       IIT_T intertranscript_outer_iit, IIT_T intertranscript_inner_iit) {
  Genomicpos_T low, high;
  int divnoO;
  int sign = donor->sign;
  Genomicpos_T donorpos = donor->chrpos;
  Genomicpos_T acceptorpos = acceptor->chrpos;

  if (intertranscript_outer_iit == NULL || intertranscript_inner_iit == NULL) {
    return "Deletion";
    
  } else {
    if (donorpos < acceptorpos) {
      low = donorpos;
      high = acceptorpos;
    } else {
      low = acceptorpos;
      high = donorpos;
    }

    divnoO = IIT_divint(intertranscript_outer_iit,chr); /* Should be the same as for inner_iit */

    if (donor_novel5_p(donor,acceptor) == true && acceptor_novel3_p(donor,acceptor) == false) {
      return "Altexon5";

    } else if (donor_novel5_p(donor,acceptor) == false && acceptor_novel3_p(donor,acceptor) == true) {
      return "Altexon3";

    } else if (IIT_contains_region_with_divno_signed(intertranscript_outer_iit,divnoO,
					      low,high,sign) == false) {
      /* Splice is not within a readthrough interval, so must be a deletion */
      return "Deletion";

    } else if (IIT_contained_by_region_with_divno_signed(intertranscript_inner_iit,divnoO,
							 low,high,sign) == true) {
      /* Splice spans between two transcripts */
      if (Site_intragenic_pair_p(donor,acceptor) == true) {
	/* The two transcripts must overlap and we picked sites from the same accession */
	return "Unknown";
      } else if (donor->lastp == true && acceptor->firstp == true) {
	return "Readthroughstd";
      } else {
	return "Readthroughalt";
      }
 
    } else if (intragenicp == true) {
      /* In between: Must have been constructed from two novel sites within a transcript */
      return "Unknown";

    } else {
      /* In between: Altexon involved in a readthrough */
      return "AltexonR";
    }
  }
}


static bool
fetch_junctions (char *type, char *chr1, char *chr2,
		 Table_T donor_acc_table, Table_T acceptor_acc_table,
		 Table_T donor_fragment_table, Table_T acceptor_fragment_table,
		 Table_T donor_continuation_table, Table_T acceptor_continuation_table,
		 Genomicpos_T donor_minlow, Genomicpos_T donor_maxhigh,
		 Genomicpos_T acceptor_minlow, Genomicpos_T acceptor_maxhigh,
		 IIT_T map_iit, IIT_T chromosome_iit, Genome_T genome,
		 int max_insertlength, Genomicpos_T overhang,
		 bool exact_junction_p) {
  bool fetchp = false;
  Site_T *donors, *donors_cont, *acceptors, *acceptors_cont, donor, donor_cont, acceptor, acceptor_cont;
  int ndonors, nacceptors, ndonors_cont, nacceptors_cont, i, j;

  int *matches;
  int nmatches, k;
  char *acc, *gene_sequence;
  bool allocp, found_exact_p;


  ndonors = Table_length(donor_fragment_table);
  donors = (Site_T *) Table_values(donor_fragment_table,NULL);

  nacceptors = Table_length(acceptor_fragment_table);
  acceptors = (Site_T *) Table_values(acceptor_fragment_table,NULL);

  ndonors_cont = Table_length(donor_continuation_table);
  donors_cont = (Site_T *) Table_values(donor_continuation_table,NULL);

  nacceptors_cont = Table_length(acceptor_continuation_table);
  acceptors_cont = (Site_T *) Table_values(acceptor_continuation_table,NULL);

  /* Intergenic */
  for (i = 0; i < ndonors; i++) {
    donor = donors[i];
    for (j = 0; j < nacceptors; j++) {
      acceptor = acceptors[j];
      if (junction_okay_p(type,/*intragenicp*/false,donor->novelp,donor->sign,chr1,donor->chrpos,
			  acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			  donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr1,donor->chrpos,donor->chrpos,/*sortp*/false);
	}
	found_exact_p = false;
	for (k = 0; k < nmatches; k++) {
	  if (donor->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(donor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_donors(&found_exact_p,&donor_minlow,&donor_maxhigh,
					     donor_fragment_table,donor_continuation_table,
					     acc,chr1,donor->chrpos,donor->prob,map_iit,chromosome_iit,genome,
					     max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"3.  Putting %s into donor_acc_table\n",acc); */
	      Table_put(donor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);

	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr2,acceptor->chrpos,acceptor->chrpos,/*sortp*/false);
	}
	found_exact_p = false;
	for (k = 0; k < nmatches; k++) {
	  if (acceptor->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(acceptor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_acceptors(&found_exact_p,&acceptor_minlow,&acceptor_maxhigh,
						acceptor_fragment_table,acceptor_continuation_table,
						acc,chr2,acceptor->chrpos,acceptor->prob,map_iit,chromosome_iit,genome,
						max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"4.  Putting %s into acceptor_acc_table\n",acc); */
	      Table_put(acceptor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }
  }

  /* Intragenic donor, alternate splices */
  found_exact_p = false;
  for (i = 0; i < ndonors; i++) {
    donor = donors[i];
    for (j = 0; j < ndonors_cont; j++) {
      donor_cont = donors_cont[j];
      if (junction_okay_p("Intragenic",/*intragenicp*/false,
			  donor->novelp,donor->sign,chr1,donor->chrpos,
			  donor_cont->novelp,donor_cont->sign,chr1,donor_cont->chrpos,
			  donor_minlow,donor_maxhigh,donor_minlow,donor_maxhigh) == true) {
	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr1,donor->chrpos,donor->chrpos,/*sortp*/false);
	}
	for (k = 0; k < nmatches; k++) {
	  if (donor->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(donor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_donors(&found_exact_p,&donor_minlow,&donor_maxhigh,
					     donor_fragment_table,donor_continuation_table,
					     acc,chr1,donor->chrpos,donor->prob,map_iit,chromosome_iit,genome,
					     max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"5.  Putting %s into donor_acc_table\n",acc); */
	      Table_put(donor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);

	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr1,donor_cont->chrpos,donor_cont->chrpos,/*sortp*/false);
	}
	for (k = 0; k < nmatches; k++) {
	  if (donor_cont->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(donor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_donors(&found_exact_p,&donor_minlow,&donor_maxhigh,
					     donor_fragment_table,donor_continuation_table,
					     acc,chr1,donor_cont->chrpos,donor_cont->prob,map_iit,chromosome_iit,genome,
					     max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"6.  Putting %s into donor_acc_table\n",acc); */
	      Table_put(donor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }
  }


  /* Intragenic acceptor, alternate splices */
  found_exact_p = false;
  for (i = 0; i < nacceptors_cont; i++) {
    acceptor_cont = acceptors_cont[i];
    for (j = 0; j < nacceptors; j++) {
      acceptor = acceptors[j];
      if (junction_okay_p("Intragenic",/*intragenicp*/false,
			  acceptor_cont->novelp,acceptor_cont->sign,chr2,acceptor_cont->chrpos,
			  acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			  acceptor_minlow,acceptor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr2,acceptor_cont->chrpos,acceptor_cont->chrpos,/*sortp*/false);
	}
	for (k = 0; k < nmatches; k++) {
	  if (acceptor_cont->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(acceptor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_acceptors(&found_exact_p,&acceptor_minlow,&acceptor_maxhigh,
						acceptor_fragment_table,acceptor_continuation_table,
						acc,chr2,acceptor_cont->chrpos,acceptor_cont->prob,map_iit,chromosome_iit,genome,
						max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"7.  Putting %s into acceptor_acc_table\n",acc); */
	      Table_put(acceptor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);

	if (map_iit == NULL) {
	  matches = NULL;
	  nmatches = 0;
	} else {
	  matches = IIT_get(&nmatches,map_iit,chr2,acceptor->chrpos,acceptor->chrpos,/*sortp*/false);
	}
	for (k = 0; k < nmatches; k++) {
	  if (acceptor->sign == IIT_interval_sign(map_iit,matches[k])) {
	    acc = IIT_label(map_iit,matches[k],&allocp);
	    if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	      /* Skip readthrough acc */
	    } else if (Table_get(acceptor_acc_table,(void *) acc) == NULL) {
	      gene_sequence = process_acceptors(&found_exact_p,&acceptor_minlow,&acceptor_maxhigh,
						acceptor_fragment_table,acceptor_continuation_table,
						acc,chr2,acceptor->chrpos,acceptor->prob,map_iit,chromosome_iit,genome,
						max_insertlength,overhang,exact_junction_p);
	      /* fprintf(stderr,"8.  Putting %s into acceptor_acc_table\n",acc); */
	      Table_put(acceptor_acc_table,(void *) acc,(void *) gene_sequence);
	      fetchp = true;
	    }
	    if (allocp) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }
  }

  FREE(acceptors_cont);
  FREE(donors_cont);
  FREE(acceptors);
  FREE(donors);

  return fetchp;
}


static int
ngenes (char *genename) {
  int n = 1;
  char *p;

  for (p = genename; *p != '\0'; p++) {
    if (*p == ',') {
      n++;
    }
  }
  return n;
}


static void
print_site_name (FILE *fp, char *type, Site_T site, Site_T site1, Site_T site2,
		 char *chr, char *genename, IIT_T map_iit, bool exact_junction_p) {
  List_T p;
  int k;

  fprintf(fp,"%s-",type);
  fprintf(fp,"%c%s@%u..%u",
	  site->sign > 0 ? '+' : '_',chr,site1->chrpos,site2->chrpos);
  fprintf(fp,"-");
  
  fprintf(fp,"%.2f-%.2f-",site1->prob,site2->prob);

  if (genename == NULL || exact_junction_p == false || ngenes(genename) > MAX_GENENAMES) {
    genename = get_genename(chr,site->chrpos,site->sign,map_iit);
    fprintf(fp,"%s",genename);
    FREE(genename);
  } else {
    fprintf(fp,"%s",genename);
  }
  fprintf(fp,"-");

  fprintf(fp,"%s",(char *) List_head(site1->exonnames));
  for (p = List_next(site1->exonnames), k = 1; p != NULL && k < MAX_EXONNAMES; p = List_next(p), k++) {
    fprintf(fp,"|%s",(char *) List_head(p));
  }

  fprintf(fp,"-");
  fprintf(fp,"%s",(char *) List_head(site2->exonnames));
  for (p = List_next(site2->exonnames), k = 1; p != NULL && k < MAX_EXONNAMES; p = List_next(p), k++) {
    fprintf(fp,"|%s",(char *) List_head(p));
  }

  return;
}


static void
print_junction_name (FILE *fp, char *type, char *frametype, bool intragenicp,
		     int similarity, char *sim_donor_acc, int sim_donor_pos, char *sim_acceptor_acc, int sim_acceptor_pos,
		     Genomicpos_T distance, Site_T donor, Site_T acceptor, List_T donor_exonnames, List_T acceptor_exonnames,
		     char *chr1, char *chr2, char *genename1, char *genename2, IIT_T map_iit,
		     bool exact_junction_p) {
  List_T p;
  int k;

  /* AltexonR is not possible with intragenicp == true */
  if (!strcmp(type,"Altexon5") || !strcmp(type,"Altexon3")) {
    /* Already assumed to be intragenic */
    fprintf(fp,"%s-%u-",type,distance);
  } else if (intragenicp == true) {
    fprintf(fp,"%sIntra-%u-",type,distance);
  } else {
    fprintf(fp,"%s-%u-",type,distance);
  }

  if (frametype != NULL) {
    fprintf(fp,"%s-",frametype);
  }

  fprintf(fp,"%c%s@%u..%c%s@%u",
	  donor->sign > 0 ? '+' : '_',chr1,donor->chrpos,
	  acceptor->sign > 0 ? '+' : '_',chr2,acceptor->chrpos);
  fprintf(fp,"-");

  fprintf(fp,"%.2f-%.2f-",donor->prob,acceptor->prob);

  if (genename1 == NULL || exact_junction_p == false || ngenes(genename1) > MAX_GENENAMES) {
    genename1 = get_genename(chr1,donor->chrpos,donor->sign,map_iit);
    if (ngenes(genename1) > 5) {
      abort();
    }
    fprintf(fp,"%s",genename1);
    FREE(genename1);
  } else {
    if (ngenes(genename1) > 5) {
      abort();
    }
    fprintf(fp,"%s",genename1);
  }
  fprintf(fp,"-");
  if (genename2 == NULL || exact_junction_p == false || ngenes(genename2) > MAX_GENENAMES) {
    genename2 = get_genename(chr2,acceptor->chrpos,acceptor->sign,map_iit);
    if (ngenes(genename2) > 5) {
      abort();
    }
    fprintf(fp,"%s",genename2);
    FREE(genename2);
  } else {
    if (ngenes(genename2) > 5) {
      abort();
    }
    fprintf(fp,"%s",genename2);
  }
  fprintf(fp,"-");

  fprintf(fp,"%s",(char *) List_head(donor_exonnames));
  for (p = List_next(donor_exonnames), k = 1; p != NULL && k < MAX_EXONNAMES; p = List_next(p), k++) {
    fprintf(fp,"|%s",(char *) List_head(p));
  }

  fprintf(fp,"-");
  fprintf(fp,"%s",(char *) List_head(acceptor_exonnames));
  for (p = List_next(acceptor_exonnames), k = 1; p != NULL && k < MAX_EXONNAMES; p = List_next(p), k++) {
    fprintf(fp,"|%s",(char *) List_head(p));
  }

  if (intragenicp == true) {
    /* No similarity needed on intragenic scrambles or inversions (similarity is 25/25) */
  } else if (donor->novelp == true || acceptor->novelp == true) {
    /* No similarity possible on novel splicing (don't know corresponding transcript) */
  } else {
    fprintf(fp,"-");
    fprintf(fp,"%d/%d-%s@%d-%s@%d",similarity,SIMILARITY_WIDTH,sim_donor_acc,sim_donor_pos,sim_acceptor_acc,sim_acceptor_pos);
  }

  return;
}


/* Fill donor_acc_used_table and acceptor_acc_used_table */
static void
test_junctions (char *orig_type, char *chr1, char *chr2,
		Table_T donor_acc_used_table, Table_T acceptor_acc_used_table,
		Table_T donor_fragment_table, Table_T acceptor_fragment_table,
		Genomicpos_T donor_minlow, Genomicpos_T donor_maxhigh,
		Genomicpos_T acceptor_minlow, Genomicpos_T acceptor_maxhigh,
		IIT_T intertranscript_outer_iit, IIT_T intertranscript_inner_iit) {
  Site_T *donors, *acceptors, donor, acceptor;
  int ndonors, nacceptors, i, j;
  char *type;

  ndonors = Table_length(donor_fragment_table);
  donors = (Site_T *) Table_values(donor_fragment_table,NULL);

  nacceptors = Table_length(acceptor_fragment_table);
  acceptors = (Site_T *) Table_values(acceptor_fragment_table,NULL);
  /* fprintf(stderr,"%d donors, %d acceptors\n",ndonors,nacceptors); */

  /* Set intragenicp to be false, because we don't know yet */
  for (i = 0; i < ndonors; i++) {
    donor = donors[i];
    for (j = 0; j < nacceptors; j++) {
      acceptor = acceptors[j];

      type = orig_type;
      if (!strcmp(orig_type,"Deletion")) {
	type = deletion_type(chr1,donor,acceptor,/*intragenicp*/false,
			     intertranscript_outer_iit,intertranscript_inner_iit);
      }

      if (junction_okay_p(type,/*intragenicp*/false,donor->novelp,donor->sign,chr1,donor->chrpos,
			  acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			  donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	Site_store_accs(donor,donor_acc_used_table);
	Site_store_accs(acceptor,acceptor_acc_used_table);
      }
    }
  }

  FREE(acceptors);
  FREE(donors);

  return;
}


static void
print_junctions (FILE *counts_fp, char source, int count_spliced, int count_unpaired,
		 char *orig_type, bool intragenicp, unsigned int similarity,
		 char *sim_donor_acc, int sim_donor_pos, char *sim_acceptor_acc, int sim_acceptor_pos,
		 char *chr1, char *chr2, char *genename1, char *genename2,
		 Table_T donor_fragment_table, Table_T acceptor_fragment_table,
		 Table_T donor_continuation_table, Table_T acceptor_continuation_table,
		 Genomicpos_T donor_minlow, Genomicpos_T donor_maxhigh,
		 Genomicpos_T acceptor_minlow, Genomicpos_T acceptor_maxhigh,
		 IIT_T map_iit, IIT_T intertranscript_outer_iit, IIT_T intertranscript_inner_iit,
		 bool exact_junction_p, bool controlsp) {
  Site_T *donors, *donors_cont, *acceptors, *acceptors_cont, donor, donor_cont, acceptor, acceptor_cont;
  int ndonors, nacceptors, ndonors_cont, nacceptors_cont, i, j;
  Genomicpos_T distance;
  char *type;

  ndonors = Table_length(donor_fragment_table);
  donors = (Site_T *) Table_values(donor_fragment_table,NULL);
  qsort(donors,ndonors,sizeof(Site_T),Site_cmp);

  nacceptors = Table_length(acceptor_fragment_table);
  acceptors = (Site_T *) Table_values(acceptor_fragment_table,NULL);
  qsort(acceptors,nacceptors,sizeof(Site_T),Site_cmp);
  /* fprintf(stderr,"%d donors, %d acceptors\n",ndonors,nacceptors); */

  ndonors_cont = Table_length(donor_continuation_table);
  donors_cont = (Site_T *) Table_values(donor_continuation_table,NULL);
  qsort(donors_cont,ndonors_cont,sizeof(Site_T),Site_cmp);

  nacceptors_cont = Table_length(acceptor_continuation_table);
  acceptors_cont = (Site_T *) Table_values(acceptor_continuation_table,NULL);
  qsort(acceptors_cont,nacceptors_cont,sizeof(Site_T),Site_cmp);

  for (i = 0; i < ndonors; i++) {
    donor = donors[i];
    donor->exonnames = list_sort(donor->exonnames);
  }

  for (j = 0; j < nacceptors; j++) {
    acceptor = acceptors[j];
    acceptor->exonnames = list_sort(acceptor->exonnames);
  }

  for (i = 0; i < ndonors_cont; i++) {
    donor_cont = donors_cont[i];
    donor_cont->exonnames = list_sort(donor_cont->exonnames);
  }

  for (j = 0; j < nacceptors_cont; j++) {
    acceptor_cont = acceptors_cont[j];
    acceptor_cont->exonnames = list_sort(acceptor_cont->exonnames);
  }


  if (true /*intergenicp == true*/) {
    /* Intergenic */
    for (i = 0; i < ndonors; i++) {
      donor = donors[i];
      for (j = 0; j < nacceptors; j++) {
	acceptor = acceptors[j];

	type = orig_type;
	if (rename_types_p == false) {
	  /* Don't rename types */
	  
	} else if (!strcmp(orig_type,"Deletion")) {
	  type = deletion_type(chr1,donor,acceptor,intragenicp,intertranscript_outer_iit,
			       intertranscript_inner_iit);
	}

	if (junction_okay_p(type,intragenicp,donor->novelp,donor->sign,chr1,donor->chrpos,
			    acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			    donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	  if (strcmp(chr1,chr2)) {
	    distance = 0U;
	  } else if (donor->chrpos < acceptor->chrpos) {
	    distance = acceptor->chrpos - donor->chrpos;
	  } else {
	    distance = donor->chrpos - acceptor->chrpos;
	  }

	  printf(">");
	  print_junction_name(stdout,type,/*frametype*/NULL,intragenicp,similarity,sim_donor_acc,sim_donor_pos,
			      sim_acceptor_acc,sim_acceptor_pos,
			      distance,donor,acceptor,donor->exonnames,acceptor->exonnames,
			      chr1,chr2,genename1,genename2,map_iit,exact_junction_p);
	  printf("\n");
	  printf("%s\n%s\n",donor->ntfragment,acceptor->ntfragment);

	  if (counts_fp != NULL) {
	    fprintf(counts_fp,"%c",source);
	    if (exact_junction_p == true) {
	      if (donor->origp == true && acceptor->origp == true) {
		/* Don't print anything else */
	      } else {
		fprintf(counts_fp,"A"); /* Alt */
	      }
	    }
	    if (bamreader == NULL) {
	      fprintf(counts_fp,"\t%d\t%d\t",count_spliced,count_unpaired);
	    } else {
	      fprintf(counts_fp,"\t%ld\t%ld\t",donor->tally,acceptor->tally);
	    }
	    print_junction_name(counts_fp,type,/*frametype*/NULL,intragenicp,similarity,sim_donor_acc,sim_donor_pos,
				sim_acceptor_acc,sim_acceptor_pos,
				distance,donor,acceptor,donor->exonnames,acceptor->exonnames,
				chr1,chr2,genename1,genename2,map_iit,exact_junction_p);
	    fprintf(counts_fp,"\n");
	  }
	}
      }
    }
  }


  if (controlsp == true) {
    /* Intragenic donor controls, alternate splices */
    for (i = 0; i < ndonors; i++) {
      donor = donors[i];
      for (j = 0; j < ndonors_cont; j++) {
	donor_cont = donors_cont[j];
	if (junction_okay_p("Intragenic",/*intragenicp*/false,
			    donor->novelp,donor->sign,chr1,donor->chrpos,
			    donor_cont->novelp,donor_cont->sign,chr1,donor_cont->chrpos,
			    donor_minlow,donor_maxhigh,donor_minlow,donor_maxhigh) == true) {
	  printf(">");
	  print_site_name(stdout,"Intragenic",/*site*/donor,/*site1*/donor,/*site2*/donor_cont,
			  chr1,genename1,map_iit,exact_junction_p);
	  printf("\n");

	  printf("%s\n%s\n",donor->ntfragment,donor_cont->ntfragment);
	}
      }
    }


    /* Intragenic acceptor controls, alternate splices */
    for (i = 0; i < nacceptors_cont; i++) {
      acceptor_cont = acceptors_cont[i];
      for (j = 0; j < nacceptors; j++) {
	acceptor = acceptors[j];
	if (junction_okay_p("Intragenic",/*intragenicp*/false,
			    acceptor_cont->novelp,acceptor_cont->sign,chr2,acceptor_cont->chrpos,
			    acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			    acceptor_minlow,acceptor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {

	  printf(">");
	  print_site_name(stdout,"Intragenic",/*site*/acceptor,/*site1*/acceptor_cont,/*site2*/acceptor,
			  chr2,genename2,map_iit,exact_junction_p);
	  printf("\n");

	  printf("%s\n%s\n",acceptor_cont->ntfragment,acceptor->ntfragment);
	}
      }
    }
  }

  FREE(acceptors_cont);
  FREE(donors_cont);
  FREE(acceptors);
  FREE(donors);

  return;
}


static void
print_donor_aafragment (FILE *fp, char *aafragment) {
  char *p, c;

  p = aafragment;
  while ((c = *p++) != '\0') {
    if (c == '5') {
      /* Don't print */
    } else if (c == '3') {
      return;
    } else if (c == '1') {
      fputc('-',fp);
    } else if (c == '2') {
      fputc('-',fp);
    } else {
      fputc(c,fp);
    }
  }
  return;
}

static void
print_acceptor_aafragment (FILE *fp, char *aafragment) {
  char *p, c;

  p = aafragment;
  while ((c = *p++) != '\0') {
    if (c == '5') {
      return;
    } else if (c == '3') {
      return;
    } else if (c == '1') {
      fputc('-',fp);
    } else if (c == '2') {
      fputc('-',fp);
    } else {
      fputc(c,fp);
    }
  }
  return;
}


static void
print_proteins (char *orig_type, bool intragenicp, unsigned int similarity,
		char *sim_donor_acc, int sim_donor_pos, char *sim_acceptor_acc, int sim_acceptor_pos,
		char *chr1, char *chr2, char *genename1, char *genename2,
		Table_T donor_fragment_table, Table_T acceptor_fragment_table,
		Table_T donor_continuation_table, Table_T acceptor_continuation_table,
		Genomicpos_T donor_minlow, Genomicpos_T donor_maxhigh,
		Genomicpos_T acceptor_minlow, Genomicpos_T acceptor_maxhigh,
		IIT_T map_iit, IIT_T intertranscript_outer_iit, IIT_T intertranscript_inner_iit,
		bool exact_junction_p, bool controlsp) {
  Site_T *donors, *donors_cont, *acceptors, *acceptors_cont, donor, donor_cont, acceptor, acceptor_cont;
  List_T p, q;
  AA_binding_T donor_aabinding, acceptor_aabinding;
  char donor_aa, acceptor_aa, *donor_aafragment, *acceptor_aafragment;
  int donor_phase, acceptor_phase;
  char *frametype;
  int ndonors, nacceptors, ndonors_cont, nacceptors_cont, i, j;
  Genomicpos_T distance;
  char *type;

  ndonors = Table_length(donor_fragment_table);
  donors = (Site_T *) Table_values(donor_fragment_table,NULL);
  qsort(donors,ndonors,sizeof(Site_T),Site_cmp);

  nacceptors = Table_length(acceptor_fragment_table);
  acceptors = (Site_T *) Table_values(acceptor_fragment_table,NULL);
  qsort(acceptors,nacceptors,sizeof(Site_T),Site_cmp);
  /* fprintf(stderr,"%d donors, %d acceptors\n",ndonors,nacceptors); */

  ndonors_cont = Table_length(donor_continuation_table);
  donors_cont = (Site_T *) Table_values(donor_continuation_table,NULL);
  qsort(donors_cont,ndonors_cont,sizeof(Site_T),Site_cmp);

  nacceptors_cont = Table_length(acceptor_continuation_table);
  acceptors_cont = (Site_T *) Table_values(acceptor_continuation_table,NULL);
  qsort(acceptors_cont,nacceptors_cont,sizeof(Site_T),Site_cmp);

  for (i = 0; i < ndonors; i++) {
    donor = donors[i];
    donor->exonnames = list_sort(donor->exonnames);
  }

  for (j = 0; j < nacceptors; j++) {
    acceptor = acceptors[j];
    acceptor->exonnames = list_sort(acceptor->exonnames);
  }

  for (i = 0; i < ndonors_cont; i++) {
    donor_cont = donors_cont[i];
    donor_cont->exonnames = list_sort(donor_cont->exonnames);
  }

  for (j = 0; j < nacceptors_cont; j++) {
    acceptor_cont = acceptors_cont[j];
    acceptor_cont->exonnames = list_sort(acceptor_cont->exonnames);
  }


  if (true /*intergenicp == true*/) {
    /* Intergenic */
    for (i = 0; i < ndonors; i++) {
      donor = donors[i];
      for (j = 0; j < nacceptors; j++) {
	acceptor = acceptors[j];

	type = orig_type;
	if (rename_types_p == false) {
	  /* Don't rename types */

	} else if (!strcmp(orig_type,"Deletion")) {
	  type = deletion_type(chr1,donor,acceptor,intragenicp,intertranscript_outer_iit,
			       intertranscript_inner_iit);
	}

	if (junction_okay_p(type,intragenicp,donor->novelp,donor->sign,chr1,donor->chrpos,
			    acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			    donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	  if (strcmp(chr1,chr2)) {
	    distance = 0U;
	  } else if (donor->chrpos < acceptor->chrpos) {
	    distance = acceptor->chrpos - donor->chrpos;
	  } else {
	    distance = donor->chrpos - acceptor->chrpos;
	  }

	  for (p = donor->aabindings; p != NULL; p = List_next(p)) {
	    donor_aabinding = (AA_binding_T) List_head(p);
	    donor_aafragment = donor_aabinding->aafragment;
	    donor_aa = donor_aafragment[strlen(donor_aafragment)-1];
	    switch (donor_aa) {
	    case ' ': donor_phase = -1; break;
	    case '1': donor_phase = 1; break;
	    case '2': donor_phase = 2; break;
	    default: donor_phase = 0;
	    }

	    for (q = acceptor->aabindings; q != NULL; q = List_next(q)) {
	      acceptor_aabinding = (AA_binding_T) List_head(q);
	      acceptor_aafragment = acceptor_aabinding->aafragment;
	      acceptor_aa = acceptor_aafragment[0];
	      switch (acceptor_aa) {
	      case ' ': acceptor_phase = -1; break;
	      case '1': acceptor_phase = 1; break;
	      case '2': acceptor_phase = 2; break;
	      default: acceptor_phase = 3; /* To capture donor_phase of 2 */
	      }

	      if (donor_phase == -1 && acceptor_phase == -1) {
		frametype = "UTRboth";
	      } else if (donor_phase == -1) {
		frametype = "UTR5";
	      } else if (acceptor_phase == -1) {
		frametype = "UTR3";
	      } else if (acceptor_phase == donor_phase + 1) {
		frametype = "Inframe";
	      } else {
		frametype = "Frameshift";
	      }

	      if (inframe_only_p == true && strcmp(frametype,"Inframe")) {
		/* Don't print */
	      } else {
		printf(">");
		print_junction_name(stdout,type,frametype,intragenicp,similarity,sim_donor_acc,sim_donor_pos,
				    sim_acceptor_acc,sim_acceptor_pos,
				    distance,donor,acceptor,donor_aabinding->exonnames,acceptor_aabinding->exonnames,
				    chr1,chr2,genename1,genename2,map_iit,exact_junction_p);
		printf("\n");
		print_donor_aafragment(stdout,donor_aafragment);
		printf("\n");
		print_acceptor_aafragment(stdout,acceptor_aafragment);
		printf("\n");
	      }
	    }
	  }
	}
      }
    }
  }

#if 0
  if (controlsp == true) {
    /* Intragenic donor controls, alternate splices */
    for (i = 0; i < ndonors; i++) {
      donor = donors[i];
      for (j = 0; j < ndonors_cont; j++) {
	donor_cont = donors_cont[j];
	if (junction_okay_p("Intragenic",/*intragenicp*/false,
			    donor->novelp,donor->sign,chr1,donor->chrpos,
			    donor_cont->novelp,donor_cont->sign,chr1,donor_cont->chrpos,
			    donor_minlow,donor_maxhigh,donor_minlow,donor_maxhigh) == true) {
	  printf(">");
	  print_site_name(stdout,"Intragenic",/*site*/donor,/*site1*/donor,/*site2*/donor_cont,
			  chr1,genename1,map_iit,exact_junction_p);
	  printf("\n");

	  printf("%s\n%s\n",donor->ntfragment,donor_cont->ntfragment);
	}
      }
    }


    /* Intragenic acceptor controls, alternate splices */
    for (i = 0; i < nacceptors_cont; i++) {
      acceptor_cont = acceptors_cont[i];
      for (j = 0; j < nacceptors; j++) {
	acceptor = acceptors[j];
	if (junction_okay_p("Intragenic",/*intragenicp*/false,
			    acceptor_cont->novelp,acceptor_cont->sign,chr2,acceptor_cont->chrpos,
			    acceptor->novelp,acceptor->sign,chr2,acceptor->chrpos,
			    acceptor_minlow,acceptor_maxhigh,acceptor_minlow,acceptor_maxhigh) == true) {
	  
	  printf(">");
	  print_site_name(stdout,"Intragenic",/*site*/acceptor,/*site1*/acceptor_cont,/*site2*/acceptor,
			  chr2,genename2,map_iit,exact_junction_p);
	  printf("\n");

	  printf("%s\n%s\n",acceptor_cont->ntfragment,acceptor->ntfragment);
	}
      }
    }
  }
#endif

  FREE(acceptors_cont);
  FREE(donors_cont);
  FREE(acceptors);
  FREE(donors);

  return;
}


static void
free_table_contents (Table_T table) {
  char **fragments;
  Site_T *sites;
  int n, i;

  n = Table_length(table);

  fragments =  (char **) Table_keys(table,NULL);
  for (i = 0; i < n; i++) {
    FREE(fragments[i]);
  }
  FREE(fragments);

  sites = (Site_T *) Table_values(table,NULL);
  for (i = 0; i < n; i++) {
    Site_free(&(sites[i]));
  }
  FREE(sites);

  return;
}



/************************************************************************
 *   Computation of sequence similarity
 ************************************************************************/

#if 0
static void
print_similarity (FILE *fp, char *x, char *y, int n) {
  int nmatches = 0, i;

  for (i = 0; i < n; i++) {
    if (x[i] == y[i]) {
      nmatches++;
    }
  }

}
#endif


static unsigned int
sequence_similarity (char **donor_acc, char **acceptor_acc,
#if 0
		     char **donor_segment, char **acceptor_segment,
#endif
		     int *maxposi, int *maxposj, Table_T donor_acc_table, Table_T acceptor_acc_table,
		     Table_T donor_acc_used_table, Table_T acceptor_acc_used_table,
		     Tableuint_T accession_pair_table, Tableuint_T accession_pair_posi_table, Tableuint_T accession_pair_posj_table) {
  unsigned int max_similarity = 0, similarity;
  int maxi = -1, maxj = -1, posi, posj;
  int n1, n2, i, j;
  char **keys1, **keys2;
  char *donor_sequence, *acceptor_sequence;
  char *accession_pair;

  /* Get accessions from acc_used_table */
  n1 = Table_length(donor_acc_used_table);
  keys1 = (char **) Table_keys(donor_acc_used_table,NULL);
  n2 = Table_length(acceptor_acc_used_table);
  keys2 = (char **) Table_keys(acceptor_acc_used_table,NULL);

  if (n1 == 0 || n2 == 0) {
    /* Must have novel splice sites */
    *donor_acc = (char *) NULL;
    *acceptor_acc = (char *) NULL;
    FREE(keys2);
    FREE(keys1);
    return 0;
  }

  /* Get sequences from acc_table */
  for (i = 0; i < n1; i++) {
    if ((donor_sequence = Table_get(donor_acc_table,keys1[i])) != NULL) {
      for (j = 0; j < n2; j++) {
	if ((acceptor_sequence = Table_get(acceptor_acc_table,keys2[j])) != NULL) {

	  accession_pair = (char *) CALLOC(strlen(keys1[i]) + strlen("_") + strlen(keys2[j]) + 1,sizeof(char));
	  sprintf(accession_pair,"%s_%s",keys1[i],keys2[j]);
	  if (accession_pair_table != NULL && (similarity = Tableuint_get(accession_pair_table,accession_pair)) != 0) {
	    posi = Tableuint_get(accession_pair_posi_table,accession_pair);
	    posj = Tableuint_get(accession_pair_posj_table,accession_pair);
	    FREE(accession_pair);
	  } else {
	    similarity = Oligoindex_compare_seqs(&posi,&posj,donor_sequence,acceptor_sequence,/*width*/SIMILARITY_WIDTH);
	    if (similarity == 0) {
	      /* Storing a similarity value of 0 causes a memory leak, because it looks like non-storage */
	      Tableuint_put(accession_pair_table,accession_pair,/*similarity*/-1U);
	      Tableuint_put(accession_pair_posi_table,accession_pair,posi);
	      Tableuint_put(accession_pair_posj_table,accession_pair,posj);
	    } else {
	      Tableuint_put(accession_pair_table,accession_pair,similarity);
	      Tableuint_put(accession_pair_posi_table,accession_pair,posi);
	      Tableuint_put(accession_pair_posj_table,accession_pair,posj);
	    }
	  }

	  if (similarity == -1U) {
	    similarity = 0;
	  }

	  if (similarity > max_similarity) {
	    maxi = i;
	    maxj = j;
	    *maxposi = posi;
	    *maxposj = posj;
	    max_similarity = similarity;
	  }
	}
      }
    }
  }

  if (maxi < 0 || maxj < 0) {
    /* No 8-mer found to match between sequences */
    maxi = maxj = 0;
    *maxposi = *maxposj = 0;
    max_similarity = 0;
  }

  *donor_acc = (char *) MALLOC((strlen(keys1[maxi])+1) * sizeof(char));
  strcpy(*donor_acc,keys1[maxi]);
  *acceptor_acc = (char *) MALLOC((strlen(keys2[maxj])+1) * sizeof(char));
  strcpy(*acceptor_acc,keys2[maxj]);

#if 0
  *donor_segment = (char *) MALLOC((SIMILARITY_WIDTH+1)*sizeof(char));
  donor_sequence = (char *) Table_get(donor_acc_table,keys1[maxi]);
  strncpy(*donor_segment,&(donor_sequence[*maxposi]),SIMILARITY_WIDTH);
  (*donor_segment)[SIMILARITY_WIDTH] = '\0';

  *acceptor_segment = (char *) MALLOC((SIMILARITY_WIDTH+1)*sizeof(char));
  acceptor_sequence = (char *) Table_get(acceptor_acc_table,keys2[maxj]);
  strncpy(*acceptor_segment,&(acceptor_sequence[*maxposj]),SIMILARITY_WIDTH);
  (*acceptor_segment)[SIMILARITY_WIDTH] = '\0';
#endif

  FREE(keys2);
  FREE(keys1);

  return max_similarity;
}


static bool
intragenic_acc_p (Table_T donor_acc_table, Table_T acceptor_acc_table) {
  int n1, n2, i, j;
  char **keys1, **keys2;

  n1 = Table_length(donor_acc_table);
  keys1 = (char **) Table_keys(donor_acc_table,NULL);
  n2 = Table_length(acceptor_acc_table);
  keys2 = (char **) Table_keys(acceptor_acc_table,NULL);

  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      /* fprintf(stderr,"Comparing %s and %s\n",keys1[i],keys2[j]); */
      if (!strcmp(keys1[i],keys2[j])) {
	/* fprintf(stderr,"Intragenic because %s == %s\n",keys1[i],keys2[j]); */
	FREE(keys1);
	FREE(keys2);
	return true;
      }
    }
  }

  FREE(keys1);
  FREE(keys2);
  return false;
}



/* Store accessions into donor_acc_table and acceptor_acc_table
   because of overlap with the given chrpos.  Store accessions into
   donor_acc_used_table and acceptor_acc_used_table when a junction
   will be printed. */

static void
make_junctions (FILE *counts_fp, char sign1, char *chr1, Genomicpos_T chrpos1, double prob1,
		char sign2, char *chr2, Genomicpos_T chrpos2, double prob2,
		char *type, char *genename1, char *genename2, char *exons1, char *exons2,
		char source, int count_spliced, int count_unpaired,
		Genome_T genome, IIT_T chromosome_iit, IIT_T map_iit,
		int max_insertlength, Genomicpos_T overhang,
		IIT_T intertranscript_outer_iit, IIT_T intertranscript_inner_iit,
		Tableuint_T accession_pair_table,
		Tableuint_T accession_pair_posi_table, Tableuint_T accession_pair_posj_table,
		bool exact_junction_p, bool controlsp) {
  Table_T donor_acc_table, acceptor_acc_table;
  Table_T donor_acc_used_table, acceptor_acc_used_table;
  Table_T donor_fragment_table, acceptor_fragment_table;
  Table_T donor_continuation_table, acceptor_continuation_table;
  Genomicpos_T donor_minlow, donor_maxhigh, acceptor_minlow, acceptor_maxhigh;
  int *matches;
  int nmatches, n, i;
  char *acc, *gene_sequence;
  char **keys;
  bool allocp, fetchp, found_exact_p;

  bool intragenicp;
  unsigned int similarity;
  char *sim_donor_acc = NULL, *sim_acceptor_acc = NULL;
#if 0
  char *donor_segment, *acceptor_segment;
#endif
  int sim_donor_pos, sim_acceptor_pos;

  donor_acc_table = Table_new(100,Table_string_compare,Table_string_hash);
  acceptor_acc_table = Table_new(100,Table_string_compare,Table_string_hash);
  donor_acc_used_table = Table_new(100,Table_string_compare,Table_string_hash);
  acceptor_acc_used_table = Table_new(100,Table_string_compare,Table_string_hash);
  donor_fragment_table = Table_new(100,Table_string_compare,Table_string_hash);
  acceptor_fragment_table = Table_new(100,Table_string_compare,Table_string_hash);
  donor_continuation_table = Table_new(100,Table_string_compare,Table_string_hash);
  acceptor_continuation_table = Table_new(100,Table_string_compare,Table_string_hash);

  donor_minlow = -1U;
  donor_maxhigh = 0U;
  acceptor_minlow = -1U;
  acceptor_maxhigh = 0U;

  if (map_iit == NULL) {
    matches = NULL;
    nmatches = 0;
  } else {
    matches = IIT_get(&nmatches,map_iit,chr1,chrpos1,chrpos1,/*sortp*/false);
  }
  found_exact_p = false;
  for (i = 0; i < nmatches; i++) {
    if ((sign1 == '+' && IIT_interval_sign(map_iit,matches[i]) > 0) ||
	(sign1 == '-' && IIT_interval_sign(map_iit,matches[i]) < 0)) {
      acc = IIT_label(map_iit,matches[i],&allocp);
      if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	/* Skip readthrough acc */
      } else {
	gene_sequence = process_donors(&found_exact_p,&donor_minlow,&donor_maxhigh,
				       donor_fragment_table,donor_continuation_table,
				       acc,chr1,chrpos1,prob1,map_iit,chromosome_iit,genome,
				       max_insertlength,overhang,exact_junction_p);
	/* fprintf(stderr,"1.  Putting acc %s into donor_acc_table\n",acc); */
	Table_put(donor_acc_table,(void *) acc,(void *) gene_sequence);
      }
      if (allocp) {
	FREE(acc);
      }
    }
  }
  FREE(matches);
  debug(fprintf(stderr,"exact_junction_p %d, found_exact_p %d\n",exact_junction_p,found_exact_p));
  if (exact_junction_p == true && found_exact_p == false) {
    add_novel_donor(donor_fragment_table,sign1,chr1,chrpos1,prob1,map_iit,chromosome_iit,genome,exons1);
  }

  if (map_iit == NULL) {
    matches = NULL;
    nmatches = 0;
  } else {
    matches = IIT_get(&nmatches,map_iit,chr2,chrpos2,chrpos2,/*sortp*/false);
  }
  found_exact_p = false;
  for (i = 0; i < nmatches; i++) {
    if ((sign2 == '+' && IIT_interval_sign(map_iit,matches[i]) > 0) ||
	(sign2 == '-' && IIT_interval_sign(map_iit,matches[i]) < 0)) {
      acc = IIT_label(map_iit,matches[i],&allocp);
      if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	/* Skip readthrough acc */
      } else {
	gene_sequence = process_acceptors(&found_exact_p,&acceptor_minlow,&acceptor_maxhigh,
					  acceptor_fragment_table,acceptor_continuation_table,
					  acc,chr2,chrpos2,prob2,map_iit,chromosome_iit,genome,
					  max_insertlength,overhang,exact_junction_p);
	/* fprintf(stderr,"2.  Putting acc %s into acceptor_acc_table\n",acc); */
	Table_put(acceptor_acc_table,(void *) acc,(void *) gene_sequence);
      }
      if (allocp) {
	FREE(acc);
      }
    }
  }
  FREE(matches);
  debug(fprintf(stderr,"exact_junction_p %d, found_exact_p %d\n",exact_junction_p,found_exact_p));
  if (exact_junction_p == true && found_exact_p == false) {
    add_novel_acceptor(acceptor_fragment_table,sign2,chr2,chrpos2,prob2,map_iit,chromosome_iit,genome,exons2);
  }

  fetchp = true;
  while (fetchp == true) {
    fetchp = fetch_junctions(type,chr1,chr2,donor_acc_table,acceptor_acc_table,
			     donor_fragment_table,acceptor_fragment_table,
			     donor_continuation_table,acceptor_continuation_table,
			     donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh,
			     map_iit,chromosome_iit,genome,max_insertlength,overhang,
			     exact_junction_p);
  }

  test_junctions(type,chr1,chr2,donor_acc_used_table,acceptor_acc_used_table,
		 donor_fragment_table,acceptor_fragment_table,
		 donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh,
		 intertranscript_outer_iit,intertranscript_inner_iit);

  if (controlsp == true) {
    /* Print all gene sequences involved.  Take keys from
       acc_used_table, and gene_sequences from acc_table. */
    n = Table_length(donor_acc_used_table);
    keys = (char **) Table_keys(donor_acc_used_table,NULL);
    for (i = 0; i < n; i++) {
      if ((gene_sequence = Table_get(donor_acc_table,keys[i])) != NULL) {
	print_gene(gene_sequence,/*acc*/keys[i]);
      }
    }
    FREE(keys);

    n = Table_length(acceptor_acc_used_table);
    keys = (char **) Table_keys(acceptor_acc_used_table,NULL);
    for (i = 0; i < n; i++) {
      if ((gene_sequence = Table_get(acceptor_acc_table,keys[i])) != NULL) {
	print_gene(gene_sequence,/*acc*/keys[i]);
      }
    }
    FREE(keys);
  }

  if ((intragenicp = intragenic_acc_p(donor_acc_used_table,acceptor_acc_used_table)) == true) {
    /* Skip */
    similarity = SIMILARITY_WIDTH;

  } else {
    similarity = sequence_similarity(&sim_donor_acc,&sim_acceptor_acc,
#if 0
				     &donor_segment,&acceptor_segment,
#endif
				     &sim_donor_pos,&sim_acceptor_pos,donor_acc_table,acceptor_acc_table,
				     donor_acc_used_table,acceptor_acc_used_table,
				     accession_pair_table,accession_pair_posi_table,accession_pair_posj_table);
  }

  if (print_proteins_p == true) {
    print_proteins(type,intragenicp,similarity,sim_donor_acc,sim_donor_pos,sim_acceptor_acc,sim_acceptor_pos,
		   chr1,chr2,genename1,genename2,donor_fragment_table,acceptor_fragment_table,
		   donor_continuation_table,acceptor_continuation_table,
		   donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh,
		   map_iit,intertranscript_outer_iit,intertranscript_inner_iit,
		   exact_junction_p,controlsp);
  } else {
    print_junctions(counts_fp,source,count_spliced,count_unpaired,type,
		    intragenicp,similarity,sim_donor_acc,sim_donor_pos,sim_acceptor_acc,sim_acceptor_pos,
		    chr1,chr2,genename1,genename2,donor_fragment_table,acceptor_fragment_table,
		    donor_continuation_table,acceptor_continuation_table,
		    donor_minlow,donor_maxhigh,acceptor_minlow,acceptor_maxhigh,
		    map_iit,intertranscript_outer_iit,intertranscript_inner_iit,
		    exact_junction_p,controlsp);
  }

  FREE(sim_donor_acc);
  FREE(sim_acceptor_acc);
  
  free_table_contents(donor_fragment_table);
  free_table_contents(acceptor_fragment_table);
  free_table_contents(donor_continuation_table);
  free_table_contents(acceptor_continuation_table);

  n = Table_length(acceptor_acc_used_table);
  keys = (char **) Table_keys(acceptor_acc_used_table,NULL);
  for (i = 0; i < n; i++) {
    FREE(keys[i]);
  }
  FREE(keys);
  Table_free(&acceptor_acc_used_table);

  n = Table_length(donor_acc_used_table);
  keys = (char **) Table_keys(donor_acc_used_table,NULL);
  for (i = 0; i < n; i++) {
    FREE(keys[i]);
  }
  FREE(keys);
  Table_free(&donor_acc_used_table);

  n = Table_length(acceptor_acc_table);
  keys = (char **) Table_keys(acceptor_acc_table,NULL);
  for (i = 0; i < n; i++) {
    gene_sequence = Table_get(acceptor_acc_table,keys[i]);
    FREE(gene_sequence);
    /* FREE(keys[i]); */
  }
  FREE(keys);
  Table_free(&acceptor_acc_table);

  n = Table_length(donor_acc_table);
  keys = (char **) Table_keys(donor_acc_table,NULL);
  for (i = 0; i < n; i++) {
    gene_sequence = Table_get(donor_acc_table,keys[i]);
    FREE(gene_sequence);
    /* FREE(keys[i]); */
  }
  FREE(keys);
  Table_free(&donor_acc_table);

  Table_free(&donor_fragment_table);
  Table_free(&acceptor_fragment_table);
  Table_free(&donor_continuation_table);
  Table_free(&acceptor_continuation_table);

  return;
}


static bool
parse_coords (char *sign, char *chr, Genomicpos_T *chrpos, char *string) {
  char *p;
  int len;

  *sign = string[0];
  if (*sign != '+' && *sign != '_') {
    fprintf(stderr,"Expecting sign to be either + or _\n");
    return false;
  } else if (*sign == '_') {
    *sign = '-';
  }
  
  p = &(string[1]);
  while (*p != '\0' && *p != ':') {
    p++;
  }

  if (*p == '\0') {
    fprintf(stderr,"Could not find colon in string\n");
    return false;
  } else {
    len = p - &(string[1]);
    strncpy(chr,&string[1],len);
    chr[len] = '\0';
    p++;
  }

  if (sscanf(p,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Could not find chromosomal coordinate in string\n");
    return false;
  } else {
    return true;
  }
}


static Table_T
get_ignore_accessions (FILE *fp) {
  Table_T ignore_table;
  char Buffer[1024000], acc[1024], *key;


  ignore_table = Table_new(100,Table_string_compare,Table_string_hash);
  while (fgets(Buffer,1024000,fp) != NULL) {
    if (sscanf(Buffer,"%s",acc) < 1) {
      /* Skip blank line */
    } else {
      key = (char *) CALLOC(strlen(acc)+1,sizeof(char));
      strcpy(key,acc);
      Table_put(ignore_table,key,/*value*/(void *) true);
    }
  }

  return ignore_table;
}



int
main (int argc, char *argv[]) {
  char Buffer[1024000], type[50], genename1[1024], genename2[1024], exons1[1024], exons2[1024], chr1[1024], chr2[1024];
  IIT_T intertranscript_outer_iit = NULL;
  IIT_T intertranscript_inner_iit = NULL;
  int count_spliced, count_unpaired;
  char sign1, sign2;
  double prob1, prob2;
  FILE *fp, *counts_fp = NULL;
  int n, i;
  char **keys;

  Genomicpos_T chrpos1, chrpos2;
  char *genomesubdir, *fileroot, *mapdir, *iitfile;
  IIT_T map_iit, chromosome_iit;
  Genome_T genome;
  int len;

  char source;
  bool exact_junction_p = false;
  int max_insertlength = 200;
  Genomicpos_T overhang = 25;
  Tableuint_T accession_pair_table = NULL, accession_pair_posi_table, accession_pair_posj_table;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:M:m:B:P:c:I:o:i:l:5:3:A?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
#if 0
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
#endif

      } else if (!strcmp(long_name,"dont-rename-types")) {
	rename_types_p = false;

      } else if (!strcmp(long_name,"all-proteins")) {
	print_proteins_p = true;
	inframe_only_p = false;

      } else if (!strcmp(long_name,"inframe-proteins")) {
	print_proteins_p = true;
	inframe_only_p = true;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'exon-exon-junctions --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'M': user_mapdir = optarg; break;
    case 'm': 
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
      break;

    case 'B': bamfile = optarg; break;
    case 'P':
      bam_lacks_chr = optarg;
      bam_lacks_chr_length = strlen(bam_lacks_chr);
      break;

    case 'c': countsfile_name = optarg; break;
    case 'I': ignore_filename = optarg; break;
    case 'o': intertranscript_outer_iitfile = optarg; break;
    case 'i': intertranscript_inner_iitfile = optarg; break;

    case 'l': max_insertlength = atoi(optarg); break;
    case '5': halflength5 = atoi(optarg); break;
    case '3': halflength3 = atoi(optarg); break;
    case 'A': all_exons_p = true; break;

    case '?': /* print_program_usage(); */ exit(0);
    default:
      fprintf(stderr,"Cannot handle flag %c\n",opt);
      exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);


  /* Access gene_map_iit */
  if (dbroot == NULL) {
    fprintf(stderr,"Must specify genome database with the -d flag\n");
    exit(9);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  }

  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);

  if (map_iitfile == NULL) {
    /* fprintf(stderr,"Must specify genes iit with the -m flag\n"); */
    map_iit = (IIT_T) NULL;
    accession_pair_table = (Tableuint_T) NULL;
    accession_pair_posi_table = (Tableuint_T) NULL;
    accession_pair_posj_table = (Tableuint_T) NULL;
  } else {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
    if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",map_iitfile);
      fprintf(stderr,"using the -M flag to gmap.\n");
      exit(9);
    }
    FREE(mapdir);
    FREE(map_iitfile);
    FREE(iitfile);

    accession_pair_table = Tableuint_new(100000,Table_string_compare,Table_string_hash);
    accession_pair_posi_table = Tableuint_new(100000,Table_string_compare,Table_string_hash);
    accession_pair_posj_table = Tableuint_new(100000,Table_string_compare,Table_string_hash);
  }

  FREE(dbversion);
  FREE(fileroot);
  FREE(dbroot);
  FREE(genomesubdir);

  if (bamfile != NULL) {
    bamreader = Bamread_new(bamfile);
  }

  if (ignore_filename != NULL) {
    fp = fopen(ignore_filename,"r");
    if (fp == NULL) {
      fprintf(stderr,"Cannot open file %s\n",ignore_filename);
      ignore_table = (Table_T) NULL;
    } else {
      ignore_table = get_ignore_accessions(fp);
      fclose(fp);
    }
  }

  if (intertranscript_outer_iitfile != NULL) {
    intertranscript_outer_iit = IIT_read(intertranscript_outer_iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
					 /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true);
  }
  if (intertranscript_inner_iitfile != NULL) {
    intertranscript_inner_iit = IIT_read(intertranscript_inner_iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
					 /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true);
  }

  if (argc > 2) {
    /* Usage: exon-exon-junctions _chr21:40684767 +chr14:66188493 */

    strcpy(type,"NA");
    if (parse_coords(&sign1,chr1,&chrpos1,argv[1]) == false) {
      fprintf(stderr,"Could not parse coordinates from %s\n",argv[1]);
      exit(9);
    }
    if (parse_coords(&sign2,chr2,&chrpos2,argv[2]) == false) {
      fprintf(stderr,"Could not parse coordinates from %s\n",argv[2]);
      exit(9);
    }

    if (all_exons_p == true) {
      max_insertlength = 1000000;
      exact_junction_p = false;
      source = 'A';
    } else {
      max_insertlength = 0;	/* Do not extend to other splice sites */
      exact_junction_p = true;	/* Want exact junction */
      extra_exon = 0;		/* Do not want alternate splices */
      source = 'X';
      overhang = 0U;
    }

    if (argc > 4) {
      /* Usage: exon-exon-junctions _chr21:40684767 +chr14:66188493 BRWD1 FUT8 */
      make_junctions(/*counts_fp*/NULL,sign1,chr1,chrpos1,/*prob1*/0.0,sign2,chr2,chrpos2,/*prob2*/0.0,
		     /*type*/"Fusion",/*genename1*/argv[3],/*genename2*/argv[4],/*exons1*/NULL,/*exons2*/NULL,
		     source,/*count_spliced*/0,/*count_unpaired*/0,
		     genome,chromosome_iit,map_iit,
		     max_insertlength,overhang,
		     intertranscript_outer_iit,intertranscript_inner_iit,accession_pair_table,
		     accession_pair_posi_table,accession_pair_posj_table,
		     exact_junction_p,/*controlsp*/false);
    } else {
      /* Usage: exon-exon-junctions _chr21:40684767 +chr14:66188493 */
      make_junctions(/*counts_fp*/NULL,sign1,chr1,chrpos1,/*prob1*/0.0,sign2,chr2,chrpos2,/*prob2*/0.0,
		     /*type*/"Fusion",/*genename1*/"DonorGene",/*genename2*/"AcceptorGene",
		     /*exons1*/NULL,/*exons2*/NULL,
		     source,/*count_spliced*/0,/*count_unpaired*/0,
		     genome,chromosome_iit,map_iit,
		     max_insertlength,overhang,
		     intertranscript_outer_iit,intertranscript_inner_iit,accession_pair_table,
		     accession_pair_posi_table,accession_pair_posj_table,
		     exact_junction_p,/*controlsp*/false);
    }


  } else {
    /* Usage: cat <file> | exon-exon-junctions */
    if (countsfile_name != NULL) {
      if ((counts_fp = fopen(countsfile_name,"w")) == NULL) {
	fprintf(stderr,"Cannot write to file %s\n",countsfile_name);
      }
    }

    while (fgets(Buffer,1024000,stdin) != NULL) {
      if (Buffer[0] == '#') {
	/* Skip */
      } else if (sscanf(Buffer,"%c %d %d %s %c %s %u %lf %c %s %u %lf %s %s %s %s",
			&source,&count_spliced,&count_unpaired,type,
			&sign1,chr1,&chrpos1,&prob1,&sign2,chr2,&chrpos2,&prob2,
			genename1,genename2,exons1,exons2) < 12) {
	fprintf(stderr,"Cannot parse type/chr1/chrpos1/chr2/chrpos2/genename1/genename2/exons1/exons2 in line %s\n",
		Buffer);
      } else if (source == 'X' || source == 'H') {
	/* Exact splice sites (or crossing) */
	debug(printf("Making junctions for exact splice sites\n"));
	make_junctions(counts_fp,sign1,chr1,chrpos1,prob1,sign2,chr2,chrpos2,prob2,
		       type,genename1,genename2,exons1,exons2,source,count_spliced,count_unpaired,
		       genome,chromosome_iit,map_iit,
		       /*max_insertlength*/0,/*overhang*/0U,
		       intertranscript_outer_iit,intertranscript_inner_iit,accession_pair_table,
		       accession_pair_posi_table,accession_pair_posj_table,
		       /*exact_junction_p*/true,/*controlsp*/true);
      } else if (source == 'P') {
	/* Approximate splice sites from paired-end reads (or spanning) */
	debug(printf("Making junctions for approximate splice sites from paired-end reads\n"));
	make_junctions(counts_fp,sign1,chr1,chrpos1,prob1,sign2,chr2,chrpos2,prob2,
		       type,genename1,genename2,exons1,exons2,source,count_spliced,count_unpaired,
		       genome,chromosome_iit,map_iit,
		       max_insertlength,overhang,
		       intertranscript_outer_iit,intertranscript_inner_iit,accession_pair_table,
		       accession_pair_posi_table,accession_pair_posj_table,
		       /*exact_junction_p*/false,/*controlsp*/true);
      } else {
	fprintf(stderr,"Do not recognize source %c.  Needs to be X (for exact) or P (for paired-end).\n",source);
	exit(9);
      }
    }


    if (counts_fp != NULL) {
      fclose(counts_fp);
    }
  }

  if (accession_pair_table != NULL) {
    Tableuint_free(&accession_pair_posj_table);
    Tableuint_free(&accession_pair_posi_table);
    n = Tableuint_length(accession_pair_table);
    keys = (char **) Tableuint_keys(accession_pair_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
    Tableuint_free(&accession_pair_table);
  }

  if (intertranscript_outer_iit != NULL) {
    IIT_free(&intertranscript_outer_iit);
  }
  if (intertranscript_inner_iit != NULL) {
    IIT_free(&intertranscript_inner_iit);
  }

  if (ignore_table != NULL) {
    n = Table_length(ignore_table);
    keys = (char **) Table_keys(ignore_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(keys[i]);
    }
    FREE(keys);
    Table_free(&ignore_table);
  }

  if (bamreader != NULL) {
    Bamread_free(&bamreader);
  }

  IIT_free(&chromosome_iit);
  IIT_free(&map_iit);
  Genome_free(&genome);

  fprintf(stderr,"exon-exon-junctions finished successfully\n");

  return 0;
}

