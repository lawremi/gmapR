static char rcsid[] = "$Id: genecompare.c 51941 2011-11-08 01:26:01Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For isdigit */
#include <unistd.h>		/* For getopt */
#include <math.h>		/* For sqrt */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include "mem.h"
#include "fopen.h"

#include "intlist.h"
#include "uintlist.h"
#include "uinttable.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "datadir.h"
#include "parserange.h"
#include "interval.h"
#include "iit-read.h"


static int global_main_nperfect_fwd = 0;
static int global_main_nerror_fwd = 0;
static int global_main_nnomatch_fwd = 0;
static int global_main_nperfect_rev = 0;
static int global_main_nerror_rev = 0;
static int global_main_nnomatch_rev = 0;

static int global_alt_nperfect_fwd = 0;
static int global_alt_nerror_fwd = 0;
static int global_alt_nnomatch_fwd = 0;
static int global_alt_nperfect_rev = 0;
static int global_alt_nerror_rev = 0;
static int global_alt_nnomatch_rev = 0;


static bool originalp = false;

typedef enum {CORRECT, FALSE, INCOMPLETE} Resulttype_T;


static void
compare_lists_fwd (int *ncommon, int *nlist1_only, int *nlist2_only, 
		   Uintlist_T list1, Uintlist_T list2) {

  while (list1 != NULL && list2 != NULL) {
    if (Uintlist_head(list1) < Uintlist_head(list2)) {
      (*nlist1_only) += 1;
      list1 = Uintlist_next(list1);
    } else if (Uintlist_head(list2) < Uintlist_head(list1)) {
      (*nlist2_only) += 1;
      list2 = Uintlist_next(list2);
    } else {
      (*ncommon) += 1;
      list1 = Uintlist_next(list1);
      list2 = Uintlist_next(list2);
    }
  }

  while (list1 != NULL) {
    (*nlist1_only) += 1;
    list1 = Uintlist_next(list1);
  }

  while (list2 != NULL) {
    (*nlist2_only) += 1;
    list2 = Uintlist_next(list2);
  }

  return;
}


static void
compare_lists_rev (int *ncommon, int *nlist1_only, int *nlist2_only, 
		   Uintlist_T list1, Uintlist_T list2) {

  while (list1 != NULL && list2 != NULL) {
    if (Uintlist_head(list1) > Uintlist_head(list2)) {
      (*nlist1_only) += 1;
      list1 = Uintlist_next(list1);
    } else if (Uintlist_head(list2) > Uintlist_head(list1)) {
      (*nlist2_only) += 1;
      list2 = Uintlist_next(list2);
    } else {
      (*ncommon) += 1;
      list1 = Uintlist_next(list1);
      list2 = Uintlist_next(list2);
    }
  }

  while (list1 != NULL) {
    (*nlist1_only) += 1;
    list1 = Uintlist_next(list1);
  }

  while (list2 != NULL) {
    (*nlist2_only) += 1;
    list2 = Uintlist_next(list2);
  }

  return;
}




static void
compare_splicesites (int *ncommon, int *npredicted_only, int *ngoldstandard_only,
		     Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		     Uintlist_T goldstandard_donors, Uintlist_T goldstandard_acceptors,
		     int sign) {

  *ncommon = *npredicted_only = *ngoldstandard_only = 0;
  if (sign > 0) {
    compare_lists_fwd(&(*ncommon),&(*npredicted_only),&(*ngoldstandard_only),
		      predicted_donors,goldstandard_donors);
    compare_lists_fwd(&(*ncommon),&(*npredicted_only),&(*ngoldstandard_only),
		      predicted_acceptors,goldstandard_acceptors);
  } else {
    compare_lists_rev(&(*ncommon),&(*npredicted_only),&(*ngoldstandard_only),
		      predicted_donors,goldstandard_donors);
    compare_lists_rev(&(*ncommon),&(*npredicted_only),&(*ngoldstandard_only),
		      predicted_acceptors,goldstandard_acceptors);
  }

  return;
}



static void
print_comparison_fwd (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		      Uintlist_T goldstandard_donors, Uintlist_T goldstandard_acceptors) {
  Uintlist_T p, q, u, v;
  unsigned int pvalue, qvalue, uvalue, vvalue, smallest;
  bool at_smallest[4];

  p = predicted_donors;
  q = goldstandard_donors;

  u = predicted_acceptors;
  v = goldstandard_acceptors;

  while (p != NULL || q != NULL || u != NULL || v != NULL) {
    pvalue = (p == NULL) ? (unsigned int) -1 : Uintlist_head(p);
    qvalue = (q == NULL) ? (unsigned int) -1 : Uintlist_head(q);
    uvalue = (u == NULL) ? (unsigned int) -1 : Uintlist_head(u);
    vvalue = (v == NULL) ? (unsigned int) -1 : Uintlist_head(v);

    smallest = pvalue;
    smallest = qvalue < smallest ? qvalue : smallest;
    smallest = uvalue < smallest ? uvalue : smallest;
    smallest = vvalue < smallest ? vvalue : smallest;

    at_smallest[0] = at_smallest[1] = at_smallest[2] = at_smallest[3] = false;
    if (pvalue == smallest) { at_smallest[0] = true; }
    if (qvalue == smallest) { at_smallest[1] = true; }
    if (uvalue == smallest) { at_smallest[2] = true; }
    if (vvalue == smallest) { at_smallest[3] = true; }

    if (at_smallest[0] == true || at_smallest[1] == true) {
      /* Deal with donor */
      if (pvalue < qvalue) {
	printf("%u donor false\n",pvalue);
	p = Uintlist_next(p);
      } else if (qvalue < pvalue) {
	printf("%u donor missed\n",qvalue);
	q = Uintlist_next(q);
      } else {
	printf("%u donor\n",pvalue);
	p = Uintlist_next(p);
	q = Uintlist_next(q);
      }
    } else {
      /* Deal with acceptor */
      if (uvalue < vvalue) {
	printf("%u acceptor false\n",uvalue);
	u = Uintlist_next(u);
      } else if (vvalue < uvalue) {
	printf("%u acceptor missed\n",vvalue);
	v = Uintlist_next(v);
      } else {
	printf("%u acceptor\n",uvalue);
	u = Uintlist_next(u);
	v = Uintlist_next(v);
      }
    }
  }

  return;
}


static void
print_comparison_rev (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		      Uintlist_T goldstandard_donors, Uintlist_T goldstandard_acceptors) {
  Uintlist_T p, q, u, v;
  unsigned int pvalue, qvalue, uvalue, vvalue, largest;
  bool at_largest[4];

  p = predicted_donors;
  q = goldstandard_donors;

  u = predicted_acceptors;
  v = goldstandard_acceptors;

  while (p != NULL || q != NULL || u != NULL || v != NULL) {
    pvalue = (p == NULL) ? (unsigned int) 0 : Uintlist_head(p);
    qvalue = (q == NULL) ? (unsigned int) 0 : Uintlist_head(q);
    uvalue = (u == NULL) ? (unsigned int) 0 : Uintlist_head(u);
    vvalue = (v == NULL) ? (unsigned int) 0 : Uintlist_head(v);

    largest = pvalue;
    largest = qvalue > largest ? qvalue : largest;
    largest = uvalue > largest ? uvalue : largest;
    largest = vvalue > largest ? vvalue : largest;

    at_largest[0] = at_largest[1] = at_largest[2] = at_largest[3] = false;
    if (pvalue == largest) { at_largest[0] = true; }
    if (qvalue == largest) { at_largest[1] = true; }
    if (uvalue == largest) { at_largest[2] = true; }
    if (vvalue == largest) { at_largest[3] = true; }

    if (at_largest[0] == true || at_largest[1] == true) {
      /* Deal with donor */
      if (pvalue > qvalue) {
	printf("%u donor false\n",pvalue);
	p = Uintlist_next(p);
      } else if (qvalue > pvalue) {
	printf("%u donor missed\n",qvalue);
	q = Uintlist_next(q);
      } else {
	printf("%u donor\n",pvalue);
	p = Uintlist_next(p);
	q = Uintlist_next(q);
      }
    } else {
      /* Deal with acceptor */
      if (uvalue > vvalue) {
	printf("%u acceptor false\n",uvalue);
	u = Uintlist_next(u);
      } else if (vvalue > uvalue) {
	printf("%u acceptor missed\n",vvalue);
	v = Uintlist_next(v);
      } else {
	printf("%u acceptor\n",uvalue);
	u = Uintlist_next(u);
	v = Uintlist_next(v);
      }
    }
  }

  return;
}


static void
print_onegene_fwd (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		   List_T predicted_info, char *tag) {
  Uintlist_T p, u;
  List_T a;
  unsigned int pvalue, uvalue, smallest;
  bool at_smallest[2];

  p = predicted_donors;
  u = predicted_acceptors;
  a = predicted_info;

  while (p != NULL || u != NULL) {
    pvalue = (p == NULL) ? (unsigned int) -1 : Uintlist_head(p);
    uvalue = (u == NULL) ? (unsigned int) -1 : Uintlist_head(u);

    smallest = pvalue;
    smallest = uvalue < smallest ? uvalue : smallest;

    at_smallest[0] = at_smallest[1] = false;
    if (pvalue == smallest) { at_smallest[0] = true; }
    if (uvalue == smallest) { at_smallest[1] = true; }

    if (at_smallest[0] == true) {
      /* Deal with donor */
      printf("%u donor",pvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      p = Uintlist_next(p);

    } else {
      /* Deal with acceptor */
      printf("%u acceptor",uvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      u = Uintlist_next(u);

      if (a != NULL) {
	a = List_next(a);
      }
    }
  }

  return;
}


static void
print_onegene_wresults_fwd (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
			    List_T predicted_info, Intlist_T results_donors, Intlist_T results_acceptors,
			    char *tag) {
  Uintlist_T p, u;
  List_T a;
  unsigned int pvalue, uvalue, smallest;
  bool at_smallest[2];

  p = predicted_donors;
  u = predicted_acceptors;
  a = predicted_info;

  while (p != NULL || u != NULL) {
    pvalue = (p == NULL) ? (unsigned int) -1 : Uintlist_head(p);
    uvalue = (u == NULL) ? (unsigned int) -1 : Uintlist_head(u);

    smallest = pvalue;
    smallest = uvalue < smallest ? uvalue : smallest;

    at_smallest[0] = at_smallest[1] = false;
    if (pvalue == smallest) { at_smallest[0] = true; }
    if (uvalue == smallest) { at_smallest[1] = true; }


    if (at_smallest[0] == true) {
      /* Deal with donor */
      printf("%u donor",pvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (Intlist_head(results_donors) == FALSE) {
	printf(" false");
      } else if (Intlist_head(results_donors) == INCOMPLETE) {
	printf(" incomplete");
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      p = Uintlist_next(p);
      results_donors = Intlist_next(results_donors);

    } else {
      /* Deal with acceptor */
      printf("%u acceptor",uvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (Intlist_head(results_acceptors) == FALSE) {
	printf(" false");
      } else if (Intlist_head(results_acceptors) == INCOMPLETE) {
	printf(" incomplete");
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      u = Uintlist_next(u);
      results_acceptors = Intlist_next(results_acceptors);

      if (a != NULL) {
	a = List_next(a);
      }
    }

  }


  return;
}


static void
print_onegene_rev (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		   List_T predicted_info, char *tag) {
  Uintlist_T p, u;
  List_T a;
  unsigned int pvalue, uvalue, largest;
  bool at_largest[2];

  p = predicted_donors;
  u = predicted_acceptors;
  a = predicted_info;

  while (p != NULL || u != NULL) {
    pvalue = (p == NULL) ? (unsigned int) 0 : Uintlist_head(p);
    uvalue = (u == NULL) ? (unsigned int) 0 : Uintlist_head(u);
 
   largest = pvalue;
    largest = uvalue > largest ? uvalue : largest;

    at_largest[0] = at_largest[1] = false;
    if (pvalue == largest) { at_largest[0] = true; }
    if (uvalue == largest) { at_largest[1] = true; }

    if (at_largest[0] == true) {
      /* Deal with donor */
      printf("%u donor",pvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      p = Uintlist_next(p);

    } else {
      /* Deal with acceptor */
      printf("%u acceptor",uvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      u = Uintlist_next(u);

      if (a != NULL) {
	a = List_next(a);
      }
    }

  }

  return;
}

static void
print_onegene_wresults_rev (Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
			    List_T predicted_info, Intlist_T results_donors, Intlist_T results_acceptors,
			    char *tag) {
  Uintlist_T p, u;
  List_T a;
  unsigned int pvalue, uvalue, largest;
  bool at_largest[2];

  p = predicted_donors;
  u = predicted_acceptors;
  a = predicted_info;

  while (p != NULL || u != NULL) {
    pvalue = (p == NULL) ? (unsigned int) 0 : Uintlist_head(p);
    uvalue = (u == NULL) ? (unsigned int) 0 : Uintlist_head(u);

    largest = pvalue;
    largest = uvalue > largest ? uvalue : largest;

    at_largest[0] = at_largest[1] = false;
    if (pvalue == largest) { at_largest[0] = true; }
    if (uvalue == largest) { at_largest[1] = true; }

    if (at_largest[0] == true) {
      /* Deal with donor */
      printf("%u donor",pvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (Intlist_head(results_donors) == FALSE) {
	printf(" false");
      } else if (Intlist_head(results_donors) == INCOMPLETE) {
	printf(" incomplete");
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      p = Uintlist_next(p);
      results_donors = Intlist_next(results_donors);

    } else {
      /* Deal with acceptor */
      printf("%u acceptor",uvalue);
      if (a != NULL) {
	printf(" %s",(char *) List_head(a));
      }
      if (Intlist_head(results_acceptors) == FALSE) {
	printf(" false");
      } else if (Intlist_head(results_acceptors) == INCOMPLETE) {
	printf(" incomplete");
      }
      if (tag) {
	printf(" %s",tag);
      }
      printf("\n");
      u = Uintlist_next(u);
      results_acceptors = Intlist_next(results_acceptors);

      if (a != NULL) {
	a = List_next(a);
      }
    }

  }

  return;
}



/* Returns true if main */
static bool
parse_gene (Uintlist_T *donor_list, Uintlist_T *acceptor_list, List_T *info_list, char *gene) {
  bool mainp;
  int exoni;
  Genomicpos_T exonstart, exonend;
  char *p, *q;
  char *info;


  *donor_list = (Uintlist_T) NULL;
  *acceptor_list = (Uintlist_T) NULL;
  *info_list = (List_T) NULL;

  p = gene;
  if (!strncmp(p,"Alt",3)) {
    mainp = false;
  } else {
    mainp = true;
  }
   

  /* Skip gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }

  exoni = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) == 2) {
      if (exoni > 0) {
	*acceptor_list = Uintlist_push(*acceptor_list,exonstart);
      }
      *donor_list = Uintlist_push(*donor_list,exonend);
      exoni++;
    }

    /* Skip part just read */
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && *p != '\n' && isspace(*p)) { p++; }

    /* Read to end of line */
    if (*p == '\0' || *p == '\n') {
      info = (char *) CALLOC(1,sizeof(char));
      info[0] = '\0';
    } else {
      q = p;
      while (*p != '\0' && *p != '\n') { p++; }

      info = (char *) CALLOC(p - q + 1,sizeof(char));
      strncpy(info,q,p-q);
      if (*p == '\n') { p++; }
    }

    *info_list = List_push(*info_list,(void *) info);
  }

  *donor_list = Uintlist_pop(*donor_list,&exonend);

  *donor_list = Uintlist_reverse(*donor_list);
  *acceptor_list = Uintlist_reverse(*acceptor_list);
  *info_list = List_reverse(*info_list);

  return mainp;
}


static int
evaluate_genes (int *best_ncommon, int *best_npredicted_only, int *best_ngoldstandard_only,
		Uintlist_T *best_goldstandard_donors, Uintlist_T *best_goldstandard_acceptors,
		Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		int *matches, int nmatches, IIT_T goldstandard_iit, int sign, bool mainp) {
  int i;
  int bestindex = -1;
  int bestscore = -1000, score, ncommon, npredicted_only, ngoldstandard_only;
  char *goldstandard_gene, *restofheader;
  Uintlist_T goldstandard_donors, goldstandard_acceptors;
  List_T goldstandard_info;
  bool allocp;

  for (i = 0; i < nmatches; i++) {
    goldstandard_gene = IIT_annotation(&restofheader,goldstandard_iit,matches[i],&allocp);
    parse_gene(&goldstandard_donors,&goldstandard_acceptors,&goldstandard_info,goldstandard_gene);
    compare_splicesites(&ncommon,&npredicted_only,&ngoldstandard_only,
			predicted_donors,predicted_acceptors,
			goldstandard_donors,goldstandard_acceptors,sign);
    if (i == 0) {
      if (mainp == true) {
	bestscore = ncommon - npredicted_only - ngoldstandard_only;
      } else {
	bestscore = ncommon - npredicted_only;
      }
      bestindex = matches[i];
      *best_ncommon = ncommon;
      *best_npredicted_only = npredicted_only;
      *best_ngoldstandard_only = ngoldstandard_only;
      *best_goldstandard_donors = goldstandard_donors;
      *best_goldstandard_acceptors = goldstandard_acceptors;

    } else {
      if (mainp == true) {
	score = ncommon - npredicted_only - ngoldstandard_only;
      } else {
	score = ncommon - npredicted_only;
      }
      if (score > bestscore) {
	bestscore = score;
	bestindex = matches[i];
	*best_ncommon = ncommon;
	*best_npredicted_only = npredicted_only;
	*best_ngoldstandard_only = ngoldstandard_only;
	Uintlist_free(&(*best_goldstandard_donors));
	Uintlist_free(&(*best_goldstandard_acceptors));
	*best_goldstandard_donors = goldstandard_donors;
	*best_goldstandard_acceptors = goldstandard_acceptors;
      }
    }

    if (allocp) {
      FREE(restofheader);
    }
  }

  return bestindex;
}


static bool
dest_presentp (Uintlist_T list, Genomicpos_T pos) {
  Uintlist_T p;

  for (p = list; p != NULL; p = Uintlist_next(p)) {
    if (Uintlist_head(p) == pos) {
      return true;
    }
  }
  return false;
}


static void
enter_introns (Uinttable_T table, Uintlist_T src_sites, Uintlist_T dest_sites) {
  Uintlist_T list, p, q;
  Genomicpos_T src, dest;

  p = src_sites;
  q = dest_sites;
  while (p != NULL && q != NULL) {
    src = Uintlist_head(p);
    dest = Uintlist_head(q);
    list = Uinttable_get(table,src);
    if (dest_presentp(list,dest) == false) {
      Uinttable_put(table,src,(void *) Uintlist_push(list,dest));
    }
    p = Uintlist_next(p);
    q = Uintlist_next(q);
  }
  return;
}

static void
evaluate_splices (int *ncommon, int *npredicted_only, int *nincomplete,
		  Intlist_T *results_donors, Intlist_T *results_acceptors,
		  Uintlist_T predicted_donors, Uintlist_T predicted_acceptors,
		  int *matches, int nmatches, IIT_T goldstandard_iit, bool mainp) {
  int n, i;
  char *goldstandard_gene, *restofheader;
  Uintlist_T goldstandard_donors, goldstandard_acceptors, list, p, q;
  List_T goldstandard_info;
  Uinttable_T donor_table, acceptor_table, genestart_table, geneend_table;
  Genomicpos_T donor, acceptor, *keys;
  bool allocp;

  /* Assume 10 introns per gene */
  donor_table = Uinttable_new(nmatches*10);
  acceptor_table = Uinttable_new(nmatches*10);
  genestart_table = Uinttable_new(nmatches);
  geneend_table = Uinttable_new(nmatches);

  for (i = 0; i < nmatches; i++) {
    goldstandard_gene = IIT_annotation(&restofheader,goldstandard_iit,matches[i],&allocp);
    parse_gene(&goldstandard_donors,&goldstandard_acceptors,&goldstandard_info,goldstandard_gene);
    if (Uintlist_length(goldstandard_donors) != Uintlist_length(goldstandard_acceptors)) {
      fprintf(stderr,"For gene %s, have %d donors and %d acceptors\n",
	      goldstandard_gene,Uintlist_length(goldstandard_donors),Uintlist_length(goldstandard_acceptors));
      abort();
    } else if (Uintlist_length(goldstandard_donors) > 0) {
      Uinttable_put(genestart_table,Uintlist_head(goldstandard_donors),(void *) true);
      Uinttable_put(geneend_table,Uintlist_last_value(goldstandard_acceptors),(void *) true);
      enter_introns(donor_table,goldstandard_donors,goldstandard_acceptors);
      enter_introns(acceptor_table,goldstandard_acceptors,goldstandard_donors);
      Uintlist_free(&goldstandard_acceptors);
      Uintlist_free(&goldstandard_donors);
    }
    if (allocp) {
      FREE(restofheader);
    }
  }

  *ncommon = *npredicted_only = *nincomplete = 0;
  *results_donors = (Intlist_T) NULL;
  *results_acceptors = (Intlist_T) NULL;

  p = predicted_donors;
  q = predicted_acceptors;
  while (p != NULL && q != NULL) {
    donor = Uintlist_head(p);
    acceptor = Uintlist_head(q);
    if ((list = Uinttable_get(donor_table,donor)) != NULL) {
      if (mainp == true && p == predicted_donors && Uinttable_get(genestart_table,donor) == false) {
	*nincomplete += 1;
	*results_donors = Intlist_push(*results_donors,INCOMPLETE);
      } else {
	*ncommon += 1;
	*results_donors = Intlist_push(*results_donors,CORRECT);
      }
      if (dest_presentp(list,acceptor) == false) {
	*npredicted_only += 1;
	*results_acceptors = Intlist_push(*results_acceptors,FALSE);
      } else if (mainp == true && Uintlist_next(q) == NULL && Uinttable_get(geneend_table,acceptor) == false) {
	*nincomplete += 1;
	*results_acceptors = Intlist_push(*results_acceptors,INCOMPLETE);
      } else {
	*ncommon += 1;
	*results_acceptors = Intlist_push(*results_acceptors,CORRECT);
      }
    } else if ((list = Uinttable_get(acceptor_table,acceptor)) != NULL) {
      *npredicted_only += 1;
      *results_donors = Intlist_push(*results_donors,FALSE);

      if (mainp == true && Uintlist_next(q) == NULL && Uinttable_get(geneend_table,acceptor) == false) {
	*nincomplete += 1;
	*results_acceptors = Intlist_push(*results_acceptors,INCOMPLETE);
      } else {
	*ncommon += 1;
	*results_acceptors = Intlist_push(*results_acceptors,CORRECT);
      }
      if (dest_presentp(list,donor) == true) {
	abort();
      }
    } else {
      *npredicted_only += 2;
      *results_donors = Intlist_push(*results_donors,FALSE);
      *results_acceptors = Intlist_push(*results_acceptors,FALSE);
    }
    p = Uintlist_next(p);
    q = Uintlist_next(q);
  }
      

  if ((n = Uinttable_length(donor_table)) > 0) {
    keys = Uinttable_keys(donor_table,/*sortp*/false);
    for (i = 0; i < n; i++) {
      list = Uinttable_get(donor_table,keys[i]);
      Uintlist_free(&list);
    }
    FREE(keys);
  }
  Uinttable_free(&donor_table);

  if ((n = Uinttable_length(acceptor_table)) > 0) {
    keys = Uinttable_keys(acceptor_table,/*sortp*/false);
    for (i = 0; i < n; i++) {
      list = Uinttable_get(acceptor_table,keys[i]);
      Uintlist_free(&list);
    }
    FREE(keys);
  }
  Uinttable_free(&acceptor_table);
  
  Uinttable_free(&genestart_table);
  Uinttable_free(&geneend_table);

  *results_donors = Intlist_reverse(*results_donors);
  *results_acceptors = Intlist_reverse(*results_acceptors);

  return;
}



static void
print_gene_fwd (char *gene, Uintlist_T goldstandard_donors, Uintlist_T goldstandard_acceptors) {
  int exoni;
  Genomicpos_T exonstart, exonend;
  char *p = gene, *q;
  bool bad_donor_p, bad_acceptor_p;

#if 0
  /* Print gene comment line */
  while (*p != '\0' && *p != '\n') { putchar(*p++); }
  if (*p == '\n') { putchar(*p++); }
#else
  /* Skip gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }
#endif

  exoni = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u\n",&exonstart,&exonend) != 2) {
      /* Print line */
      while (*p != '\0' && *p != '\n') { putchar(*p++); }
      if (*p == '\n') { putchar(*p++); }
    } else {
      /* Peek at next line.  This works only if last line is an exon. */
      q = p;
      while (*q != '\0' && *q != '\n') { q++; }
      if (*q == '\n') { q++; }

      bad_donor_p = false;
      if (*q != '\0') {
	/* Compare exonend with goldstandard donor */
	while (goldstandard_donors != NULL && Uintlist_head(goldstandard_donors) < exonend) {
	  printf("# Missed donor: %u\n",Uintlist_head(goldstandard_donors));
	  goldstandard_donors = Uintlist_next(goldstandard_donors);
	}
	if (goldstandard_donors == NULL) {
	  bad_donor_p = true;
	} else if (exonend < Uintlist_head(goldstandard_donors)) {
	  bad_donor_p = true;
	} else {
	  goldstandard_donors = Uintlist_next(goldstandard_donors);
	}
      }

      bad_acceptor_p = false;
      if (exoni > 0) {
	/* Compare exonstart with goldstandard acceptor */
	while (goldstandard_acceptors != NULL && Uintlist_head(goldstandard_acceptors) < exonstart) {
	  printf("# Missed acceptor: %u\n",Uintlist_head(goldstandard_acceptors));
	  goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
	}
	if (goldstandard_acceptors == NULL) {
	  bad_acceptor_p = true;
	} else if (exonstart < Uintlist_head(goldstandard_acceptors)) {
	  bad_acceptor_p = true;
	} else {
	  goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
	}
      }

      /* Print exon */
      while (*p != '\0' && *p != '\n') { putchar(*p++); }
      if (*p == '\n') { p++; }

      if (bad_donor_p == true) {
	printf(" BADDONOR");
      }
      if (bad_acceptor_p == true) {
	printf(" BADACCEPTOR");
      }
      printf("\n");

      exoni++;
    }
  }

  while (goldstandard_donors != NULL) {
    printf("# Missed donor: %u\n",Uintlist_head(goldstandard_donors));
    goldstandard_donors = Uintlist_next(goldstandard_donors);
  }

  while (goldstandard_acceptors != NULL) {
    printf("# Missed acceptor: %u\n",Uintlist_head(goldstandard_acceptors));
    goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
  }

  return;
}


static void
print_gene_rev (char *gene, Uintlist_T goldstandard_donors, Uintlist_T goldstandard_acceptors) {
  int exoni;
  Genomicpos_T exonstart, exonend;
  char *p = gene, *q;
  bool bad_donor_p, bad_acceptor_p;

#if 0
  /* Print gene comment line */
  while (*p != '\0' && *p != '\n') { putchar(*p++); }
  if (*p == '\n') { putchar(*p++); }
#else
  /* Skip gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }
#endif

  exoni = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u\n",&exonstart,&exonend) != 2) {
      while (*p != '\0' && *p != '\n') { putchar(*p++); }
      if (*p == '\n') { putchar(*p++); }
    } else {
      /* Peek at next line.  This works only if last line is an exon. */
      q = p;
      while (*q != '\0' && *q != '\n') { q++; }
      if (*q == '\n') { q++; }

      bad_donor_p = false;
      if (*q != '\0') {
	/* Compare exonend with goldstandard donor */
	while (goldstandard_donors != NULL && Uintlist_head(goldstandard_donors) > exonend) {
	  printf("# Missed donor: %u\n",Uintlist_head(goldstandard_donors));
	  goldstandard_donors = Uintlist_next(goldstandard_donors);
	}
	if (goldstandard_donors == NULL) {
	  bad_donor_p = true;
	} else if (exonend > Uintlist_head(goldstandard_donors)) {
	  bad_donor_p = true;
	} else {
	  goldstandard_donors = Uintlist_next(goldstandard_donors);
	}
      }

      bad_acceptor_p = false;
      if (exoni > 0) {
	/* Compare exonstart with goldstandard acceptor */
	while (goldstandard_acceptors != NULL && Uintlist_head(goldstandard_acceptors) > exonstart) {
	  printf("# Missed acceptor: %u\n",Uintlist_head(goldstandard_acceptors));
	  goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
	}
	if (goldstandard_acceptors == NULL) {
	  bad_acceptor_p = true;
	} else if (exonstart > Uintlist_head(goldstandard_acceptors)) {
	  bad_acceptor_p = true;
	} else {
	  goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
	}
      }

      /* Print exon */
      while (*p != '\0' && *p != '\n') { putchar(*p++); }
      if (*p == '\n') { p++; }

      if (bad_donor_p == true) {
	printf(" BADDONOR");
      }
      if (bad_acceptor_p == true) {
	printf(" BADACCEPTOR");
      }
      printf("\n");

      exoni++;
    }
  }

  while (goldstandard_donors != NULL) {
    printf("# Missed donor: %u\n",Uintlist_head(goldstandard_donors));
    goldstandard_donors = Uintlist_next(goldstandard_donors);
  }

  while (goldstandard_acceptors != NULL) {
    printf("# Missed acceptor: %u\n",Uintlist_head(goldstandard_acceptors));
    goldstandard_acceptors = Uintlist_next(goldstandard_acceptors);
  }

  return;
}


static void
print_gene_nocomment (char *gene) {
  char *p = gene;

  /* Skip gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }

  printf("%s",p);
  return;
}


static void
remove_from_table (Uinttable_T donortable, Uintlist_T donors, Uintlist_T acceptors) {
  Uintlist_T list, p, q;
  Genomicpos_T donor, acceptor;

  p = donors;
  q = acceptors;
  while (p != NULL && q != NULL) {
    donor = Uintlist_head(p);
    acceptor = Uintlist_head(q);
    list = Uinttable_get(donortable,donor);
    Uinttable_put(donortable,donor,Uintlist_remove(list,acceptor));
    p = Uintlist_next(p);
    q = Uintlist_next(q);
  }

  return;
}


static bool
has_intron_p (Uintlist_T donors, Uintlist_T acceptors,
	      Genomicpos_T donor, Genomicpos_T acceptor) {
  Uintlist_T p, q;

  p = donors;
  q = acceptors;
  while (p != NULL && q != NULL) {
    if (Uintlist_head(p) == donor && Uintlist_head(q) == acceptor) {
      return true;
    }
    p = Uintlist_next(p);
    q = Uintlist_next(q);
  }

  return false;
}


static void
evaluate_chromosome (IIT_T predicted_iit, IIT_T goldstandard_iit, char *chr,
		     bool originalp) {
  int nmissed = 0;
  int main_nperfect_fwd = 0, main_nerror_fwd = 0, main_nnomatch_fwd = 0;
  int main_nperfect_rev = 0, main_nerror_rev = 0, main_nnomatch_rev = 0;
  int alt_nperfect_fwd = 0, alt_nerror_fwd = 0, alt_nnomatch_fwd = 0;
  int alt_nperfect_rev = 0, alt_nerror_rev = 0, alt_nnomatch_rev = 0;

  Uinttable_T donortable;
  Genomicpos_T predicted_low, predicted_high, *donors, donor, acceptor;
  int *predicted_matches, predicted_nmatches, n, i, j;
  int *goldstandard_matches, goldstandard_nmatches, bestindex;
  int ncommon = 0, npredicted_only = 0, ngoldstandard_only = 0, nincomplete = 0;
  int goldstandard_divno;
  int sign;
  Interval_T predicted_interval;
  char *predicted_label, *predicted_gene, *goldstandard_label, *goldstandard_gene, *restofheader;
  Uintlist_T predicted_donors, predicted_acceptors, goldstandard_donors, goldstandard_acceptors;
  List_T predicted_info, goldstandard_info;
  Uintlist_T list, p, q;
  Intlist_T results_donors, results_acceptors;
  bool mainp, allocp, allocp0, allocp1, allocp2;

  if ((goldstandard_divno = IIT_divint(goldstandard_iit,/*divstring*/chr)) <= 0) {
    fprintf(stderr,"Chromosome %s in predicted IIT not found in goldstandard IIT\n",chr);
    return;
  } else {
    goldstandard_matches = IIT_get(&goldstandard_nmatches,goldstandard_iit,/*divstring*/chr,
				   /*x*/0,/*y*/-1U,/*sortp*/false);
    donortable = Uinttable_new(goldstandard_nmatches*10);

    for (i = 0; i < goldstandard_nmatches; i++) {
      goldstandard_gene = IIT_annotation(&restofheader,goldstandard_iit,goldstandard_matches[i],&allocp);
      parse_gene(&goldstandard_donors,&goldstandard_acceptors,&goldstandard_info,goldstandard_gene);
      p = goldstandard_donors;
      q = goldstandard_acceptors;
      while (p != NULL && q != NULL) {
	donor = Uintlist_head(p);
	acceptor = Uintlist_head(q);
	list = Uinttable_get(donortable,donor);
	if (Uintlist_find(list,acceptor) == false) {
	  Uinttable_put(donortable,donor,Uintlist_push(list,acceptor));
	}
	p = Uintlist_next(p);
	q = Uintlist_next(q);
      }
      if (allocp) FREE(restofheader);
      Uintlist_free(&goldstandard_donors);
      Uintlist_free(&goldstandard_acceptors);
    }
    FREE(goldstandard_matches);
  }

  predicted_matches = IIT_get(&predicted_nmatches,predicted_iit,/*divstring*/chr,/*x*/0,/*y*/-1U,/*sortp*/false);
  for (i = 0; i < predicted_nmatches; i++) {
    predicted_label = IIT_label(predicted_iit,predicted_matches[i],&allocp0);
    predicted_gene = IIT_annotation(&restofheader,predicted_iit,predicted_matches[i],&allocp1);
    mainp = parse_gene(&predicted_donors,&predicted_acceptors,&predicted_info,predicted_gene);

    remove_from_table(donortable,predicted_donors,predicted_acceptors);

    predicted_interval = IIT_interval(predicted_iit,predicted_matches[i]);
    sign = Interval_sign(predicted_interval);
    predicted_low = Interval_low(predicted_interval);
    predicted_high = Interval_high(predicted_interval);

    goldstandard_matches = IIT_get_signed_with_divno(&goldstandard_nmatches,goldstandard_iit,goldstandard_divno,
						     predicted_low,predicted_high,/*sortp*/false,sign);

    if (sign >= 0) {
      printf(">%s %s:%u..%u ",predicted_label,chr,predicted_low,predicted_high);
    } else {
      printf(">%s %s:%u..%u ",predicted_label,chr,predicted_high,predicted_low);
    }

    if (goldstandard_nmatches == 0) {
      if (mainp == true) {
	if (sign > 0) main_nnomatch_fwd++; else main_nnomatch_rev++;
      } else {
	if (sign > 0) alt_nnomatch_fwd++; else alt_nnomatch_rev++;
      }
      printf("NOMATCH\n");
      printf("X 0 %d 0\n",Uintlist_length(predicted_donors)+Uintlist_length(predicted_acceptors));
      if (originalp) {
	print_gene_nocomment(predicted_gene);
      } else if (sign > 0) {
	print_onegene_fwd(predicted_donors,predicted_acceptors,predicted_info,/*tag*/"false");
      } else {
	print_onegene_rev(predicted_donors,predicted_acceptors,predicted_info,/*tag*/"false");
      }

    } else if (1) {
      evaluate_splices(&ncommon,&npredicted_only,&nincomplete,
		       &results_donors,&results_acceptors,
		       predicted_donors,predicted_acceptors,
		       goldstandard_matches,goldstandard_nmatches,goldstandard_iit,
		       mainp);
      if (npredicted_only == 0 && nincomplete == 0) {
	if (mainp == true) {
	  if (sign > 0) main_nperfect_fwd++; else main_nperfect_rev++;
	} else {
	  if (sign > 0) alt_nperfect_fwd++; else alt_nperfect_rev++;
	}
	printf("PERFECT\n");
	printf("%d %d %d",ncommon,npredicted_only,nincomplete);
	printf("\n");
	if (originalp) {
	  print_gene_nocomment(predicted_gene);
	} else if (sign > 0) {
	  print_onegene_fwd(predicted_donors,predicted_acceptors,predicted_info,/*tag*/NULL);
	} else {
	  print_onegene_rev(predicted_donors,predicted_acceptors,predicted_info,/*tag*/NULL);
	}

      } else {
	if (mainp == true) {
	  if (sign > 0) main_nerror_fwd++; else main_nerror_rev++;
	} else {
	  if (sign > 0) alt_nerror_fwd++; else alt_nerror_rev++;
	}
	printf("ERROR\n");
	printf("%d %d %d",ncommon,npredicted_only,nincomplete);
	printf("\n");
	if (originalp) {
	  print_gene_nocomment(predicted_gene);
	} else if (sign > 0) {
	  print_onegene_wresults_fwd(predicted_donors,predicted_acceptors,predicted_info,
				     results_donors,results_acceptors,/*tag*/NULL);
	} else {
	  print_onegene_wresults_rev(predicted_donors,predicted_acceptors,predicted_info,
				     results_donors,results_acceptors,/*tag*/NULL);
	}
      }

      Intlist_free(&results_donors);
      Intlist_free(&results_acceptors);

    } else {
      bestindex = evaluate_genes(&ncommon,&npredicted_only,&ngoldstandard_only,
				 &goldstandard_donors,&goldstandard_acceptors,
				 predicted_donors,predicted_acceptors,
				 goldstandard_matches,goldstandard_nmatches,goldstandard_iit,
				 sign,mainp);
      goldstandard_label = IIT_label(goldstandard_iit,bestindex,&allocp2);
      if (mainp == false) {
	ngoldstandard_only = 0;
      }

      if (npredicted_only == 0 && ngoldstandard_only == 0) {
	if (mainp == true) {
	  if (sign > 0) main_nperfect_fwd++; else main_nperfect_rev++;
	} else {
	  if (sign > 0) alt_nperfect_fwd++; else alt_nperfect_rev++;
	}
	printf("PERFECT\n");
	printf("%s %d %d",goldstandard_label,ncommon,npredicted_only);
	if (mainp == true) {
	  printf(" %d",ngoldstandard_only);
	}
	printf("\n");
	if (originalp) {
	  print_gene_nocomment(predicted_gene);
	} else if (sign > 0) {
	  print_onegene_fwd(predicted_donors,predicted_acceptors,predicted_info,/*tag*/NULL);
	} else {
	  print_onegene_rev(predicted_donors,predicted_acceptors,predicted_info,/*tag*/NULL);
	}

      } else {
	if (mainp == true) {
	  if (sign > 0) main_nerror_fwd++; else main_nerror_rev++;
	} else {
	  if (sign > 0) alt_nerror_fwd++; else alt_nerror_rev++;
	}
	printf("ERROR\n");
	printf("%s %d %d",goldstandard_label,ncommon,npredicted_only);
	if (mainp == true) {
	  printf(" %d",ngoldstandard_only);
	}
	printf("\n");
	if (originalp) {
	  print_gene_nocomment(predicted_gene);
	} else if (sign > 0) {
	  print_comparison_fwd(predicted_donors,predicted_acceptors,goldstandard_donors,goldstandard_acceptors);
	} else {
	  print_comparison_rev(predicted_donors,predicted_acceptors,goldstandard_donors,goldstandard_acceptors);
	}
      }

      Uintlist_free(&goldstandard_donors);
      Uintlist_free(&goldstandard_acceptors);
      if (allocp2) {
	FREE(goldstandard_label);
      }
    }

#if 0
    if (sign > 0) {
      print_gene_fwd(predicted_gene,goldstandard_donors,goldstandard_acceptors);
    } else {
      print_gene_rev(predicted_gene,goldstandard_donors,goldstandard_acceptors);
    }
#endif

    Uintlist_free(&predicted_donors);
    Uintlist_free(&predicted_acceptors);
    if (allocp1) {
      FREE(restofheader);
    }
    if (allocp0) {
      FREE(predicted_label);
    }
  }

  FREE(predicted_matches);


  if ((n = Uinttable_length(donortable)) > 0) {
    donors = Uinttable_keys(donortable,/*sortp*/true);
    for (i = 0; i < n; i++) {
      donor = donors[i];
      list = Uinttable_get(donortable,donor);
      for (p = list; p != NULL; p = Uintlist_next(p)) {
	acceptor = Uintlist_head(p);
	printf(">Missed %s:%u..%u\n",chr,donor,acceptor);
	nmissed++;
	if (donor < acceptor) {
	  goldstandard_matches = IIT_get(&goldstandard_nmatches,goldstandard_iit,/*divstring*/chr,
					 /*x*/donor,/*y*/acceptor,/*sortp*/false);
	} else {
	  goldstandard_matches = IIT_get(&goldstandard_nmatches,goldstandard_iit,/*divstring*/chr,
					 /*x*/acceptor,/*y*/donor,/*sortp*/false);
	}
	for (j = 0; j < goldstandard_nmatches; j++) {
	  goldstandard_gene = IIT_annotation(&restofheader,goldstandard_iit,goldstandard_matches[j],&allocp);
	  parse_gene(&goldstandard_donors,&goldstandard_acceptors,&goldstandard_info,goldstandard_gene);
	  if (has_intron_p(goldstandard_donors,goldstandard_acceptors,donor,acceptor) == true) {
	    goldstandard_label = IIT_label(goldstandard_iit,goldstandard_matches[j],&allocp2);
	    printf("%s\n",goldstandard_label);
	    if (allocp2) FREE(goldstandard_label);
	  }
	  Uintlist_free(&goldstandard_donors);
	  Uintlist_free(&goldstandard_acceptors);
	  if (allocp) FREE(restofheader);
	}
      }
      Uintlist_free(&list);
    }
    FREE(donors);
  }
  Uinttable_free(&donortable);

  global_main_nperfect_fwd += main_nperfect_fwd;
  global_main_nerror_fwd += main_nerror_fwd;
  global_main_nnomatch_fwd += main_nnomatch_fwd;
  global_main_nperfect_rev += main_nperfect_rev;
  global_main_nerror_rev += main_nerror_rev;
  global_main_nnomatch_rev += main_nnomatch_rev;

  global_alt_nperfect_fwd += alt_nperfect_fwd;
  global_alt_nerror_fwd += alt_nerror_fwd;
  global_alt_nnomatch_fwd += alt_nnomatch_fwd;
  global_alt_nperfect_rev += alt_nperfect_rev;
  global_alt_nerror_rev += alt_nerror_rev;
  global_alt_nnomatch_rev += alt_nnomatch_rev;

  fprintf(stderr,"Chromosome +%s, splices missed: %d\n",chr,nmissed);
  fprintf(stderr,"Chromosome +%s, main: %d perfect, %d error, %d nomatch\n",
	  chr,main_nperfect_fwd,main_nerror_fwd,main_nnomatch_fwd);
  fprintf(stderr,"Chromosome -%s, main: %d perfect, %d error, %d nomatch\n",
	  chr,main_nperfect_rev,main_nerror_rev,main_nnomatch_rev);
  fprintf(stderr,"Chromosome +%s, alt: %d perfect, %d error, %d nomatch\n",
	  chr,alt_nperfect_fwd,alt_nerror_fwd,alt_nnomatch_fwd);
  fprintf(stderr,"Chromosome -%s, alt: %d perfect, %d error, %d nomatch\n",
	  chr,alt_nperfect_rev,alt_nerror_rev,alt_nnomatch_rev);

  return;
}


int
main (int argc, char *argv[]) {
  char *user_genomedir = NULL, *dbroot = NULL, *genomesubdir;
  char *mapdir, *goldstandard = NULL, *iitfile;
  char *fileroot, *dbversion;

  Genomicpos_T genomicstart, genomiclength, chroffset, chrlength, chrstart, chrend;
  char *chr = NULL;
  bool revcomp;

  IIT_T goldstandard_iit, predicted_iit;
  int predicted_divno;

  int coordi, predi;
  int total;

  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"Dd:m:O")) != -1) {
    switch (c) {
    case 'D': user_genomedir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'm': goldstandard = optarg; break;
    case 'O': originalp = true; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);


  if (dbroot != NULL && goldstandard != NULL) {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(goldstandard)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,goldstandard);
    if ((goldstandard_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				     /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Could not open known sites IIT file %s\n",goldstandard);
      exit(9);
    }
    FREE(iitfile);

    predi = 1;
    coordi = 2;

  } else {

    if ((goldstandard_iit = IIT_read(argv[1],/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				     /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/false)) == NULL) {
      fprintf(stderr,"Could not open IIT file %s\n",argv[1]);
      exit(9);
    }

    predi = 2;
    coordi = 3;
  }


  if ((predicted_iit = IIT_read(argv[predi],/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/false)) == NULL) {
    fprintf(stderr,"Could not open IIT file %s\n",argv[predi]);
    exit(9);
  }

  if (argc > coordi) {
    Parserange_universal(&chr,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			 &chroffset,&chrlength,argv[coordi],genomesubdir,fileroot);
  }


  if (chr != NULL) {
    evaluate_chromosome(predicted_iit,goldstandard_iit,chr,originalp);
  } else {
    for (predicted_divno = 1; predicted_divno < IIT_ndivs(predicted_iit); predicted_divno++) {
      chr = IIT_divstring(predicted_iit,predicted_divno);
      evaluate_chromosome(predicted_iit,goldstandard_iit,chr,originalp);
    }

    total = global_main_nperfect_fwd + global_main_nerror_fwd + global_main_nnomatch_fwd;
    fprintf(stderr,"Main, total fwd: %d perfect (%.1f%%), %d error (%.1f%%), %d nomatch (%.1f%%)\n",
	    global_main_nperfect_fwd,100.0 * (double) global_main_nperfect_fwd/(double) total,
	    global_main_nerror_fwd,100.0 * (double) global_main_nerror_fwd/(double) total,
	    global_main_nnomatch_fwd,100.0 * (double) global_main_nnomatch_fwd/(double) total);
    total = global_main_nperfect_rev + global_main_nerror_rev + global_main_nnomatch_rev;
    fprintf(stderr,"Main, total rev: %d perfect (%.1f%%), %d error (%.1f%%), %d nomatch (%.1f%%)\n",
	    global_main_nperfect_rev,100.0 * (double) global_main_nperfect_rev/(double) total,
	    global_main_nerror_rev,100.0 * (double) global_main_nerror_rev/(double) total,
	    global_main_nnomatch_rev,100.0 * (double) global_main_nnomatch_rev/(double) total);

    total = global_alt_nperfect_fwd + global_alt_nerror_fwd + global_alt_nnomatch_fwd;
    fprintf(stderr,"Alt, total fwd: %d perfect (%.1f%%), %d error (%.1f%%), %d nomatch (%.1f%%)\n",
	    global_alt_nperfect_fwd,100.0 * (double) global_alt_nperfect_fwd/(double) total,
	    global_alt_nerror_fwd,100.0 * (double) global_alt_nerror_fwd/(double) total,
	    global_alt_nnomatch_fwd,100.0 * (double) global_alt_nnomatch_fwd/(double) total);
    total = global_alt_nperfect_rev + global_alt_nerror_rev + global_alt_nnomatch_rev;
    fprintf(stderr,"Alt, total rev: %d perfect (%.1f%%), %d error (%.1f%%), %d nomatch (%.1f%%)\n",
	    global_alt_nperfect_rev,100.0 * (double) global_alt_nperfect_rev/(double) total,
	    global_alt_nerror_rev,100.0 * (double) global_alt_nerror_rev/(double) total,
	    global_alt_nnomatch_rev,100.0 * (double) global_alt_nnomatch_rev/(double) total);
  }
    
  IIT_free(&predicted_iit);
  IIT_free(&goldstandard_iit);

  return 0;
}

