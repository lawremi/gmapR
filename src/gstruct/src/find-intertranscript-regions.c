#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <ctype.h>		/* For isspace */
#include <math.h>		/* For qsort */

#include "mem.h"
#include "bool.h"
#include "types.h"
#include "uintlist.h"
#include "uinttable.h"
#include "list.h"
#include "interval.h"
#include "iit-read.h"
#include "table.h"
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif




typedef unsigned int Chrpos_T;


static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static char *knowngenes_iitfile = NULL;

static char *ignore_filename = NULL;
static Table_T ignore_table = NULL;


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"genome", required_argument, 0, 'd'}, /* dbroot */

  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */

  {"ignore", required_argument, 0, 'I'},  /* ignore_filename */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"find-intertranscript-regions\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}



/* Modified from gstruct.c */
static char *
parse_gene (Chrpos_T *first_donor, Chrpos_T *last_donor,
	    Chrpos_T *first_acceptor, Chrpos_T *last_acceptor,
	    char *gene) {
  char *genename;
  Chrpos_T exonstart, exonend;
  char *annot, *p, *q;
  Uintlist_T exonstarts = NULL, exonends = NULL;
  int gene_namelength, k;

  /* Process gene comment line */
  p = annot = gene;
  while (*p != '\0' && *p != '\n' && !isspace(*p)) {
    p++;
  }

  gene_namelength = p - annot;
  genename = (char *) CALLOC(gene_namelength + 1,sizeof(char));
  strncpy(genename,annot,gene_namelength);
  for (q = genename, k = 0; k < gene_namelength; q++, k++) {
    if (*q == '-') {
      *q = '.';
    }
  }

  /* Read to end of gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }

  /* Process exon lines */
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) == 2) {
      exonstarts = Uintlist_push(exonstarts,exonstart);
      exonends = Uintlist_push(exonends,exonend);
    }

    /* Skip part just read */
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && *p != '\n' && isspace(*p)) { p++; }

    /* Read to end of line */
    if (*p == '\0' || *p == '\n') {
      /* info = (char *) CALLOC(1,sizeof(char)); */
      /* info[0] = '\0'; */
    } else {
      /* q = p; */
      while (*p != '\0' && *p != '\n') { p++; }

      /* info = (char *) CALLOC(p - q + 1,sizeof(char)); */
      /* strncpy(info,q,p-q); */
      if (*p == '\n') { p++; }
    }

    /* *info_list = List_push(*info_list,(void *) info); */
  }

  if (Uintlist_length(exonstarts) < 2 || Uintlist_length(exonends) < 2) {
    Uintlist_free(&exonstarts);
    Uintlist_free(&exonends);
    FREE(genename);
    return (char *) NULL;

  } else {
    *last_acceptor = Uintlist_head(exonstarts);
    exonstarts = Uintlist_reverse(exonstarts);
    *first_acceptor = Uintlist_head(Uintlist_next(exonstarts));

    *last_donor = Uintlist_head(Uintlist_next(exonends));
    exonends = Uintlist_reverse(exonends);
    *first_donor = Uintlist_head(exonends);
    Uintlist_free(&exonstarts);
    Uintlist_free(&exonends);
    return genename;
  }
}


static int
gene_nexons (char *gene) {
  int nexons = 0;
  Chrpos_T exonstart, exonend;
  char *p;

  p = gene;

  /* Read to end of gene comment line */
  while (*p != '\0' && *p != '\n') { p++; }
  if (*p != '\0') { p++; }

  /* Process exon lines */
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) == 2) {
      nexons++;
    }

    /* Skip part just read */
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && isspace(*p)) { p++; }
    while (*p != '\0' && !isspace(*p)) { p++; }
    while (*p != '\0' && *p != '\n' && isspace(*p)) { p++; }

    /* Read to end of line */
    if (*p == '\0' || *p == '\n') {
      /* info = (char *) CALLOC(1,sizeof(char)); */
      /* info[0] = '\0'; */
    } else {
      /* q = p; */
      while (*p != '\0' && *p != '\n') { p++; }

      /* info = (char *) CALLOC(p - q + 1,sizeof(char)); */
      /* strncpy(info,q,p-q); */
      if (*p == '\n') { p++; }
    }

    /* *info_list = List_push(*info_list,(void *) info); */
  }

  return nexons;
}



static int
string_cmp (const void *a, const void *b) {
  char *x = * (char **) a;
  char *y = * (char **) b;

  return strcmp(x,y);
}


static List_T
remove_duplicates (List_T names) {
  List_T unique = NULL;
  char **array;
  bool *duplicatep;
  int n, i;

  array = (char **) List_to_array_n(&n,names);
  duplicatep = (bool *) CALLOC(n,sizeof(bool));

  qsort(array,n,sizeof(char *),string_cmp);
  unique = List_push(NULL,(void *) array[0]);
  for (i = 1; i < n; i++) {
    if (strcmp(array[i],array[i-1])) {
      unique = List_push(unique,(void *) array[i]);
    } else {
      duplicatep[i] = true;
    }
  }

  for (i = 1; i < n; i++) {
    if (duplicatep[i] == true) {
      free(array[i]);		/* created using strdup */
    }
  }
  FREE(duplicatep);
  FREE(array);

  List_free(&names);
  return List_reverse(unique);
}

static bool
overlapp (List_T names1, List_T names2) {
  List_T p, q;

  for (p = names1; p != NULL; p = List_next(p)) {
    for (q = names2; q != NULL; q = List_next(q)) {
      if (!strcmp(List_head(p),List_head(q))) {
	return true;
      }
    }
  }

  return false;
}


static void
print_list_as_string (List_T list) {
  List_T p;

  if (list != NULL) {
    printf("%s",(char *) List_head(list));

    for (p = List_next(list); p != NULL; p = List_next(p)) {
      printf(",%s",(char *) List_head(p));
    }
  }

  return;
}


static void
table_stringlist_gc (Uinttable_T *table) {
  Chrpos_T *keys;
  int n, i;
  char *string;
  List_T strings, p;

  keys = Uinttable_keys(*table,/*sortp*/false);
  n = Uinttable_length(*table);
  for (i = 0; i < n; i++) {
    strings = (List_T) Uinttable_get(*table,keys[i]);
    for (p = strings; p != NULL; p = List_next(p)) {
      string = (char *) List_head(p);
      free(string);		/* created using strdup */
    }
    List_free(&strings);
  }
  FREE(keys);
  Uinttable_free(&(*table));
  return;
}


/* Previously called
   IIT_contained_by_region_with_divno_signed(knowngenes_iit,divno,low,high,sign),
   but that allows single exon genes to disrupt a readthrough */

static bool
multiexon_gene_contained_by_region_p (IIT_T knowngenes_iit, int divno, Chrpos_T low, Chrpos_T high, int sign) {
  int *matches;
  int nmatches, i;
  Interval_T interval;
  char *gene, *restofheader;
  bool alloc_annot_p;

  matches = IIT_get_signed_with_divno(&nmatches,knowngenes_iit,divno,low,high,/*sortp*/false,sign);
  for (i = 0; i < nmatches; i++) {
    interval = IIT_interval(knowngenes_iit,matches[i]);
    if (low < Interval_low(interval) && Interval_high(interval) < high) {
      gene = IIT_annotation(&restofheader,knowngenes_iit,matches[i],&alloc_annot_p);
      if (gene_nexons(gene) > 1) {
	if (alloc_annot_p) FREE(restofheader);
	FREE(matches);
	return true;
      }
      if (alloc_annot_p) FREE(restofheader);
    }
  }
  FREE(matches);
  return false;
}



typedef struct Chrpos_obj_T *Chrpos_obj_T;
struct Chrpos_obj_T {
  Chrpos_T pos;
};

static Chrpos_obj_T
Chrpos_obj_new (Chrpos_T pos) {
  Chrpos_obj_T new = MALLOC(sizeof(*new));

  new->pos = pos;
  return new;
}

static void
Chrpos_obj_free (Chrpos_obj_T *old) {
  FREE(*old);
  return;
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
  char *genomesubdir, *fileroot, *mapdir, *iitfile;
  IIT_T knowngenes_iit;
  int len;

  FILE *fp;
  int divno;
  char *chr, *acc, *gene, *restofheader, *genename;
  int *matches;
  int nmatches, n, i, j, k;
  Uinttable_T plus_donor_accs, plus_acceptor_accs, minus_donor_accs, minus_acceptor_accs,
    plus_donor_genes, plus_acceptor_genes, minus_donor_genes, minus_acceptor_genes,
    plus_outside, minus_outside;
  Interval_T interval;
  bool alloc_label_p, alloc_annot_p;
  Chrpos_T first_donor, last_donor, first_acceptor, last_acceptor;
  Chrpos_obj_T outside;
  int sign;

  char **accessions;
  Chrpos_T *keys, *donors, *acceptors, donor, acceptor;
  int ndonors, nacceptors;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:M:m:I:?",
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

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'find-intertranscript-regions --help'",long_name);
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
      knowngenes_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(knowngenes_iitfile,optarg);
      if ((len = strlen(knowngenes_iitfile)) > 4 && strcmp(&(knowngenes_iitfile[len-4]),".iit") == 0) {
	knowngenes_iitfile[len-4] = '\0';
      }
      break;

    case 'I': ignore_filename = optarg; break;

    case '?': /* print_program_usage(); */ exit(0);
    default:
      fprintf(stderr,"Cannot handle flag %c\n",opt);
      exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);


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


  /* Access gene_knowngenes_iit */
  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  if (knowngenes_iitfile == NULL) {
    /* fprintf(stderr,"Must specify genes iit with the -m flag\n"); */
    knowngenes_iit = (IIT_T) NULL;
  } else {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(knowngenes_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,knowngenes_iitfile);
    if ((knowngenes_iit = IIT_read(iitfile,/*name*/knowngenes_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",knowngenes_iitfile,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",knowngenes_iitfile);
      fprintf(stderr,"using the -M flag to gmap.\n");
      exit(9);
    }
    FREE(mapdir);
    FREE(knowngenes_iitfile);
    FREE(iitfile);
  }

  FREE(dbversion);
  FREE(fileroot);
  FREE(dbroot);
  FREE(genomesubdir);


  for (divno = 1; divno < IIT_ndivs(knowngenes_iit); divno++) {
    plus_donor_accs = Uinttable_new(1000);
    plus_acceptor_accs = Uinttable_new(1000);
    minus_donor_accs = Uinttable_new(1000);
    minus_acceptor_accs = Uinttable_new(1000);

    plus_donor_genes = Uinttable_new(1000);
    plus_acceptor_genes = Uinttable_new(1000);
    minus_donor_genes = Uinttable_new(1000);
    minus_acceptor_genes = Uinttable_new(1000);
    
    plus_outside = Uinttable_new(1000);
    minus_outside = Uinttable_new(1000);


    chr = IIT_divstring(knowngenes_iit,divno);
    matches = IIT_get(&nmatches,knowngenes_iit,chr,/*x*/0,/*y*/-1U,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      acc = IIT_label(knowngenes_iit,matches[i],&alloc_label_p);

      if (ignore_table != NULL && Table_get(ignore_table,acc) != NULL) {
	/* Skip readthrough acc */
      } else {
	gene = IIT_annotation(&restofheader,knowngenes_iit,matches[i],&alloc_annot_p);
	interval = IIT_interval(knowngenes_iit,matches[i]);
	if ((genename = parse_gene(&first_donor,&last_donor,&first_acceptor,&last_acceptor,gene)) != NULL) {
	  if ((sign = Interval_sign(interval)) > 0) {
	    Uinttable_put(plus_donor_accs,last_donor,
			  (void *) List_push((List_T) Uinttable_get(plus_donor_accs,last_donor),(void *) strdup(acc)));
	    Uinttable_put(plus_donor_genes,last_donor,
			  (void *) List_push((List_T) Uinttable_get(plus_donor_genes,last_donor),(void *) strdup(genename)));
	    Uinttable_put(plus_acceptor_accs,first_acceptor,
			  (void *) List_push((List_T) Uinttable_get(plus_acceptor_accs,first_acceptor),(void *) strdup(acc)));
	    Uinttable_put(plus_acceptor_genes,first_acceptor,
			  (void *) List_push((List_T) Uinttable_get(plus_acceptor_genes,first_acceptor),(void *) strdup(genename)));
	    if ((outside = (Chrpos_obj_T) Uinttable_get(plus_outside,last_donor)) == NULL) {
	      Uinttable_put(plus_outside,last_donor,(void *) Chrpos_obj_new(first_donor));
	    } else if (first_donor < outside->pos) {
	      outside->pos = first_donor;
	    }
	    if ((outside = (Chrpos_obj_T) Uinttable_get(plus_outside,first_acceptor)) == NULL) {
	      Uinttable_put(plus_outside,first_acceptor,(void *) Chrpos_obj_new(last_acceptor));
	    } else if (last_acceptor > outside->pos) {
	      outside->pos = last_acceptor;
	    }

	  } else if (sign < 0) {
	    Uinttable_put(minus_donor_accs,last_donor,
			  (void *) List_push((List_T) Uinttable_get(minus_donor_accs,last_donor),(void *) strdup(acc)));
	    Uinttable_put(minus_donor_genes,last_donor,
			  (void *) List_push((List_T) Uinttable_get(minus_donor_genes,last_donor),(void *) strdup(genename)));
	    Uinttable_put(minus_acceptor_accs,first_acceptor,
			  (void *) List_push((List_T) Uinttable_get(minus_acceptor_accs,first_acceptor),(void *) strdup(acc)));
	    Uinttable_put(minus_acceptor_genes,first_acceptor,
			  (void *) List_push((List_T) Uinttable_get(minus_acceptor_genes,first_acceptor),(void *) strdup(genename)));
	    if ((outside = (Chrpos_obj_T) Uinttable_get(minus_outside,last_donor)) == NULL) {
	      Uinttable_put(minus_outside,last_donor,(void *) Chrpos_obj_new(first_donor));
	    } else if (first_donor > outside->pos) {
	      outside->pos = first_donor;
	    }
	    if ((outside = Uinttable_get(minus_outside,first_acceptor)) == NULL) {
	      Uinttable_put(minus_outside,first_acceptor,(void *) Chrpos_obj_new(last_acceptor));
	    } else if (last_acceptor < outside->pos) {
	      outside->pos = first_acceptor;
	    }

	  } else {
	    fprintf(stderr,"Interval sign is 0\n");
	    abort();
	  }
	  FREE(genename);
	}

	if (alloc_annot_p) FREE(restofheader);
      }

      if (alloc_label_p) FREE(acc);
    }
    FREE(matches);


    /* Plus direction */
    donors = Uinttable_keys(plus_donor_accs,/*sortp*/true);
    acceptors = Uinttable_keys(plus_acceptor_accs,/*sortp*/true);
    ndonors = Uinttable_length(plus_donor_accs);
    nacceptors = Uinttable_length(plus_acceptor_accs);

    for (i = 0; i < ndonors; i++) {
      Uinttable_put(plus_donor_genes,donors[i],
		    (void *) remove_duplicates(Uinttable_get(plus_donor_genes,donors[i])));
    }
    for (j = 0; j < nacceptors; j++) {
      Uinttable_put(plus_acceptor_genes,acceptors[j],
		    (void *) remove_duplicates(Uinttable_get(plus_acceptor_genes,acceptors[j])));
    }

    i = j = 0;
    while (i < ndonors) {
      donor = donors[i];
      while (j < nacceptors && acceptors[j] <= donor) {
	j++;
      }
      k = j;
      while (k < nacceptors &&
	     multiexon_gene_contained_by_region_p(knowngenes_iit,divno,donor,acceptors[k],/*sign*/+1) == false) {
	acceptor = acceptors[k];
	if (overlapp(Uinttable_get(plus_donor_accs,donor),Uinttable_get(plus_acceptor_accs,acceptor)) == false &&
	    overlapp(Uinttable_get(plus_donor_genes,donor),Uinttable_get(plus_acceptor_genes,acceptor)) == false) {
	  print_list_as_string((List_T) Uinttable_get(plus_donor_genes,donor));
	  printf("..");
	  print_list_as_string((List_T) Uinttable_get(plus_acceptor_genes,acceptor));
	  printf("\t");

	  printf("%s:%u..%u\t",chr,donor,acceptor);

	  printf("%s:",chr);
	  outside = (Chrpos_obj_T) Uinttable_get(plus_outside,donor);
	  printf("%u..",outside->pos);
	  outside = (Chrpos_obj_T) Uinttable_get(plus_outside,acceptor);
	  printf("%u",outside->pos);
	  printf("\t");

	  print_list_as_string((List_T) Uinttable_get(plus_donor_accs,donor));
	  printf("\t");
	  print_list_as_string((List_T) Uinttable_get(plus_acceptor_accs,acceptor));
	  printf("\n");
	}
	k++;
      }
      i++;
    }
    FREE(acceptors);
    FREE(donors);


    /* Minus direction */
    donors = Uinttable_keys(minus_donor_accs,/*sortp*/true);
    acceptors = Uinttable_keys(minus_acceptor_accs,/*sortp*/true);
    ndonors = Uinttable_length(minus_donor_accs);
    nacceptors = Uinttable_length(minus_acceptor_accs);

    for (i = ndonors - 1; i >= 0; i--) {
      Uinttable_put(minus_donor_genes,donors[i],
		    (void *) remove_duplicates(Uinttable_get(minus_donor_genes,donors[i])));
    }
    for (j = nacceptors - 1; j >= 0; j--) {
      Uinttable_put(minus_acceptor_genes,acceptors[j],
		    (void *) remove_duplicates(Uinttable_get(minus_acceptor_genes,acceptors[j])));
    }

    i = ndonors - 1;
    j = nacceptors - 1;
    while (i >= 0) {
      donor = donors[i];
      while (j >= 0 && acceptors[j] >= donor) {
	j--;
      }
      k = j;
      while (k >= 0 &&
	     multiexon_gene_contained_by_region_p(knowngenes_iit,divno,acceptors[k],donor,/*sign*/-1) == false) {
	acceptor = acceptors[k];
	if (overlapp(Uinttable_get(minus_donor_accs,donor),Uinttable_get(minus_acceptor_accs,acceptor)) == false &&
	    overlapp(Uinttable_get(minus_donor_genes,donor),Uinttable_get(minus_acceptor_genes,acceptor)) == false) {
	  print_list_as_string((List_T) Uinttable_get(minus_donor_genes,donor));
	  printf("..");
	  print_list_as_string((List_T) Uinttable_get(minus_acceptor_genes,acceptor));
	  printf("\t");

	  printf("%s:%u..%u\t",chr,donor,acceptor);

	  printf("%s:",chr);
	  outside = (Chrpos_obj_T) Uinttable_get(minus_outside,donor);
	  printf("%u..",outside->pos);
	  outside = (Chrpos_obj_T) Uinttable_get(minus_outside,acceptor);
	  printf("%u",outside->pos);
	  printf("\t");

	  print_list_as_string((List_T) Uinttable_get(minus_donor_accs,donor));
	  printf("\t");
	  print_list_as_string((List_T) Uinttable_get(minus_acceptor_accs,acceptor));
	  printf("\n");
	}
	k--;
      }
      i--;
    }
    FREE(acceptors);
    FREE(donors);

    keys = Uinttable_keys(minus_outside,/*sortp*/false);
    n = Uinttable_length(minus_outside);
    for (i = 0; i < n; i++) {
      outside = (Chrpos_obj_T) Uinttable_get(minus_outside,keys[i]);
      Chrpos_obj_free(&outside);
    }
    FREE(keys);
    Uinttable_free(&minus_outside);

    keys = Uinttable_keys(plus_outside,/*sortp*/false);
    n = Uinttable_length(plus_outside);
    for (i = 0; i < n; i++) {
      outside = (Chrpos_obj_T) Uinttable_get(plus_outside,keys[i]);
      Chrpos_obj_free(&outside);
    }
    FREE(keys);
    Uinttable_free(&plus_outside);

    table_stringlist_gc(&minus_acceptor_genes);
    table_stringlist_gc(&minus_donor_genes);
    table_stringlist_gc(&plus_acceptor_genes);
    table_stringlist_gc(&plus_donor_genes);

    table_stringlist_gc(&minus_acceptor_accs);
    table_stringlist_gc(&minus_donor_accs);
    table_stringlist_gc(&plus_acceptor_accs);
    table_stringlist_gc(&plus_donor_accs);

  }

  IIT_free(&knowngenes_iit);

  if (ignore_table != NULL) {
    n = Table_length(ignore_table);
    accessions = (char **) Table_keys(ignore_table,NULL);
    for (i = 0; i < n; i++) {
      FREE(accessions[i]);
    }
    FREE(accessions);
    Table_free(&ignore_table);
  }


  return 0;
}


