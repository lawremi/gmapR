/* $Id: samread.h 137858 2014-06-02 21:58:56Z twu $ */
#ifndef SAMREAD_INCLUDED
#define SAMREAD_INCLUDED
#include <stdio.h>
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"

extern char *
Samread_get_acc (unsigned int *flag, char *line);

/* Returns pointer to auxinfo */
extern char *
Samread_parse_line (char **acc, unsigned int *flag, int *mapq, char **chr, Genomicpos_T *chrpos, char **cigar,
		    char **mate_chr, Genomicpos_T *mate_chrpos, int *readlength, char **read, char **quality_string,
		    char *line);

extern char *
Samread_chrinfo (Genomicpos_T *chrpos, char **cigar, char *line);

extern void
Samread_print_altered_single (FILE *fp, unsigned int newflag, int mapq, char *line);

extern void
Samread_print_altered_paired (FILE *fp, unsigned int newflag, int mapq, char *mate_chr, Genomicpos_T mate_chrpos,
			      int insert_length, char *line);
extern void
Samread_print_altered_mate (FILE *fp, char *chr, Genomicpos_T chrpos, char *mate_chr, Genomicpos_T mate_chrpos,
			    int insert_length, char *line);
extern void
Samread_print_altered (FILE *fp, unsigned int newflag, int mapq, int insert_length, char *line);

extern char
Samread_splice_strand (char *auxinfo);

extern Intlist_T
Samread_parse_cigar (Uintlist_T *npositions, int *readlength, char *cigar);

extern void
Samread_print_cigar (Intlist_T types, Uintlist_T npositions);

extern bool
Samread_cigar_hardclipped_p (char *cigar);

extern Genomicpos_T
Samread_chrpos_high (Intlist_T types, Uintlist_T npositions, Genomicpos_T chrpos_low);

extern int
Samread_get_query_coordinates (int *query5, int *query3, Intlist_T types, Uintlist_T npositions,
			       int readlength, char *cigar);

extern int
Samread_compute_insert_length (int *querylength5, int *querylength3,
			       char *cigar5, Genomicpos_T chrpos_low_5, char *cigar3, Genomicpos_T chrpos_low_3);

#endif

