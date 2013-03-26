#include <stdlib.h>
#include <string.h>

#include <gstruct/bool.h>
#include <gstruct/datadir.h>
#include <gstruct/iit-read.h>
#include <gstruct/genome.h>
#include <gstruct/interval.h>
#include <gstruct/genomicpos.h>
#include <gstruct/complement.h>

#include "gmapR.h"
#include "genome.h"

IIT_T
readChromosomeIIT (const char *genome_dir, const char *db) {
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
  free(iitfile);
  free(fileroot);
  free(dbversion);
  free(genomesubdir);

  return chromosome_iit;
}

Genome_T
createGenome (const char *genome_dir, const char *db) {
  char *genomesubdir, *fileroot, *dbversion;
  genomesubdir = Datadir_find_genomesubdir(&fileroot, &dbversion,
                                           (char *)genome_dir, (char *)db);
  
  Genome_T genome = Genome_new(genomesubdir, fileroot, /*snps_root*/NULL,
                               /*uncompressedp*/false, /*access*/USE_MMAP_ONLY);
  free(fileroot);
  free(dbversion);
  free(genomesubdir);

  return genome;
}

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_inplace (char *sequence, Genomicpos_T length) {
  char temp;
  int i, j;

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

SEXP
R_Genome_getSeq (SEXP genome_dir_R, SEXP db_R,
                 SEXP seqnames_R, SEXP start_R, SEXP width_R, SEXP strand_R)
{
  const char *genome_dir =
    genome_dir_R == R_NilValue ? NULL : CHAR(asChar(genome_dir_R));
  const char *db = CHAR(asChar(db_R));
  Genome_T genome = createGenome(genome_dir, db);
  IIT_T chromosome_iit = readChromosomeIIT(genome_dir, db);
  int *start = INTEGER(start_R);
  int *width = INTEGER(width_R);

  int max_width = 0;
  for (int i = 0; i < length(width_R); i++) {
    if (width[i] > max_width)
      max_width = width[i];
  }
  char *buffer = (char *) R_alloc(sizeof(char), max_width + 1);

  int seqname_index;
  const char *old_seqname = NULL; // premature optimization
  SEXP result;
  PROTECT(result = allocVector(STRSXP, length(start_R)));
  for (int i = 0; i < length(start_R); i++) {
    const char *seqname = CHAR(STRING_ELT(seqnames_R, i));
    if (old_seqname == NULL || strcmp(seqname, old_seqname)) {
      seqname_index = IIT_find_linear(chromosome_iit, (char *)seqname);
      if (seqname_index < 0) {
        error("Cannot find chromosome %s in genome", seqname);
      }
    }
    old_seqname = seqname;
    int chroffset = Interval_low(IIT_interval(chromosome_iit, seqname_index));
    Genome_fill_buffer_simple(genome, chroffset + start[i] - 1, width[i],
                              buffer);
    if (CHAR(STRING_ELT(strand_R, i))[0] == '-')
      make_complement_inplace(buffer, width[i]);
    SET_STRING_ELT(result, i, mkChar(buffer));
  }

  IIT_free(&chromosome_iit);
  Genome_free(&genome);

  UNPROTECT(1);
  return result;
}

