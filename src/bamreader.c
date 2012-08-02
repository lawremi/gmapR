#include <sys/types.h>
#include <gstruct/bamread.h>

#include "gmapR.h"

void R_Bamread_free(SEXP bamreader_R) {
  Bamreader_T bamreader = R_ExternalPtrAddr(bamreader_R);
  Bamread_free(&bamreader);
}

SEXP
R_Bamread_new (SEXP bamfile_R) {
  SEXP bamreader_R;
  Bamreader_T bamreader;
  const char *bamfile;

  bamfile = CHAR(asChar(bamfile_R));
  
  bamreader = Bamread_new((char *)bamfile);
  if (bamreader == NULL) {
    error("Could not open bamfile\n");
  }

  bamreader_R = R_MakeExternalPtr((void *) bamreader, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(bamreader_R, R_Bamread_free);

  return bamreader_R;
}
