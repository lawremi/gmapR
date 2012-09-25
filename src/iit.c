#include <gstruct/iit-read.h>

#include "iit.h"

static void R_IIT_free(SEXP iit_R) {
  IIT_T iit = R_ExternalPtrAddr(iit_R);
  IIT_free(&iit);
}

SEXP R_IIT_new(IIT_T iit) {
  SEXP iit_R;
  iit_R = R_MakeExternalPtr((void *) iit, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(iit_R, R_IIT_free);
  return iit_R;
}

SEXP R_iit_read(SEXP iitfile_R) {
  char *iitfile = (char *) CHAR(asChar(iitfile_R));
  IIT_T iit = IIT_read(iitfile, /*name*/NULL, /*readonlyp*/true,
                       /*divread*/READ_ALL, /*divstring*/NULL,
                       /*add_iit_p*/false, /*labels_read_p*/true);
  return R_IIT_new(iit);
}

/* 'IIT_output_direct' is not included in the libgstruct binary */

/*
  #include <gstruct/iit-write.h>

  SEXP R_iit_write(SEXP tally_iit_R, SEXP iitfile_R) {
  IIT_T tally_iit = (IIT_T) R_ExternalPtrAddr(tally_iit_R);
  char *iitfile = (char *) CHAR(asChar(iitfile_R));
  IIT_output_direct(iitfile, tally_iit, IIT_LATEST_VERSION);
  return iitfile_R;
  }
*/

