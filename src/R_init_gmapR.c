#include <R_ext/Rdynload.h>

#include "gmapR.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

  /* bamtally.c */
  CALLMETHOD_DEF(R_Bamtally_iit, 19),

  /* bamreader.c */
  CALLMETHOD_DEF(R_Bamread_new, 1),
  
  {NULL, NULL, 0}
};

void R_init_gmapR(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}