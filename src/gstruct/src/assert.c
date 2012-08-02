static char rcsid[] = "$Id: assert.c 46991 2011-09-12 17:36:30Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"

const Except_T Assert_Failed = { "Assertion Failed" };

/*
void
(assert) (int e) {
  assert(e);
}
*/

