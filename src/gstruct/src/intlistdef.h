/* $Id: intlistdef.h 46991 2011-09-12 17:36:30Z twu $ */
#ifndef INTLISTDEF_INCLUDED
#define INTLISTDEF_INCLUDED

#define T Intlist_T
struct T {
  int first;
  struct T *rest;
};

#undef T
#endif
