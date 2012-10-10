/* $Id: listdef.h 46991 2011-09-12 17:36:30Z twu $ */
#ifndef LISTDEF_INCLUDED
#define LISTDEF_INCLUDED

#define T List_T
struct T {
  void *first;
  struct T *rest;
};

#undef T
#endif

