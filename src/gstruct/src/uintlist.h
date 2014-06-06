/* $Id: uintlist.h 137605 2014-05-30 00:33:41Z twu $ */
#ifndef UINTLIST_INCLUDED
#define UINTLIST_INCLUDED
#include "bool.h"

#define T Uintlist_T
typedef struct T *T;

extern T 
Uintlist_push (T list, unsigned int x);
extern T 
Uintlist_pop (T list, unsigned int *x);
extern unsigned int 
Uintlist_head (T list);
extern T 
Uintlist_next (T list);
extern void 
Uintlist_head_set (T list, unsigned int x);
extern void 
Uintlist_free (T *list);
extern T 
Uintlist_reverse (T list);
extern int 
Uintlist_length (T list);
extern bool
Uintlist_equal (T list1, T list2);
extern unsigned int *
Uintlist_to_array (int *n, T list);
extern T 
Uintlist_copy (T list);
extern T 
Uintlist_append (T list, T tail);
extern unsigned int 
Uintlist_last_value (T this);
extern unsigned int 
Uintlist_index (T this, int index);
extern T
Uintlist_remove (T list, unsigned int value);
extern bool
Uintlist_find (T this, unsigned int value);
extern char *
Uintlist_to_string (T this);
#undef T
#endif
