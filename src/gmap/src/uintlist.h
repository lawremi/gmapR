/* $Id: uintlist.h 140158 2014-06-28 23:42:12Z twu $ */
#ifndef UINTLIST_INCLUDED
#define UINTLIST_INCLUDED
#include "types.h"
#include "bool.h"

#define T Uintlist_T
typedef struct T *T;

extern T 
Uintlist_push (T list, UINT4 x);
extern T 
Uintlist_pop (T list, UINT4 *x);
extern UINT4 
Uintlist_head (T list);
extern T 
Uintlist_next (T list);
extern void 
Uintlist_head_set (T list, UINT4 x);
extern void 
Uintlist_free (T *list);
extern T 
Uintlist_reverse (T list);
extern int 
Uintlist_length (T list);
extern UINT4 *
Uintlist_to_array (int *n, T list);
extern UINT4 *
Uintlist_to_array_out (int *n, T list);
extern T
Uintlist_from_array (UINT4 *array, int n);
extern T 
Uintlist_copy (T list);
extern T 
Uintlist_append (T list, T tail);
extern UINT4 
Uintlist_last_value (T this);
extern UINT4 
Uintlist_index (T this, int index);
extern bool
Uintlist_find (T this, UINT4 value);
extern char *
Uintlist_to_string (T this);
#undef T
#endif
