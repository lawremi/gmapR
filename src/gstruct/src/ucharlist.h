/* $Id: ucharlist.h 66857 2012-06-19 21:05:57Z twu $ */
#ifndef UCHARLIST_INCLUDED
#define UCHARLIST_INCLUDED
#include "bool.h"

#define T Ucharlist_T
typedef struct T *T;

extern T 
Ucharlist_push (T list, unsigned char x);
extern T 
Ucharlist_pop (T list, unsigned char *x);
extern unsigned char 
Ucharlist_head (T list);
extern T 
Ucharlist_next (T list);
extern void 
Ucharlist_head_set (T list, unsigned char x);
extern void 
Ucharlist_free (T *list);
extern T 
Ucharlist_reverse (T list);
extern int 
Ucharlist_length (T list);
extern unsigned char *
Ucharlist_to_array (int *n, T list);
extern T 
Ucharlist_copy (T list);
extern T 
Ucharlist_append (T list, T tail);
extern unsigned char 
Ucharlist_last_value (T this);
extern unsigned char 
Ucharlist_index (T this, int index);
extern T
Ucharlist_remove (T list, unsigned char value);
extern bool
Ucharlist_find (T this, unsigned char value);
#undef T
#endif
