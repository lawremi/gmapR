static char rcsid[] = "$Id: ucharlist.c 66857 2012-06-19 21:05:57Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ucharlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T Ucharlist_T
struct T {
  unsigned char first;
  T rest;
};

T
Ucharlist_push (T list, unsigned char x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
Ucharlist_pop (T list, unsigned char *x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}
  
unsigned char
Ucharlist_head (T list) {
  return list->first;
}

T
Ucharlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Ucharlist_head_set (T this, unsigned char x) {
  this->first = x;
  return;
}

void
Ucharlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE(prev);
  }
}

T
Ucharlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Ucharlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

unsigned char *
Ucharlist_to_array (int *n, T list) {
  unsigned char *array;
  int i;

  *n = Ucharlist_length(list);
  array = (unsigned char *) CALLOC(*n,sizeof(unsigned char));
  for (i = 0; i < *n; i++) {
    array[i] = list->first;
    list = list->rest;
  }
  return array;
}

T
Ucharlist_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

T
Ucharlist_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

unsigned char
Ucharlist_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

unsigned char
Ucharlist_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


T
Ucharlist_remove (T list, unsigned char value) {
  T *p, x;

  if (list == NULL) {
    return (T) NULL;
  } else {
    p = &list;
    while (*p) {
      if ((*p)->first == value) {
	x = *p;
	*p = (*p)->rest;
	FREE(x);
      } else {
	p = &(*p)->rest;
      }
    }
    return list;
  }
}


bool
Ucharlist_find (T this, unsigned char value) {
  T r;

  for (r = this; r != NULL; r = r->rest) {
    if (r->first == value) {
      return true;
    }
  }
  return false;
}

