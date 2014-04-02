static char rcsid[] = "$Id: doublelist.c 106198 2013-08-28 23:07:34Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "doublelist.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

#define T Doublelist_T
struct T {
  double first;
  T rest;
};


T
Doublelist_push (T list, double elt) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = elt;
  new->rest = list;
  return new;
}

T
Doublelist_pop (T list, double *elt) {
  T head;

  if (list) {
    head = list->rest;
    *elt = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}
  
double
Doublelist_head (T list) {
  return list->first;
}

T
Doublelist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Doublelist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE(prev);
  }
}

T
Doublelist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Doublelist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

double *
Doublelist_to_array (int *n, T list) {
  double *array;
  int i;

  *n = Doublelist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (double *) CALLOC(*n,sizeof(double));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

T
Doublelist_from_string (char *string) {
  T doublelist = NULL;
  char *p = string;
  double x;

  while (sscanf(p,"%lf",&x) > 0) {
    doublelist = Doublelist_push(doublelist,x);
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == ',') {
      p++;
    }
  }
  return doublelist;
}

double
Doublelist_max (T this) {
  T p;
  double maxvalue = 0.0;

  if (this != NULL) {
    maxvalue = this->first;
  }
  for (p = this; p; p = p->rest) {
    if (p->first > maxvalue) {
      maxvalue = p->first;
    }
  }
  return maxvalue;
}  

double
Doublelist_min (T this) {
  T p;
  double minvalue = 0.0;

  if (this != NULL) {
    minvalue = this->first;
  }
  for (p = this; p; p = p->rest) {
    if (p->first < minvalue) {
      minvalue = p->first;
    }
  }
  return minvalue;
}  

