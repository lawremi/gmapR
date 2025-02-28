/* $Id: bool.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef BOOL_INCLUDED
#define BOOL_INCLUDED

/* typedef enum{false,true} bool; */

#if __STDC_VERSION__ < 202000
typedef unsigned char bool;
#define false 0
#define true 1
#endif


#endif
