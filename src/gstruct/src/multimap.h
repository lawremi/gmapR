#ifndef MULTIMAP_INCLUDED
#define MULTIMAP_INCLUDED
#include "bool.h"
#include "table.h"
#include "bamread.h"
#include "genome.h"
#include "iit-read.h"

extern Table_T
Multimap_resolve (Bamreader_T bamreader, IIT_T tally_iit, IIT_T chromosome_iit,
		  int minimum_mapq, bool need_unique_p, bool need_primary_p);

#endif

