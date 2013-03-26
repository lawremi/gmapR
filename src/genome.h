#ifndef GENOME_H
#define GENOME_H

#include <gstruct/iit-read.h>
#include <gstruct/genome.h>

/* Internal genome helper functions */

IIT_T readChromosomeIIT (const char *genome_dir, const char *db);
Genome_T createGenome (const char *genome_dir, const char *db);

#endif
