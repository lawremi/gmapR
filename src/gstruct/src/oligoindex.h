#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED

#define T Oligoindex_T
typedef struct T *T;

#define STRAIGHT_MASK_8  0x0000FFFF /* 8-mer: 1111 1111 1111 1111 */

extern int
Oligoindex_compare_seqs (int *maxposi, int *maxposj, char *donor_sequence, char *acceptor_sequence,
			 int width);

#undef T
#endif

