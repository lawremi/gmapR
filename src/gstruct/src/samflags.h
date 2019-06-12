/* $Id: samflags.h 219289 2019-05-21 01:13:18Z twu $ */
#ifndef SAMFLAGS_INCLUDED
#define SAMFLAGS_INCLUDED

#define PAIRED_READ        0x0001 /* 1 */
#define PAIRED_MAPPING     0x0002 /* 2 */
#define QUERY_UNMAPPED     0x0004 /* 4 */
#define MATE_UNMAPPED      0x0008 /* 8 */
#define QUERY_MINUSP       0x0010 /* 16 */
#define MATE_MINUSP        0x0020 /* 32 */
#define FIRST_READ_P       0x0040 /* 64 */
#define SECOND_READ_P      0x0080 /* 128 */
#define NOT_PRIMARY        0x0100 /* 256 */
#define BAD_READ_QUALITY   0x0200 /* 512 */
#define DUPLICATE_READ     0x0400 /* 1024 */
#define SUPPLEMENTARY      0x0800 /* 2048 */

/* 83 = first read, minus strand for paired */
/* 99 = first read, plus strand for paired */
/* 147 = second read, minus strand for paired */
/* 163 = second read, plus strand for paired */

/* For forcing a read to be primary */
#define SET_PRIMARY        0xFEFF /* do a logical-and (&) with this */

#endif

