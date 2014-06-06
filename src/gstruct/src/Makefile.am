
if MAKE_BINARIES
bin_PROGRAMS = 

if SAMTOOLS
bin_PROGRAMS += bam_tally
endif

endif


lib_LTLIBRARIES = libgstruct-@LIBGSTRUCT_API_VERSION@.la



BAM_TALLY_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 uinttable.c uinttable.h table.c table.h \
 chrnum.c chrnum.h \
 bamread.c bamread.h samread.c samread.h parserange.c parserange.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h \
 matchdef.h matchpool.c matchpool.h mismatchdef.h mismatchpool.c mismatchpool.h tally.c tally.h bamtally.c bamtally.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h bamtally_main.c

if SAMTOOLS
bam_tally_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SIMD_FLAGS) $(POPCNT_FLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
bam_tally_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_tally_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_tally_SOURCES = $(BAM_TALLY_FILES)
endif


LIBGSTRUCT_LA_FILES = config.h fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 complement.h chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h sequence.c sequence.h genome.c genome.h \
 uinttable.c uinttable.h table.c table.h \
 chrnum.c chrnum.h \
 bamread.c bamread.h samread.c samread.h parserange.c parserange.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h \
 matchdef.h matchpool.c matchpool.h mismatchdef.h mismatchpool.c mismatchpool.h tally.c tally.h bamtally.c bamtally.h \
 datadir.c datadir.h

# Do not add SAMTOOLS LIBS here
if SAMTOOLS
libgstruct_@LIBGSTRUCT_API_VERSION@_la_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SIMD_FLAGS) $(POPCNT_FLAGS) $(SAMTOOLS_CFLAGS) -fPIC -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
libgstruct_@LIBGSTRUCT_API_VERSION@_la_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) -version-info $(LIBGSTRUCT_SO_VERSION)
libgstruct_@LIBGSTRUCT_API_VERSION@_la_LIBADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
dist_libgstruct_@LIBGSTRUCT_API_VERSION@_la_SOURCES = $(LIBGSTRUCT_LA_FILES)
endif


#gstruct_includedir = $(includedir)/gstruct-$(LIBGSTRUCT_API_VERSION)
#gstruct_include_HEADERS = bamtally.h

gstruct_includedir = $(includedir)
gstruct_include_HEADERS = config.h fopen.h bool.h types.h \
 except.h assert.h mem.h \
 intlistdef.h intlist.h listdef.h list.h \
 littleendian.h bigendian.h \
 interval.h uintlist.h \
 stopwatch.h access.h iitdef.h iit-read.h \
 chrom.h genomicpos.h \
 md5.h complement.h sequence.h genome.h \
 uinttable.h table.h \
 chrnum.h \
 bamread.h samread.h parserange.h \
 dynprog.h dynprog_simd.h dynprog_cigar.h dynprog_single.h \
 tableuint.h iit-write.h ucharlist.h \
 matchdef.h matchpool.h mismatchdef.h mismatchpool.h bamtally.h datadir.h

# May want nodist_bamtally_libinclude_HEADERS = ../lib/bamtally.h
#bamtally_libincludedir = $(libdir)/bamtally-$(LIBGSTRUCT_API_VERSION)/include
#nodist_bamtally_libinclude_HEADERS = bamtally-config.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ../lib/gstruct-$(LIBGSTRUCT_API_VERSION).pc



