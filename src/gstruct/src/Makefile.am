
if MAKE_BINARIES
bin_PROGRAMS = 

if SAMTOOLS
bin_PROGRAMS += bam_tally bam_get bam_fastq
endif

endif


lib_LTLIBRARIES = libgstruct-@LIBGSTRUCT_API_VERSION@.la



EXON_EXON_JUNCTIONS_FILES = bool.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 fopen.h types.h littleendian.c littleendian.h bigendian.c bigendian.h \
 stopwatch.c stopwatch.h access.c access.h \
 listdef.h list.c list.h \
 uintlist.c uintlist.h \
 interval.c interval.h iitdef.h iit-read.c iit-read.h \
 genomicpos.c genomicpos.h \
 intlistdef.h intlist.c intlist.h md5.c md5.h \
 chrnum.c chrnum.h \
 complement.h sequence.c sequence.h genome.c genome.h \
 table.c table.h \
 datadir.c datadir.h \
 exon-exon-junctions.c

exon_exon_junctions_CC = $(PTHREAD_CC)
exon_exon_junctions_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" 
exon_exon_junctions_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
exon_exon_junctions_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
dist_exon_exon_junctions_SOURCES = $(EXON_EXON_JUNCTIONS_FILES)



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
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h dynprog_single.c dynprog_single.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
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
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h dynprog_single.c dynprog_single.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
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
 dynprog.h dynprog_simd.h dynprog_single.h \
 tableuint.h iit-write.h ucharlist.h bamtally.h datadir.h

# May want nodist_bamtally_libinclude_HEADERS = ../lib/bamtally.h
#bamtally_libincludedir = $(libdir)/bamtally-$(LIBGSTRUCT_API_VERSION)/include
#nodist_bamtally_libinclude_HEADERS = bamtally-config.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ../lib/gstruct-$(LIBGSTRUCT_API_VERSION).pc



BAM_PILEUP_FILES = fopen.h bool.h types.h \
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
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h bam_pileup.c

if SAMTOOLS
bam_pileup_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
bam_pileup_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_pileup_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_pileup_SOURCES = $(BAM_PILEUP_FILES)
endif


BAM_GET_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 interval.c interval.h iitdef.h iit-read.c iit-read.h \
 chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 uinttable.c uinttable.h table.c table.h \
 chrnum.c chrnum.h \
 bamread.c bamread.h samread.c samread.h \
 parserange.c parserange.h \
 getopt.c getopt1.c getopt.h bam_get.c

bam_get_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DBAM_INPUT=1
bam_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_get_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_get_SOURCES = $(BAM_GET_FILES)


BAM_PRINT_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 uinttable.c uinttable.h table.c table.h \
 chrnum.c chrnum.h \
 bamread.c bamread.h samread.c samread.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
 orderstat.c orderstat.h lgamma.c lgamma.h tally.c tally.h \
 tableint.c tableint.h splice.c splice.h maxent_hr.c maxent_hr.h gstruct.c gstruct.h \
 datadir.c datadir.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h bam_print.c

bam_print_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBAM_INPUT=1
bam_print_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_print_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_print_SOURCES = $(BAM_PRINT_FILES)


BAM_FASTQ_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 listdef.h list.c list.h \
 genomicpos.c genomicpos.h complement.h \
 intlistdef.h intlist.c intlist.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 md5.c md5.h sequence.c sequence.h \
 genome.c genome.h chrom.c chrom.h \
 uinttable.c uinttable.h \
 interval.c interval.h iit-read.c iit-read.h \
 bamread.c bamread.h samread.c samread.h \
 table.c table.h \
 getopt.c getopt1.c getopt.h bam_fastq.c

bam_fastq_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
bam_fastq_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_fastq_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_fastq_SOURCES = $(BAM_FASTQ_FILES)



SPLICECLEAN_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 stopwatch.c stopwatch.h access.c access.h \
 listdef.h list.c list.h \
 genomicpos.c genomicpos.h \
 chrom.c chrom.h \
 uinttable.c uinttable.h table.c table.h \
 orderstat.c orderstat.h lgamma.c lgamma.h tally.c tally.h \
 bigendian.c bigendian.h interval.c interval.h iitdef.h iit-read.c iit-read.h \
 getopt.c getopt1.c getopt.h spliceclean.c

spliceclean_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS)
spliceclean_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
spliceclean_LDADD = $(PTHREAD_LIBS)
dist_spliceclean_SOURCES = $(SPLICECLEAN_FILES)



#SPLICEGENE_FILES = fopen.h bool.h types.h \
# except.c except.h assert.c assert.h mem.c mem.h \
# intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
# littleendian.c littleendian.h bigendian.c bigendian.h \
# interval.c interval.h uintlist.c uintlist.h \
# stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
# uinttable.c uinttable.h genomicpos.c genomicpos.h \
# getopt.c getopt1.c getopt.h splicegene.c
#
#splicegene_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
#splicegene_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
#splicegene_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
#dist_splicegene_SOURCES = $(SPLICEGENE_FILES)


#SAM_GSTRUCT_FILES = fopen.h bool.h types.h \
# except.c except.h assert.c assert.h mem.c mem.h \
# intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
# littleendian.c littleendian.h bigendian.c bigendian.h \
# interval.c interval.h uintlist.c uintlist.h \
# stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
# complement.h chrom.c chrom.h genomicpos.c genomicpos.h \
# md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
# table.c table.h uinttable.c uinttable.h \
# samflags.h samread.c samread.h \
# orderstat.c orderstat.h lgamma.c lgamma.h tally.c tally.h \
# iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
# splice.c splice.h gstruct.c gstruct.h cappaths.c cappaths.h splicegraph.c splicegraph.h \
# datadir.c datadir.h parserange.c parserange.h \
# getopt.c getopt1.c getopt.h gstruct_main.c
#
#sam_gstruct_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
#sam_gstruct_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
#sam_gstruct_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
#dist_sam_gstruct_SOURCES = $(SAM_GSTRUCT_FILES)


BAM_SINGLEMAP_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 table.c table.h uinttable.c uinttable.h \
 chrnum.c chrnum.h \
 samflags.h samread.c samread.h bamread.c bamread.h \
 orderstat.c orderstat.h lgamma.c lgamma.h tally.c tally.h \
 tableint.c tableint.h tableuint.c tableuint.h \
 iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h singlemap.c
bam_singlemap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBAM_INPUT=1
bam_singlemap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_singlemap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_singlemap_SOURCES = $(BAM_SINGLEMAP_FILES)


BAM_GSTRUCT_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 table.c table.h uinttable.c uinttable.h \
 chrnum.c chrnum.h \
 samflags.h samread.c samread.h bamread.c bamread.h \
 orderstat.c orderstat.h lgamma.c lgamma.h tally.c tally.h \
 tableuint.c tableuint.h iit-write.c iit-write.h ucharlist.c ucharlist.h bamtally.c bamtally.h \
 tableint.c tableint.h splice.c splice.h maxent_hr.c maxent_hr.h gstruct.c gstruct.h \
 cappaths.c cappaths.h splicegraph.c splicegraph.h \
 datadir.c datadir.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h gstruct_main.c

bam_gstruct_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(SAMTOOLS_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBAM_INPUT=1
bam_gstruct_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(SAMTOOLS_LDFLAGS)
bam_gstruct_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(SAMTOOLS_LIBS)
dist_bam_gstruct_SOURCES = $(BAM_GSTRUCT_FILES)


