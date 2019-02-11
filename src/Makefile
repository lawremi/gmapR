SUBDIRS = gmap gstruct

.PHONY: all clean $(SUBDIRS)

OBJECTS = bamreader.o bamtally.o iit.o variantsummary.o genome.o R_init_gmapR.o

R_SRC_DIR = ${CURDIR}
PREFIX = ${R_SRC_DIR}/../inst/usr
INCLUDE_DIR = $(PREFIX)/include

GSTRUCT_LIB = $(PREFIX)/$(LIBnn)/libgstruct-1.0.a
GSTRUCT_INCLUDE_DIR = $(INCLUDE_DIR)/gstruct

SAMTOOLS_LIB = samtools/libbam.a

PATCH_O = samtools_patch.o
KNETFILE_O = knetfile.o
BAMOBJ_0 = \
  bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
  bam_pileup.o bam_lpileup.o bam_md.o razf.o faidx.o \
  $(KNETFILE_O) bam_sort.o sam_header.o bam_reheader.o kprobaln.o $(PATCH_O)
BAMOBJ=$(BAMOBJ_0:%=samtools/%)
DFLAGS = -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 \
  -U_FORTIFY_SOURCE -DBGZF_CACHE \
  -Dfprintf=_samtools_fprintf \
  -Dexit=_samtools_exit \
  -Dabort=_samtools_abort

PKG_CPPFLAGS += -I$(INCLUDE_DIR)
PKG_CFLAGS += -g -O3 $(DFLAGS)
PKG_LIBS += $(GSTRUCT_LIB) $(SAMTOOLS_LIB) -lz

SHLIB = gmapR.so

all: $(SHLIB) gmap

$(SHLIB): $(OBJECTS) $(GSTRUCT_LIB)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(ALL_LIBS)

$(OBJECTS): $(GSTRUCT_INCLUDE_DIR) $(OBJECTS:%.o=%.c)

$(GSTRUCT_LIB) $(GSTRUCT_INCLUDE_DIR): gstruct

$(SAMTOOLS_LIB): $(BAMOBJ)
	$(AR) -crus $@ $(BAMOBJ)

$(SUBDIRS): %: %/Makefile $(SAMTOOLS_LIB)
	cd $@; \
	$(MAKE) install

gstruct/Makefile: gstruct/configure
	cd $(dir $@); \
        CFLAGS="-g -O3" \
	./configure --enable-static --disable-shared \
	            --prefix=${PREFIX} --includedir=${GSTRUCT_INCLUDE_DIR} \
	            --libdir=${PREFIX}/${LIBnn} \
	            --with-samtools-lib=${R_SRC_DIR}/samtools \
		    --disable-maintainer-mode \
                    --disable-binaries
## does not appear to be a true dependency yet
##	            --with-gmap=${PREFIX}/bin

gmap/Makefile: gmap/configure
	cd $(dir $@); \
	 ./configure --with-gmapdb=${GMAPDB} --prefix=${PREFIX} \
	            --libdir=${PREFIX}/${LIBnn} --disable-maintainer-mode

clean:
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir distclean; \
	done
	rm -rf samtools $(PREFIX)
	rm -f *.o *.so *.dll

