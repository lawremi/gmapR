SUBDIRS = gmap gstruct

.PHONY: all clean $(SUBDIRS)

OBJECTS = bamreader.o bamtally.o iit.o variantsummary.o genome.o R_init_gmapR.o

R_SRC_DIR = ${CURDIR}
PREFIX = ${R_SRC_DIR}/../inst/usr
INCLUDE_DIR = $(PREFIX)/include

GSTRUCT_LIB = $(PREFIX)/$(LIBnn)/libgstruct-1.0.a
GSTRUCT_INCLUDE_DIR = $(INCLUDE_DIR)/gstruct

SAMTOOLS_LIB = samtools/libbam.a

PKG_CPPFLAGS += -I$(INCLUDE_DIR)
PKG_CFLAGS += -g -O3
PKG_LIBS += $(GSTRUCT_LIB) $(SAMTOOLS_LIB) -lz

SHLIB = gmapR.so

all: $(SHLIB) gmap

$(SHLIB): $(OBJECTS) $(GSTRUCT_LIB)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(ALL_LIBS)

$(OBJECTS): $(GSTRUCT_INCLUDE_DIR) $(OBJECTS:%.o=%.c)

$(GSTRUCT_LIB) $(GSTRUCT_INCLUDE_DIR): gstruct

RSAMTOOLS_PATH := $(shell R_LIBS=$(R_LIBRARY_DIR) $(R_HOME)/bin/Rscript \
		    --vanilla -e 'cat(system.file(package="Rsamtools"))')
${R_SRC_DIR}/samtools: $(RSAMTOOLS_PATH)
## gmap/gstruct assume samtools headers and libs are in one directory,
## so we need to create one and populate it with links to Rsamtools.
	rm -rf samtools; mkdir samtools
	ln -sf $(RSAMTOOLS_PATH)/usrlib/$(R_ARCH)/libbam.a samtools/libbam.a
	for header in $(RSAMTOOLS_PATH)/include/samtools/*.h; do \
	  ln -sf $$header samtools/$(notdir $(header)); \
	done

$(SUBDIRS): %: %/Makefile
	cd $@; \
	$(MAKE) install

gstruct/Makefile: gstruct/configure ${R_SRC_DIR}/samtools 
	cd $(dir $@); \
	./configure --enable-static --disable-shared \
	            --prefix=${PREFIX} --includedir=${GSTRUCT_INCLUDE_DIR} \
	            --libdir=${PREFIX}/${LIBnn} \
	            --with-samtools-lib=${R_SRC_DIR}/samtools \
                    --disable-binaries  \
                    --disable-maintainer-mode \
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

