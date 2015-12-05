# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# $Id: configure.ac,v 1.46 2010-07-21 22:05:35 twu Exp $

AC_PREREQ(2.57)

## We want to get the version number from file 'VERSION' (rather than
## specifying the version info in 'configure.ac'.  Hence, we need a bit
## of M4 magic.  Note that M4 esyscmd has a trailing newline because the
## shell command output has one, hence the patsubst() trick.
#m4_define([PKG_VERSION],
#          m4_bpatsubst(m4_esyscmd([cat VERSION]),
#                       [\([0-9.]*\)\(\w\|\W\)*],
#                       [\1]))
#m4_define([PKG_VERSION],m4_esyscmd([cat VERSION]))
m4_define([PKG_VERSION],esyscmd([cat VERSION]))

AC_INIT([gstruct],PKG_VERSION,[Thomas Wu <twu@gene.com>])
AC_MSG_CHECKING(package version)
AC_MSG_RESULT(PKG_VERSION)

AM_MAINTAINER_MODE([enable])

### Read defaults

# Handle CONFIG_SITE better than Autoconf does
for ac_site_file in $CONFIG_SITE; do
  if { (eval echo $ac_site_file | grep '^\./' >/dev/null 2>&1)
       ac_status=$?
       (exit $ac_status); }; then
    # ac_site_file starts with ./
    if test -r "$ac_site_file"; then
      ac_site_file="$ac_site_file"
      echo "really loading site script $ac_site_file: file was found"
      . "$ac_site_file"
      ax_user_site_file_loaded=yes
    else
      AC_MSG_ERROR([cannot find $ac_site_file])
    fi
  else
    # ac_site_file does not start with ./
    if test -r "./$ac_site_file"; then
      echo "really loading site script ./$ac_site_file: file was found"
      . "./$ac_site_file"
      ax_user_site_file_loaded=yes
    elif test -r "$ac_site_file"; then
      echo "really loading site script $ac_site_file: file was found"
      . "$ac_site_file" 2>/dev/null
      ax_user_site_file_loaded=yes
    else
      AC_MSG_ERROR([cannot find $ac_site_file])
    fi
  fi
done

if test x"$ax_user_site_file_loaded" != x"yes"; then
   if test -r "./config.site"; then
     echo "loading default site script ./config.site"
     . "./config.site"
   fi
fi




# Set default CFLAGS if not already set by user
AC_MSG_CHECKING(CFLAGS)
AC_ARG_VAR([CFLAGS], [Compiler flags (default: -O3)])
if test x"$CFLAGS" = x; then
        AC_MSG_RESULT(not set by user so using default -O3)
        ACX_EXPAND(CFLAGS,'-O3')
else
	AC_MSG_RESULT($CFLAGS)
fi


AC_CONFIG_SRCDIR([src/bamtally.c])
AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_AUX_DIR([config])
# For autotest: AC_CONFIG_TESTDIR([tests],[tests:src:util])

AC_CONFIG_MACRO_DIR([config])

AC_SUBST([LIBGSTRUCT_SO_VERSION], [1:0:0])
AC_SUBST([LIBGSTRUCT_API_VERSION], [1.0])

AC_CANONICAL_SYSTEM
AC_SYS_LARGEFILE
AC_ARG_PROGRAM

#AM_INIT_AUTOMAKE([no-dependencies])
#AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)
AM_INIT_AUTOMAKE

AM_CONDITIONAL(FULLDIST,test "x$enable_fulldist" = xyes)
AC_ARG_ENABLE([fulldist],
	       AC_HELP_STRING([--enable-fulldist],
                              [For use by program maintainer]),
	      [enable_fulldist="yes"],
              [enable_fulldist="no"])

AM_CONDITIONAL(MAINTAINER,test "x$enable_maintainer" = xyes)
AC_ARG_ENABLE([maintainer],
	       AC_HELP_STRING([--enable-maintainer],
                              [For use by program maintainer]),
	      [enable_maintainer="yes"],
              [enable_maintainer="no"])

# For Perl scripts, because otherwise configure inserts 
# @bindir@ => "${exec_prefix}/bin";
AC_MSG_CHECKING(bindir)
ACX_EXPAND(BINDIR,$bindir)
AC_SUBST(BINDIR)
AC_MSG_RESULT($BINDIR)

# Works in conjunction with AC_PROG_LIBTOOL
#LT_INIT([disable-static])
AC_MSG_CHECKING(whether to link statically)
AC_ARG_ENABLE([static-linking],
               AC_HELP_STRING([--enable-static-linking],
                              [Link binaries statically (default=no)]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
	xyes)
	AC_MSG_RESULT(enabled)
	STATIC_LDFLAG="-all-static"
	;;
	
	xno)
	AC_MSG_RESULT(disabled)
	;;

	x)
	AC_MSG_RESULT(not specified so disabled by default)
	;;
esac	
AC_SUBST(STATIC_LDFLAG)


# Checks for programs.
#AC_PATH_PROG([PERL],[perl])
ACX_PATH_PERL

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_LIBTOOL

# Generates X_CFLAGS, X_LIBS, X_DISPLAY_MISSING, X_EXTRA_LIBS, X_PRE_LIBS
AC_PATH_XTRA
AC_PATH_PROG([QSUB],[qsub])


########################################################################
#   Checks for programs from GMAP package
########################################################################

# Genome directory
AC_MSG_CHECKING(gmap program directory, specified by user)
AC_ARG_WITH([gmap],
            AC_HELP_STRING([--with-gmap=DIR],
                           [Default GMAP program directory]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" = x; then
	ACX_EXPAND(GMAP_PATH,$PATH)
else
	GMAP_PATH=$answer
fi
AC_MSG_RESULT($answer)


AC_PATH_PROG([GET_GENOME],[get-genome],false,$GMAP_PATH)
AC_PATH_PROG([IIT_STORE],[iit_store],false,$GMAP_PATH)
AC_PATH_PROG([IIT_GET],[iit_get],false,$GMAP_PATH)
AC_PATH_PROG([GMAP_BUILD],[gmap_build],false,$GMAP_PATH)
AC_PATH_PROG([GMAP],[gmap],false,$GMAP_PATH)
AC_PATH_PROG([GSNAP],[gsnap],false,$GMAP_PATH)
AC_PATH_PROG([SAM_SORT],[sam_sort],false,$GMAP_PATH)
AC_PATH_PROG([UNIQSCAN],[uniqscan],false,$GMAP_PATH)


# Checks for libraries.
AC_CHECK_LIB(m, rint)

#AC_CHECK_LIB(popt, poptGetContext, [answer="yes"], [answer="no"])
#if test x"$answer" == xyes; then
#	AC_DEFINE([HAVE_LIBPOPT],[1],[Define to 1 if you have the 'popt' library (-lpopt)])
#	POPT_LIBS="-lpopt"
#else
#	POPT_LIBS=""
#fi
#AC_SUBST(POPT_LIBS)

# In call to ACX_PTHREAD, don't provide ACTION-IF-FOUND; otherwise,
# HAVE_PTHREAD won't be defined
AC_MSG_CHECKING(for pthreads feature)
AC_ARG_ENABLE([pthreads],
              AC_HELP_STRING([--enable-pthreads],
                             [Enable pthreads (default=yes)]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
	xyes)
        AC_MSG_RESULT(enabled)
	ACX_PTHREAD([],[pthread_warning=yes])
	;;
	
	xno)
	AC_MSG_RESULT(disabled)
	;;

	x)
	AC_MSG_RESULT([not specified so enabled by default])
	ACX_PTHREAD([],[pthread_warning=yes])
	;;
esac	


# Checks for header files.
AC_HEADER_STDC
AC_HEADER_DIRENT
AC_CHECK_HEADERS([fcntl.h limits.h stddef.h stdlib.h string.h strings.h unistd.h sys/types.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_VOLATILE
AC_TYPE_SIZE_T
AC_TYPE_OFF_T
AC_CHECK_TYPES([caddr_t])
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(unsigned long long)
AC_CHECK_SIZEOF(off_t)


# Checks for library functions.
#AC_FUNC_MALLOC -- Causes use of rpl_malloc
#AC_FUNC_REALLOC -- Causes use of rpl_realloc
#AC_FUNC_CLOSEDIR_VOID -- We are checking return value

AC_FUNC_FSEEKO

#AC_FUNC_MMAP        # Checks only private fixed mapping of already-mapped memory

ACX_MMAP_FIXED   # Checks fixed mapping
ACX_MMAP_VARIABLE   # Checks variable mapping

if test x"$acx_mmap_fixed_ok" = xyes; then
  AC_DEFINE(HAVE_MMAP,1,[Define to 1 if you have a working 'mmap' system call.])
else
  if test x"$acx_mmap_variable_ok" = xyes; then
    AC_DEFINE(HAVE_MMAP,1,[Define to 1 if you have a working 'mmap' system call.])
  fi
fi

ACX_MMAP_FLAGS
ACX_MADVISE_FLAGS

AC_CHECK_FUNCS([ceil floor index log madvise memcpy memmove memset munmap pow rint stat64 strtoul sysconf sysctl sigaction])

ACX_STRUCT_STAT64
ACX_PAGESIZE

ACX_FUNC_FOPEN_BINARY
ACX_FUNC_FOPEN_TEXT


AC_MSG_CHECKING(for builtin popcount/clz/ctz features)
AC_ARG_ENABLE([builtin-popcount],
              AC_HELP_STRING([--enable-builtin-popcount],
                             [Enable -mpopcnt if it __builtin_popcount compiles and runs (default=yes).  May want to disable if compiling on one machine and running on another.]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
	xyes)
        AC_MSG_RESULT(enabled)
	ACX_BUILTIN_POPCOUNT
	;;
	
	xno)
	AC_MSG_RESULT(disabled)
	;;

	x)
	AC_MSG_RESULT([not specified so enabled by default])
	ACX_BUILTIN_POPCOUNT
	;;
esac	
if test x"$acx_mpopcnt_ok" = x"yes"; then
  POPCNT_CFLAGS="$POPCNT_CFLAGS -mpopcnt"
fi
AC_SUBST(POPCNT_CFLAGS)

ACX_ASM_BSR


AC_MSG_CHECKING(whether sse2 is enabled)
AC_ARG_ENABLE([sse2],
	      AC_HELP_STRING([--enable-sse2],
                             [Enable sse2 simd commands if they compile and run (default=yes).]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
     xyes)
     AC_MSG_RESULT(enabled)
     ax_cv_want_sse2_ext=yes
     ;;

     xno)
     AC_MSG_RESULT(disabled by user)
     ax_cv_want_sse2_ext=no
     ;;

     x)
     AC_MSG_RESULT([not specified so enabled by default])
     ax_cv_want_sse2_ext=yes
     ;;
esac	

ax_cv_sse2_shift_defect=no
if test "$ax_cv_want_sse2_ext" = yes; then
ACX_SSE2_SHIFT_DEFECT
fi


AC_MSG_CHECKING(whether ssse3 is enabled)
AC_ARG_ENABLE([ssse3],
	      AC_HELP_STRING([--enable-ssse3],
                             [Enable ssse3 simd commands if they compile and run (default=yes).  Requires that sse2 be enabled.]),
              [answer="$enableval"],
              [answer=""])
if test "$ax_cv_want_sse2_ext" = no; then
   AC_MSG_RESULT([disabled because the user disabled sse2])
   ax_cv_want_ssse3_ext=no
else
    case x"$answer" in
	 xyes)
	 AC_MSG_RESULT(enabled)
	 ax_cv_want_ssse3_ext=yes
	 ;;

	 xno)
	 AC_MSG_RESULT(disabled by user)
	 ax_cv_want_ssse3_ext=no
	 ;;

	 x)
	 AC_MSG_RESULT([not specified so enabled by default])
	 ax_cv_want_ssse3_ext=yes
	 ;;
    esac	
fi



AC_MSG_CHECKING(whether sse4.1 is enabled)
AC_ARG_ENABLE([sse4.1],
	      AC_HELP_STRING([--enable-sse4.1],
                             [Enable sse4.1 simd commands if they compile and run (default=yes).  Requires that ssse3 be enabled.]),
              [answer="$enableval"],
              [answer=""])
if test "$ax_cv_want_ssse3_ext" = no; then
   AC_MSG_RESULT([disabled because the user disabled ssse3])
   ax_cv_want_sse41_ext=no
else
    case x"$answer" in
	 xyes)
	 AC_MSG_RESULT(enabled)
	 ax_cv_want_sse41_ext=yes
	 ;;

	 xno)
	 AC_MSG_RESULT(disabled by user)
	 ax_cv_want_sse41_ext=no
	 ;;

	 x)
	 AC_MSG_RESULT([not specified so enabled by default])
	 ax_cv_want_sse41_ext=yes
	 ;;
    esac	
fi


AC_MSG_CHECKING(whether sse4.2 is enabled)
AC_ARG_ENABLE([sse4.2],
	      AC_HELP_STRING([--enable-sse4.2],
                             [Enable sse4.2 simd commands if they compile and run (default=yes).  Requires that sse4.1 be enabled.]),
              [answer="$enableval"],
              [answer=""])
if test "$ax_cv_want_sse41_ext" = no; then
   AC_MSG_RESULT([disabled because the user disabled sse4.1])
   ax_cv_want_sse42_ext=no
else
    case x"$answer" in
	 xyes)
	 AC_MSG_RESULT(enabled)
	 ax_cv_want_sse42_ext=yes
	 ;;

	 xno)
	 AC_MSG_RESULT(disabled by user)
	 ax_cv_want_sse42_ext=no
	 ;;

	 x)
	 AC_MSG_RESULT([not specified so enabled by default])
	 ax_cv_want_sse42_ext=yes
	 ;;
    esac	
fi


AC_MSG_CHECKING(whether avx is enabled)
AC_ARG_ENABLE([avx],
	      AC_HELP_STRING([--enable-avx],
                             [Enable avx simd commands if they compile and run (default=yes).  Requires that sse4.2 be enabled.]),
              [answer="$enableval"],
              [answer=""])
if test "$ax_cv_want_sse42_ext" = no; then
   AC_MSG_RESULT([disabled because the user disabled sse4.2])
   ax_cv_want_avx_ext=no
else
    case x"$answer" in
	 xyes)
	 AC_MSG_RESULT(enabled)
	 ax_cv_want_avx_ext=yes
	 ;;

	 xno)
	 AC_MSG_RESULT(disabled by user)
	 ax_cv_want_avx_ext=no
	 ;;

	 x)
	 AC_MSG_RESULT([not specified so enabled by default])
	 ax_cv_want_avx_ext=yes
	 ;;
    esac	
fi


AC_MSG_CHECKING(whether avx2 is enabled)
AC_ARG_ENABLE([avx2],
	      AC_HELP_STRING([--enable-avx2],
                             [Enable avx2 simd commands if they compile and run (default=yes).  Requires that avx be enabled.]),
              [answer="$enableval"],
              [answer=""])
if test "$ax_cv_want_avx_ext" = no; then
   AC_MSG_RESULT([disabled because the user disabled avx])
   ax_cv_want_avx2_ext=no
else
    case x"$answer" in
	 xyes)
	 AC_MSG_RESULT(enabled)
	 ax_cv_want_avx2_ext=yes
	 ;;

	 xno)
	 AC_MSG_RESULT(disabled by user)
	 ax_cv_want_avx2_ext=no
	 ;;

	 x)
	 AC_MSG_RESULT([not specified so enabled by default])
	 ax_cv_want_avx2_ext=yes
	 ;;
    esac	
fi


AC_MSG_CHECKING(whether simd is enabled)
AC_ARG_ENABLE([simd],
	      AC_HELP_STRING([--enable-simd],
                             [Enable simd commands in general if they compile and run (default=yes).]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
     xyes)
     AC_MSG_RESULT(enabled)
     AX_EXT   # Sets SIMD_CFLAGS, and HAVE_ALTIVEC,HAVE_MMX,HAVE_SSE,HAVE_SSE2,HAVE_SSE3,HAVE_SSSE3,HAVE_SSE4_1,HAVE_SSE4_2,HAVE_AVX
     ;;

     xno)
     AC_MSG_RESULT(disabled)
     ;;

     x)
     AC_MSG_RESULT([not specified so enabled by default])
     AX_EXT   # Sets SIMD_CFLAGS, and HAVE_ALTIVEC,HAVE_MMX,HAVE_SSE,HAVE_SSE2,HAVE_SSE3,HAVE_SSSE3,HAVE_SSE4_1,HAVE_SSE4_2,HAVE_AVX
     ;;
esac	


# Genome directory
AC_MSG_CHECKING(gmapdb)
AC_ARG_WITH([gmapdb],
            AC_HELP_STRING([--with-gmapdb=DIR],
                           [Default GMAP database directory]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" = x; then
	ACX_EXPAND(GMAPDB,'${prefix}/share')
else
	GMAPDB=$answer
fi
AC_SUBST(GMAPDB)
AC_MSG_RESULT($GMAPDB)


# Determine whether to make binaries
AC_MSG_CHECKING(whether to make binaries)
AC_ARG_ENABLE([binaries],
              AC_HELP_STRING([--enable-binaries],
                             [Enable building of binaries (default=yes)]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
   xno)
   AC_MSG_RESULT(disabled)
   binaries=false
   ;;

   *)
   AC_MSG_RESULT(enabled)
   binaries=true
   ;;

esac
AM_CONDITIONAL([MAKE_BINARIES], [test "x$binaries" = xtrue])


# zlib package

AC_MSG_CHECKING(for zlib support)
AC_ARG_ENABLE([zlib],
              AC_HELP_STRING([--enable-zlib],
                             [Enable zlib support (option needed for uncompressing gzip files) (default=yes)]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
   xno)
   AC_MSG_RESULT(disabled)
   ZLIB_LIBS=""
   have_zlib=no
   ;;

   *)
   AC_MSG_RESULT(enabled)
   failed=0;
   passed=0;
   AC_CHECK_HEADER(zlib.h,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_CHECK_LIB(z,gzopen,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_CHECK_LIB(z,gzeof,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_CHECK_LIB(z,gzgetc,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_CHECK_LIB(z,gzgets,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_CHECK_LIB(z,gzclose,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
   AC_MSG_CHECKING(if zlib package is complete)
   if test $failed -gt 0
   then
      AC_MSG_RESULT(no -- some components failed test)
      ZLIB_LIBS=""
      have_zlib=no
   else
      AC_MSG_RESULT(working)
      AC_DEFINE(HAVE_ZLIB,1,[Define to 1 if you have a working zlib library.])
      ZLIB_LIBS="-lz"
      have_zlib=yes

      AC_CHECK_LIB(z,gzbuffer,have_gzbuffer=1,have_gzbuffer=0)
      if test $have_gzbuffer -gt 0
      then
        AC_DEFINE(HAVE_ZLIB_GZBUFFER,1,[Define to 1 if your zlib library has a gzbuffer function.])
      fi
   fi
   ;;

esac
AC_SUBST(ZLIB_LIBS)


# Samtools program

AC_MSG_CHECKING(for samtools program)
AC_ARG_WITH([samtools-program],
            AC_HELP_STRING([--with-samtools-program=FILE],
                           [Full path of samtools program]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" != x; then
  SAMTOOLS=$answer
  AC_MSG_RESULT(provided by --with-samtools-program)
else
  AC_ARG_WITH([samtools],
              AC_HELP_STRING([--with-samtools=DIR],
                             [Location of samtools bin, include, and lib subdirs]),
              [answer="$withval"],
              [answer=""])
  if test x"$answer" != x; then
    if test -r "$answer/bin/samtools"; then
      SAMTOOLS=$answer/bin/samtools
      AC_MSG_RESULT(provided by --with-samtools and found)
    else
      AC_PATH_PROG([SAMTOOLS],[samtools],false,$PATH)
    fi
  fi
fi


# Samtools headers

#AC_MSG_CHECKING(for samtools headers)
#AC_ARG_WITH([samtools-headers],
#            AC_HELP_STRING([--with-samtools-headers=DIR],
#                           [Directory containing samtools *.h header files]),
#            [answer="$withval"],
#            [answer=""])
#if test x"$answer" != x; then
#  SAMTOOLS_CFLAGS="-I$answer"
#  AC_MSG_RESULT(provided by --with-samtools-headers)
#else
#  AC_ARG_WITH([samtools],
#              AC_HELP_STRING([--with-samtools=DIR],
#                             [Location of samtools bin, include, and lib subdirs]),
#              [answer="$withval"],
#              [answer=""])
#  if test x"$answer" != x; then
#    if test -r "$answer/include/bam.h"; then
#      SAMTOOLS_CFLAGS="-I$answer/include"
#      AC_MSG_RESULT(provided by --with-samtools and found)
#    else
#      SAMTOOLS_CFLAGS="-I$answer/include"
#      AC_MSG_RESULT(not found)
#    fi
#  fi
#fi

SAMTOOLS_CFLAGS="-I."
AC_SUBST(SAMTOOLS_CFLAGS)


# Samtools library

AC_MSG_CHECKING(for samtools library)
AC_ARG_WITH([samtools-lib],
            AC_HELP_STRING([--with-samtools-lib=DIR],
                           [Directory containing samtools libbam.a library]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" != x; then
  AC_DEFINE(HAVE_SAMTOOLS_LIB,1,[Define to 1 if you have a working Samtools library.])
  SAMTOOLS_LDFLAGS="-L$answer -Wl,-rpath,$answer"
  SAMTOOLS_LIBS="-lbam"
  AC_MSG_RESULT(provided by --with-samtools-lib)
else
  AC_ARG_WITH([samtools],
              AC_HELP_STRING([--with-samtools=DIR],
                             [Location of samtools bin, include, and lib subdirs]),
              [answer="$withval"],
              [answer=""])
  if test x"$answer" != x; then
    if test -r "$answer/lib/libbam.a"; then
      AC_DEFINE(HAVE_SAMTOOLS_LIB,1,[Define to 1 if you have a working Samtools library.])
      SAMTOOLS_LDFLAGS="-L$answer/lib -Wl,-rpath,$answer/lib"
      SAMTOOLS_LIBS="-lbam"
      AC_MSG_RESULT(provided by --with-samtools and found)
    else
     SAMTOOLS_LDFLAGS=""
     SAMTOOLS_LIBS=""
     AC_MSG_RESULT(provided by --with-samtools but $answer/lib/libbam.a not found)
    fi
  fi
fi
AC_SUBST(SAMTOOLS_LDFLAGS)
AC_SUBST(SAMTOOLS_LIBS)
AM_CONDITIONAL([HAVE_SAMTOOLS_LIB], [test "x$SAMTOOLS_LIBS" != x])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([lib/gstruct-${LIBGSTRUCT_API_VERSION}.pc:lib/gstruct.pc.in])

AC_OUTPUT

if test x"$perl_warning" = xyes; then
  AC_MSG_WARN([

WARNING: You don't seem to have perl installed with working
versions of the IO::File and Getopt::Std modules.  If you know
where such a version of Perl exists, please specify its full
path under PERL in the config.site file and run configure again.
])
fi

if test x"$pthread_warning" = xyes; then
  AC_MSG_WARN([

WARNING: You don't seem to have headers and libraries for pthreads.
GMAP will work without these, but the program can run faster on
multiprocessor machines if POSIX threads are available.  If you wish
to enable pthreads, please specify the appropriate values for CC,
PTHREAD_CFLAGS, and PTHREAD_LIBS in the config.site file and run
configure again.  Otherwise, you may proceed with the rest of the
installation.
])
fi

if test x"$GSNAP" = xfalse; then
  AC_MSG_WARN([

WARNING: You don't seem to have gsnap specified.
])
fi

if test x"$SAMTOOLS" = xfalse; then
  AC_MSG_WARN([

WARNING: You don't seem to have samtools specified.
])
fi

AC_MSG_CHECKING(Standard compiler and flags to be used)
AC_MSG_RESULT($CC $CFLAGS)
AC_MSG_CHECKING(Standard linker flags to be used)
AC_MSG_RESULT($LD_FLAGS)

AC_MSG_CHECKING(pthread compiler flags to be used)
AC_MSG_RESULT($PTHREAD_CFLAGS)
AC_MSG_CHECKING(popcnt compiler flags to be used)
AC_MSG_RESULT($POPCNT_CFLAGS)

AC_MSG_CHECKING(SIMD features available on computer)
AC_MSG_RESULT($ax_cv_cpu_features)
AC_MSG_CHECKING(SIMD compiler flags to be used)
AC_MSG_RESULT($SIMD_CFLAGS)

