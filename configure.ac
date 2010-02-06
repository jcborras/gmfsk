# Process this file with autoconf to produce a configure script.
# 
AC_PREREQ(2.57)

AC_INIT([gmfsk], [0.7pre1], [oh2bns@sral.fi])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

AH_TEMPLATE([GETTEXT_PACKAGE], [Gettext domain name])

GETTEXT_PACKAGE=gmfsk
AC_SUBST(GETTEXT_PACKAGE)
AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE, "$GETTEXT_PACKAGE")

# Supported languages.
#
ALL_LINGUAS="fi sk"
AM_GLIB_GNU_GETTEXT

# Checks for programs.
#
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_CHECK_PROGS(RPMBUILD, rpmbuild, rpm)

# GConf (some systems hide gconftool in strange places...)
#
AC_PATH_PROG([GCONFTOOL], [gconftool-2], [], [$PATH:/opt/gnome/bin:/opt/gnome2/bin])
AM_GCONF_SOURCE_2

# Scrollkeeper.
#
AC_PATH_PROG(SK_CONFIG, scrollkeeper-config, no)
if test x$SK_CONFIG = xno; then
	AC_MSG_ERROR([Couldn't find scrollkeeper-config.])
fi
SCROLLKEEPER_BUILD_REQUIRED=0.3.5
AC_SUBST(SCROLLKEEPER_BUILD_REQUIRED)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h sys/ioctl.h unistd.h dfftw.h fftw.h linux/ppdev.h])
 
# Checks for typedefs, structures, and compiler characteristics.
#
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions.
#
AC_PROG_GCC_TRADITIONAL
AC_FUNC_STRFTIME
AC_FUNC_VPRINTF
AC_CHECK_LIB([m],[cos])
AC_CHECK_FUNCS([floor ceil fmod memcpy memmove memset pow sqrt strcasecmp strchr strdup strerror strstr mkdir log2 carg csqrt lrint lrintf])

PKG_CHECK_MODULES(gnomeui_2, [libgnomeui-2.0])
PKG_CHECK_MODULES(pulse_simple, [libpulse-simple],[AC_DEFINE([HAVE_PULSE_SIMPLE],[1],[Use PulseAudio])],[AC_MSG_NOTICE([libpulse-simple not found])])

PACKAGE_CFLAGS="$gnomeui_2_CFLAGS $pulse_simple_CFLAGS"
PACKAGE_LIBS="$gnomeui_2_LIBS $pulse_simple_LIBS"

AC_SUBST(PACKAGE_CFLAGS)
AC_SUBST(PACKAGE_LIBS)

AC_CHECK_LIB([hamlib], [rig_init],,[AC_MSG_NOTICE([hamlib not found])])

# Check for this library last so that it doesn't kill checks
# for the libc functions
AC_SEARCH_LIBS([fftw_create_plan], [dfftw fftw],,
	       [AC_MSG_ERROR([FFTW libraries not found!!!])], [-lm])

AC_ARG_ENABLE(
  [debug],
  AS_HELP_STRING([--enable-debug],[Set debug compiler options]),
  # Flags for when we are debugging
  CFLAGS="$CFLAGS -Wall -Werror -g3 -O0",
  # Flags for when we are not debugging
  CFLAGS="$CFLAGS -DNDEBUG"
)

# Output
#
AC_CONFIG_FILES([Makefile
		 gmfsk.spec
		 help/Makefile
		 help/gmfsk/Makefile
		 help/gmfsk/C/Makefile
		 m4/Makefile
		 po/Makefile.in
		 src/Makefile
		 src/cw/Makefile
		 src/feld/Makefile
		 src/mfsk/Makefile
		 src/misc/Makefile
		 src/mt63/Makefile
		 src/olivia/Makefile
		 src/psk31/Makefile
		 src/rtty/Makefile
		 src/samplerate/Makefile
		 src/throb/Makefile])

AC_OUTPUT
