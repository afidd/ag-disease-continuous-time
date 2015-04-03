#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([adct], [0.1], [ajd27@cornell.edu])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Checks for libraries.
AX_PTHREAD
AX_BOOST_BASE([1.55])
AX_BOOST_CHRONO
AX_BOOST_FILESYSTEM
AX_BOOST_LOG
AX_BOOST_LOG_SETUP
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
#AC_SEARCH_LIBS([main], [gsl])
AC_CHECK_LIB([hdf5], [main])
AC_CHECK_LIB([hdf5_hl], [main])


# Checks for header files.
AC_ARG_WITH([semimarkov],
  [AS_HELP_STRING([--with-semimarkov],
    [Include directory for Semi-Markov library])],
  [SEMIMARKOV_CXXFLAGS=-I$withval
   AC_SUBST(SEMIMARKOV_CXXFLAGS)],
  [with_semimarkov=check])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([floor])

python getgit.py src/contact_version.hpp

AC_OUTPUT
