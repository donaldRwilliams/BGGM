## -*- mode: autoconf; autoconf-indentation: 4; -*-
AC_PREREQ(2.61)

AC_INIT([BGGM], m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Use R to configure compilers
CXX=`${R_HOME}/bin/R CMD config CXX`
CXXFLAGS=`${R_HOME}/bin/R CMD config CXXFLAGS`


## Specify C++ as the programming language
AC_LANG(C++)

## Check for OpenMP support using R macros
#openmp_flag=""
#openmp_cflag=""

## Use R's OpenMP flags: Currently bggm_fast.cpp does not use any parallelization. openMP not necessary
#openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
#openmp_cflag='$(SHLIB_OPENMP_CFLAGS)'

## Substitute OpenMP flags
#AC_SUBST([OPENMP_CFLAG], ["${openmp_cflag}"])
#AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])

## Configure output files
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

