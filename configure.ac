#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(NHDS, 1.5, d.verscharen@ucl.ac.uk)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR(src/)

: ${FCFLAGS=""}

# Checks for programs.
AC_PROG_FC([h5fc f90 fort gfortran])
AC_LANG([Fortran])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_SEARCH_LIBS([h5open_f],[hdf5_fortran])

AC_ARG_WITH([hdf5],[AS_HELP_STRING([−-with−hdf5=DIR],[enable path for HDF5 library])],[WITH_HDF5="yes" HDF5PATH="$withval"], [WITH_HDF5="no"])
AM_CONDITIONAL(WITH_HDF5, test x$WITH_HDF5 = xyes)
if test -n "$WITH_HDF5"; then
  HDF5_INC="$HDF5PATH/include"
  HDF5_LIB="$HDF5PATH/lib"
  HDF5_LDFLAGS="$HDF5_LIB/libhdf5_fortran.a $HDF5_LIB/libhdf5.a -lz"

  AC_SUBST(HDF5_INC)
  LD_LIBRARY_PATH="$HDF5_LIB:$LD_LIBRARY_PATH"
  AC_SUBST(HDF5_LIB)
  AC_SUBST(HDF5_LDFLAGS)
fi




AC_OUTPUT(Makefile src/Makefile)
