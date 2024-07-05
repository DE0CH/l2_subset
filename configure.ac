AC_INIT([l2_subset],[1.0])
AM_INIT_AUTOMAKE
AC_PROG_CC
# Save the current LIBS value
save_LIBS=$LIBS

# Check for GSL and CBLAS
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [
  AC_MSG_ERROR([CBLAS library not found. Please install GSL or specify its location in CPPFLAGS and LDFLAGS])
])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm], [], [
  AC_MSG_ERROR([GSL library not found. Please install GSL or specify its location in CPPFLAGS and LDFLAGS])
], [-lgslcblas])

# Check for GSL headers
AC_CHECK_HEADERS([gsl/gsl_blas.h], [], [
  AC_MSG_ERROR([GSL headers not found. Please install GSL or specify its location in CPPFLAGS and LDFLAGS])
])

# Restore the original LIBS value
LIBS=$save_LIBS

# Add GSL libraries to LIBS for the final link
LIBS="$LIBS -lgsl -lgslcblas"
AC_CHECK_LIB([m], [pow], [], [AC_MSG_ERROR([libm (math library) is required but not found])])
AC_CHECK_LIB([gsl], [gsl_qrng_alloc], [], [AC_MSG_ERROR([libgsl is required but not found])])
AC_CHECK_HEADER([gsl/gsl_qrng.h], [], [AC_MSG_ERROR([gsl/gls_qrng.h is required but not found])])
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR([math.h is required but not found])])
AC_CHECK_HEADER([stdio.h], [], [AC_MSG_ERROR([stdio.h is required but not found])])
AC_CHECK_HEADER([stdlib.h], [], [AC_MSG_ERROR([stdlib.h is required but not found])])
AC_CHECK_HEADER([string.h], [], [AC_MSG_ERROR([string.h is required but not found])])

AC_CONFIG_FILES([Makefile])
CFLAGS="-Wall -Werror -O3"
AC_OUTPUT
