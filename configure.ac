AC_INIT([l2_subset], [1.0])
AM_INIT_AUTOMAKE
AC_PROG_CC

# Check for math library
AC_CHECK_LIB([m], [pow])

# Set CFLAGS
AC_SUBST([AM_CFLAGS], ["-Wall -O3 -fsanitize=address -fsanitize=undefined"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT