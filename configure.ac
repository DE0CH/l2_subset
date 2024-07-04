AC_INIT([l2_subset], [1.0])
AM_INIT_AUTOMAKE
AC_PROG_CC

# Check for math library
AC_CHECK_LIB([m], [pow], [], [AC_MSG_ERROR([libm (math library) is required but not found])])
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR([math.h is required but not found])])
AC_CHECK_HEADER([stdio.h], [], [AC_MSG_ERROR([stdio.h is required but not found])])
AC_CHECK_HEADER([stdlib.h], [], [AC_MSG_ERROR([stdlib.h is required but not found])])
AC_CHECK_HEADER([string.h], [], [AC_MSG_ERROR([string.h is required but not found])])

# Set CFLAGS
AC_SUBST([AM_CFLAGS], ["-Wall -O3 -fsanitize=address -fsanitize=undefined"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT