AC_INIT([Netboost], 0.0.2)
AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_MACRO_DIR([m4])

# M4 macros for checking of CPU features.
m4_include([tools/ax_gcc_x86_avx_xgetbv.m4])
m4_include([tools/ax_gcc_x86_cpuid.m4])
m4_include([tools/ax_check_compile_flag.m4])
m4_include([tools/ax_ext.m4])

# Checks for common programs using default macros
AC_PROG_CC
AC_PROG_CPP

# Current CPU 
AX_EXT

CC_ADD_FLAGS=""
AX_CHECK_COMPILE_FLAG(-Wno-deprecated, CC_ADD_FLAGS="$CC_ADD_FLAGS -Wno-deprecated", [])

AC_CHECK_TOOL([STRIP],[strip])

AC_SUBST(CC_ADD_FLAGS)

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars src/mcupgma/clustering_util/Makefile])

# AC_OUTPUT(src/Makevars src/mcupgma/clustering_util/Makefile src/netboost.h)
AC_OUTPUT
