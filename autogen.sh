#!/bin/sh
# These two commands build the configure script and Makefile.in
# from configure.ac and Makefile.am. They are run by the maintainer
# and rerun whenever configure.ac changes.
rm -f src/Makefile.in
autoreconf --install
automake --add-missing --copy
