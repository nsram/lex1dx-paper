#!/bin/sh
# Copies the C++ core from the repository (../cpp) into src/, so the package is self-contained.
# configure and configure.win run this before every install from a checkout. Run it yourself, from
# the r/ directory, before `R CMD build`. In a source tarball there is no ../cpp and it does nothing.
set -e
core=../cpp
[ -f "$core/include/lex/lex.hpp" ] || exit 0
cp "$core"/src/*.cpp "$core"/src/*.hpp src/
mkdir -p src/lex
cp "$core/include/lex/lex.hpp" src/lex/
