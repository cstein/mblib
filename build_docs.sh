#!/usr/bin/env bash

# get git-revision
sed "s/GITREVISION/`git describe --always --dirty`/" Doxyfile.in > Doxyfile

doxygen Doxyfile

rm -f Doxyfile
