#!/bin/bash
# f265 compilation script.

set -e

# Path to f265 repository.
F265_PATH=$HOME/repos/f265

cd $F265_PATH

if [ "$1" = "-c" ]; then
    scons
else
    scons
fi

