#!/bin/bash
# HM compilation script.

set -e

# Path to HM repository.
HM_PATH=$HOME/repos/hm

cd $HM_PATH/build/linux
make

