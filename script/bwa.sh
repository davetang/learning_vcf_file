#!/usr/bin/env bash

set -euo pipefail

tool=bwa
ver=0.7.17

wget https://github.com/lh3/bwa/releases/download/v${ver}/${tool}-${ver}.tar.bz2
tar -xjf ${tool}-${ver}.tar.bz2
cd ${tool}-${ver}
make
./bwa

>&2 echo Done
exit 0

