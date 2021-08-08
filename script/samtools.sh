#!/usr/bin/env bash

set -euo pipefail

ver=1.13
tool=samtools
url=https://github.com/samtools/samtools/releases/download/${ver}/${tool}-${ver}.tar.bz2
dir=/home/ngs/Analysis/tools/${tool}-${ver}

wget ${url}
tar xjf ${tool}-${ver}.tar.bz2
cd ${tool}-${ver}
./configure
make && make install
cd ..

rm -rf ${tool}-${ver} ${tool}-${ver}.tar.bz2

>&2 echo Done
exit 0

