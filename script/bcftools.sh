#!/usr/bin/env bash

set -euo pipefail

ver=1.13
tool=bcftools
url=https://github.com/samtools/${tool}/releases/download/${ver}/${tool}-${ver}.tar.bz2

wget ${url}
tar xjf ${tool}-${ver}.tar.bz2
cd ${tool}-${ver}
./configure
make && make install
cd ..

rm -rf ${tool}-${ver} ${tool}-${ver}.tar.bz2

>&2 echo Done
exit 0

