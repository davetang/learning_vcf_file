#!/usr/bin/env bash

set -euo pipefail

dir=$(pwd)/tools
if [[ ! -d ${dir}/bin ]]; then
   mkdir -p ${dir}/bin
fi

ver=1.14
for tool in htslib bcftools samtools; do
   check=${tool}
   if [[ ${tool} == htslib ]]; then
      check=bgzip
   fi
   if [[ ! -e ${dir}/bin/${check} ]]; then
      url=https://github.com/samtools/${tool}/releases/download/${ver}/${tool}-${ver}.tar.bz2
      wget ${url}
      tar xjf ${tool}-${ver}.tar.bz2
      cd ${tool}-${ver}
      ./configure --prefix=${dir}
      make && make install
      cd ..
      rm -rf ${tool}-${ver}*
   fi
done

if [[ ! -e ${dir}/bin/bwa ]]; then
   wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.tar.gz
   tar -xzf v0.7.17.tar.gz
   cd bwa-0.7.17 && make && mv bwa ${dir}/bin/
   cd ..
   rm -rf bwa-0.7.17 v0.7.17.tar.gz
fi

freebayes_ver=1.3.6
if [[ ! -e ${dir}/bin/freebayes ]]; then
   cd ${dir}/bin
   wget https://github.com/freebayes/freebayes/releases/download/v${freebayes_ver}/freebayes-${freebayes_ver}-linux-amd64-static.gz
   gunzip freebayes-${freebayes_ver}-linux-amd64-static.gz
   chmod 755 freebayes-${freebayes_ver}-linux-amd64-static
   ln -s freebayes-${freebayes_ver}-linux-amd64-static freebayes
   cd ../..
fi

gatk_ver=4.2.4.1
if [[ ! -d ${dir}/gatk-${gatk_ver} ]]; then
   cd ${dir}
   wget https://github.com/broadinstitute/gatk/releases/download/${gatk_ver}/gatk-${gatk_ver}.zip
   unzip gatk-${gatk_ver}.zip
   ln -s gatk-${gatk_ver}/gatk .
   rm gatk-${gatk_ver}.zip
   cd ..
fi

bpipe_ver=0.9.11
if [[ ! -d ${dir}/bpipe-${bpipe_ver} ]]; then
   cd ${dir}
   wget https://github.com/ssadedin/bpipe/releases/download/${bpipe_ver}/bpipe-${bpipe_ver}.tar.gz
   tar -xzf bpipe-${bpipe_ver}.tar.gz
   ln -s bpipe-${bpipe_ver}/bin/bpipe .
   rm bpipe-${bpipe_ver}.tar.gz
fi

>&2 echo Done
exit 0

