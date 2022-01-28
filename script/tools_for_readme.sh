#!/usr/bin/env bash

set -euo pipefail

# path preceeding this script
relpath=$(dirname $0)

if [[ -x $(command -v readlink) ]]; then
   root=$(readlink -f ${relpath}/..)
elif [[ -x $(command -v readpath) ]]; then
   root=$(realpath ${relpath}/..)
else
   >&2 echo Could not convert ${relpath}/.. to an absolute path
   exit 1
fi

if [[ ! -d ${root}/bin ]]; then
   mkdir ${root}/bin
fi

for tool in htslib bcftools; do
   ver=1.14
   check=${tool}
   if [[ ${tool} == htslib ]]; then
      check=bgzip
   fi
   if [[ ! -e ${root}/bin/${check} ]]; then
      url=https://github.com/samtools/${tool}/releases/download/${ver}/${tool}-${ver}.tar.bz2
      wget ${url}
      tar xjf ${tool}-${ver}.tar.bz2
      cd ${tool}-${ver}
      ./configure --prefix=${root}
      make && make install
      cd ..
      rm -rf ${tool}-${ver}*
   fi
done

if [[ ! -e ${root}/bin/vcf2bed ]]; then
   cd ${root}
   bedops_ver=2.4.40
   url=https://github.com/bedops/bedops/releases/download/v${bedops_ver}/bedops_linux_x86_64-v${bedops_ver}.tar.bz2
   wget ${url}
   tar -xjf bedops_linux_x86_64-v${bedops_ver}.tar.bz2
   rm bedops_linux_x86_64-v${bedops_ver}.tar.bz2
fi

if [[ ! -e ${root}/bin/gh-md-toc ]]; then
   cd ${root}
   toc_ver=0.8.0
   url=https://github.com/ekalinin/github-markdown-toc/archive/refs/tags/${toc_ver}.tar.gz
   wget -O github-markdown-toc-${toc_ver}.tar.gz ${url}
   tar xzf github-markdown-toc-${toc_ver}.tar.gz
   mv github-markdown-toc-${toc_ver}/gh-md-toc ${root}/bin
   rm -rf github-markdown-toc-${toc_ver}*
fi

if [[ ! -e ${root}/bin/vt ]]; then
   cd ${root}
   vt_ver=0.57721
   url=https://github.com/atks/vt/archive/refs/tags/${vt_ver}.tar.gz
   wget ${url} -O vt-${vt_ver}.tar.gz
   tar -xzf vt-${vt_ver}.tar.gz
   cd vt-${vt_ver}
   make
   mv vt ${root}/bin/
   cd ..
   rm -rf vt-${vt_ver}*
fi

for t in bcftools vcf2bed bgzip tabix gh-md-toc vt; do
   if [[ ! -L ${root}/${t} ]]; then
      ln -s ${root}/bin/${t} .
   fi
done

>&2 echo Done
exit 0

