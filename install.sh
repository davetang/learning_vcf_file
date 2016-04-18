#!/bin/bash

# just clone the Bpipe repository
wget http://download.bpipe.org/versions/bpipe-0.9.9.tar.gz
tar -xzf bpipe-0.9.9.tar.gz
ln -s bpipe-0.9.9/bin/bpipe

# download and install BWA
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

# download and install HTSlib
wget https://github.com/samtools/htslib/releases/download/1.3/htslib-1.3.tar.bz2
tar -xjf htslib-1.3.tar.bz2
cd htslib-1.3
make
make prefix=. install
cd ..

# download and install SAMTools
wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
tar -xjf samtools-1.3.tar.bz2
cd samtools-1.3
make
make prefix=. install
cd ..

# download and install BCFtools
wget https://github.com/samtools/bcftools/releases/download/1.3/bcftools-1.3.tar.bz2
tar -xjf bcftools-1.3.tar.bz2
cd bcftools-1.3
make
make prefix=. install
cd ..

