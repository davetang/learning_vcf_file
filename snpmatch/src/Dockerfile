FROM continuumio/miniconda

MAINTAINER Dave Tang <me@davetang.org>

RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y build-essential \
    && apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install -c bioconda bcftools && pip install SNPmatch

