# syntax=docker/dockerfile:1

#sudo docker run --rm -v ~/code/SBSC/test_output:/data -it test_sbsc SBSCall.py --chrom chr17 --processes 4 --cancer_pile SBSC/tests/resources/HCC1937_t_RNF157.pileup.gz --normal_pile SBSC/tests/resources/HCC1937_n_RNF157.pileup.gz --ref SBSC/tests/resources/chr17_RNF157.fa --output ../data/out
#sudo docker build --tag test_sbsc .

#https://gist.github.com/jmarshall/adc75969306354f4398896f1fc1f4865
FROM python:3.11.1-slim

USER root

RUN apt-get update && apt-get install -y \
	build-essential \
	curl \
	git \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libperl-dev \
	libssl-dev \
	zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
ARG htsversion=1.9
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig
RUN curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install)
RUN curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install)
# RUN git clone --depth 1 git://github.com/samtools/htslib-plugins && \
#     (cd htslib-plugins && make PLUGINS='hfile_cip.so hfile_mmap.so' install)

WORKDIR /app

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN mkdir SBSC
COPY . SBSC

COPY requirements_dev.txt requirements_dev.txt
RUN pip3 install -r requirements_dev.txt

RUN mkdir ../data

#CMD [ "pip3", "", ""]