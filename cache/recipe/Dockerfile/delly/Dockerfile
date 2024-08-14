# use the ubuntu base image
FROM ubuntu:22.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install delly
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make STATIC=1 PARALLEL=1 all \
    && make install

# install htslib, samtools, bcftools
RUN apt-get update && apt-get install -y wget && apt-get clean\
&& cd /opt \
&& wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
&& wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
&& wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
&& cd /opt && tar jxvf htslib-1.20.tar.bz2 && cd /opt/htslib-1.20 && ./configure && make && make install \
&& cd /opt && tar jxvf samtools-1.20.tar.bz2 && cd /opt/samtools-1.20 && ./configure && make && make install \
&& cd /opt && tar jxvf bcftools-1.20.tar.bz2 && cd /opt/bcftools-1.20 && ./configure && make && make install \
&& cd /opt && rm -rf /opt/htslib-1.20 /opt/samtools-1.20 /opt/bcftools-1.20 htslib-1.20.tar.bz2 samtools-1.20.tar.bz2 bcftools-1.20.tar.bz2


# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/delly/bin
WORKDIR /opt/delly/bin
COPY --from=0 /opt/delly/bin/delly .

# Workdir
WORKDIR /root/

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

# Set OpenMP threads
ENV OMP_NUM_THREADS 32

# by default /bin/sh is executed
CMD ["/bin/sh"]
