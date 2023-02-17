FROM debian:bullseye-20230208-slim as build
ARG BCFTOOLS_VER=1.16

# Install Dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    wget \
    tar \
    ca-certificates \
    bzip2 \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libperl-dev \
    libgsl0-dev

# Install Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    tar jxf bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VER} && \
    ./configure --enable-libgsl && \
    make && \
    make install && \
    cd .. && \
    rm bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    rm -rf bcftools-${BCFTOOLS_VER}

# Start 2nd Stage
FROM debian:bullseye-20230208-slim

# Copy Bcftools Application from 1st Build Stage
COPY --from=build /usr/local/ /usr/local

# Install Bcftools Required Packages and Clean Cache
RUN apt-get update && apt-get install --no-install-recommends -y \
    libcurl3-gnutls \
    libgsl25 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoclean

# Create Working Directory
WORKDIR /data

# ghcr.io/tgen/jetstream_containers/bcftools:1.16-optimized