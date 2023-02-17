FROM debian:bullseye-20230208-slim as build

ARG BWA_MEM2_VER=2.2.1
ARG SAMTOOLS_VER=1.16.1

# Install Dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    perl \
    automake \
    autoconf \
    wget \
    ca-certificates \
    tar \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \ 
    libncurses5-dev

# Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install

# Install Bwa-Mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VER}/bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 && \
    tar jxf bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 && \
    rm bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 && \
    mv bwa-mem2-${BWA_MEM2_VER}_x64-linux bwa-mem2 && \
    mv /bwa-mem2/bwa* /usr/local/bin/

# Starting 2nd Stage
FROM debian:bullseye-20230208-slim

# Copy Bwa-Mem2/Samtools Application from 1st Build Stage
COPY --from=build /usr/local/ /usr/local

# Install Bwa-Mem2/Samtools Required Packages and Clean Cache
RUN apt-get update && apt-get install --no-install-recommends -y \
    libcurl3-gnutls \
    libncursesw6 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoclean

# Create Working Directory
WORKDIR /data

# ghcr.io/tgen/jetstream_containers/bwa_mem2_samtools:2.2.1-optimized