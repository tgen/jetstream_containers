FROM ghcr.io/tgen/jetstream_containers/r-base:3.6.1

ENV PATH=/tgen_mutation_burden-1.2.3:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib/:$LD_LIBRARY_PATH

ARG SAMTOOLS_VER=1.10
ARG BAMTOOLS_VER=2.5.1
ARG BEDTOOLS_VER=2.29.0
ARG BCFTOOLS_VER=1.10.2
ARG HTSLIB_VER=1.10.2
ARG MUTATION_BURDEN_VER=1.2.3

COPY tgen_mutation_burden-1.2.3.tar.gz /
COPY module_install.sh /
COPY req_modules.csv /

# Installing dependencies for SAMtools
RUN dnf upgrade -y && dnf install -y \
    wget \
    tar \
    bzip2 \
    make \
    gcc \
    zlib-devel \
    ncurses-devel \
    bzip2-devel \
    xz-devel \
    libcurl-devel \
#Installing dependencies for BamTools
    cmake \
    gcc-c++ \
    json-glib-devel \
    pkgconf-pkg-config \
#Installing dependencies for BEDTools
    python38 \
# Install dependencies for BcftTools
    gsl-devel \
    perl-libs \
# Install dependencies for HTSlib
    autoconf \
    automake \
    perl-Data-Dumper \
    curl-devel \
    openssl-devel \ 
# Install dependencies for Mutation Burden
    parallel &&\
    dnf clean all

# Download and install SAMtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install

# Download and install BamTools
RUN wget https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VER}.tar.gz -q --show-progress &&\
    tar -zxf v${BAMTOOLS_VER}.tar.gz &&\
    rm v${BAMTOOLS_VER}.tar.gz &&\
    cd bamtools-${BAMTOOLS_VER} &&\
    mkdir build &&\
    cd build &&\
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. &&\
    make &&\
    make install &&\
    cd / &&\
    rm -rf bamtools-${BAMTOOLS_VER}

# Download and install BedTools
RUN alternatives --set python /usr/bin/python3.8 &&\ 
    wget https://github.com/arq5x/bedtools2/archive/v${BEDTOOLS_VER}.tar.gz -q --show-progress &&\
    tar -zxf v${BEDTOOLS_VER}.tar.gz &&\
    rm v${BEDTOOLS_VER}.tar.gz &&\
    cd bedtools2-${BEDTOOLS_VER} &&\
    make &&\
    cp bin/* /usr/local/bin &&\
    cd / &&\
    rm -rf bedtools2-${BEDTOOLS_VER}

# Install BCFTools
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2 &&\
    tar jxf bcftools-${BCFTOOLS_VER}.tar.bz2 &&\
    cd bcftools-${BCFTOOLS_VER} && ./configure --enable-libgsl && make && make install && cd / &&\
    rm bcftools-${BCFTOOLS_VER}.tar.bz2 &&\
    rm -rf bcftools-${BCFTOOLS_VER}

# Install HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 &&\
    tar jxf htslib-${HTSLIB_VER}.tar.bz2 &&\
    cd htslib-${HTSLIB_VER} &&\
    ./configure &&\
    make &&\ 
    make install &&\
    cd / &&\
    rm -rf htslib-${HTSLIB_VER} &&\
    rm htslib-${HTSLIB_VER}.tar.bz2

# Configure and Run Script
RUN chmod +x module_install.sh && \
    ./module_install.sh && \
    rm module_install.sh && \
    rm req_modules.csv

# Download and install tgen_mutation_burden
RUN tar zxf tgen_mutation_burden-${MUTATION_BURDEN_VER}.tar.gz &&\
    rm tgen_mutation_burden-${MUTATION_BURDEN_VER}.tar.gz &&\
    cd tgen_mutation_burden-${MUTATION_BURDEN_VER} &&\
    chmod 755 tgen_mutation_burden.sh

# ghcr.io/tgen/jetstream_containers/mutation-burden:1.2.3 