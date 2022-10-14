FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

ARG HTSLIB_VER=1.16
RUN dnf upgrade -y && dnf install -y \ 
    tar \
    wget \
    autoconf \
    automake \
    make \
    gcc \
    perl-Data-Dumper \
    zlib-devel \
    bzip2 \
    bzip2-devel \
    xz-devel \
    curl-devel \
    openssl-devel &&\
    dnf clean all


RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 &&\
    tar jxf htslib-${HTSLIB_VER}.tar.bz2 &&\
    cd htslib-${HTSLIB_VER} &&\
    ./configure &&\
    make &&\ 
    make install &&\
    cd / &&\
    rm -rf htslib-${HTSLIB_VER} &&\
    rm htslib-${HTSLIB_VER}.tar.bz2

# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/htslib:1.16