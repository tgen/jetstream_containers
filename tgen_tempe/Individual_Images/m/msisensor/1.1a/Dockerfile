FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

ENV PATH=/msisensor-pro-1.1.a/binary:$PATH

# Versions
ARG MSISENSOR_VER=1.1.a
ARG SAMTOOLS_VER=1.10

# Installing pre requisits for SAMtools
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
    openssl-devel &&\
    dnf clean all 

# Download and install SAMtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 &&\
    tar jxf samtools-${SAMTOOLS_VER}.tar.bz2 &&\
    rm samtools-${SAMTOOLS_VER}.tar.bz2 &&\
    cd samtools-${SAMTOOLS_VER} &&\
    ./configure &&\
    make &&\
    make install &&\
    cd / &&\
    rm -rf samtools-${SAMTOOLS_VER}

# Download and install MSIsensor
RUN wget https://github.com/xjtu-omics/msisensor-pro/archive/refs/tags/${MSISENSOR_VER}.tar.gz &&\
    tar zxf ${MSISENSOR_VER}.tar.gz &&\
    rm ${MSISENSOR_VER}.tar.gz

CMD [ "msisensor-pro","--help" ]
# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/msisensor:1.1.a