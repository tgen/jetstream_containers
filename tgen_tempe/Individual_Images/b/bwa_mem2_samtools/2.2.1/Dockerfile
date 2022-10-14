FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19 

ENV PATH=$PATH:/usr/local/bin/bwa-mem2

ARG SAMTOOLSVER=1.15
ARG BWA_MEM2_VER=2.2.1

# Update installed packages and install requirements for SAMtools
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
    openssl-devel \
    zlib \
    ncurses-libs \
    bzip2-libs \
    xz-libs \
    libcurl \
    openssl-libs &&\
    dnf clean all

# Download and install SAMtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLSVER}.tar.bz2 && \
    rm samtools-${SAMTOOLSVER}.tar.bz2 && \
    cd samtools-${SAMTOOLSVER} && \
    ./configure && \
    make && \
    make install

# Install BWA-Mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VER}/bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 &&\
    tar jxf bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 &&\
    mv bwa-mem2-${BWA_MEM2_VER}_x64-linux bwa-mem2 &&\
    mv /bwa-mem2/bwa* /usr/local/bin/


WORKDIR /bwa-mem2 

CMD [ "bwa-mem2" ]
# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/bwa_mem2_samtools:2.2.1





