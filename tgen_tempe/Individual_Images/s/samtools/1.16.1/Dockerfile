FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19 as build1

ARG SAMTOOLS_VER=1.16.1

# Update installed packages and install requirements for SAMtools
RUN dnf upgrade -y && \
    dnf install -y wget \
                   tar \
                   bzip2 \
                   make \
                   gcc \
                   zlib-devel \
                   ncurses-devel \
                   bzip2-devel \
                   xz-devel \
                   libcurl-devel \
                   openssl-devel

# Download and install SAMtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install

# Finalize container without the build tools
FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

COPY --from=build1 /usr/local/ /usr/local

RUN dnf upgrade -y && \
    dnf install -y zlib \
                    ncurses-libs \
                    bzip2-libs \
                    xz-libs \
                    libcurl \
                    openssl-libs

# ghcr.io/tgen/jetstream_containers/samtools:1.16.1 
