FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

# Versions
ARG SUBREAD_VER=2.0.0

RUN dnf upgrade -y && dnf install -y wget gzip tar bzip2 make gcc zlib zlib-devel ncurses-devel bzip2-devel xz-devel libcurl-devel openssl-devel &&\
    dnf clean all 

# Download and install Subread
RUN wget https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VER}/subread-${SUBREAD_VER}-source.tar.gz -q --show-progress && \
    tar zxf subread-${SUBREAD_VER}-source.tar.gz && \
    rm subread-${SUBREAD_VER}-source.tar.gz && \
    cd subread-${SUBREAD_VER}-source/src && \
    make -f Makefile.Linux && \
    cd / && \
    cp -R /subread-${SUBREAD_VER}-source/bin/* /usr/local/bin && \
    rm -rf subread-${SUBREAD_VER}-source

# Version Checks
RUN echo "SUBREAD_VER=${SUBREAD_VER}"

# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/subread:2.0.0.0