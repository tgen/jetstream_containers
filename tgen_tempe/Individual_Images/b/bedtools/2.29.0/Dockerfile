FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

# Versions
ARG BEDTOOLS_VER=2.29.0
ARG BAMTOOLS_VER=2.5.1

#Install dependencies BamTools
RUN dnf upgrade -y && dnf install -y \
    wget \
    tar \
    cmake \ 
    gcc-c++ \ 
    json-glib-devel \
    pkgconf-pkg-config \
    zlib-devel \
#Installing dependencies for BEDTools
    python38 \
    bzip2-devel \
    xz-devel &&\
    dnf clean all 

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

# Version Checks
RUN echo "### BEDTOOLS_VER=${BEDTOOLS_VER} ### " &&\
    echo "### BAMTOOLS_VER=${BAMTOOLS_VER} ###"

CMD [ "bedtools", "--help" ]
# ghcr.io/tgen/jetstream_containers/bedtools:2.29.0
