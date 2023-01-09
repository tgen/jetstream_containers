FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

# Versions
ARG MANTA_VER=1.6.0

# Install dependencies
RUN dnf upgrade -y && dnf install -y python2 make cmake wget tar bzip2 boost-devel gcc gcc-c++ libstdc++-static zlib-devel htslib-tools &&\
    dnf clean all

RUN alternatives --set python /usr/bin/python2

# Download and install manta
RUN wget https://github.com/Illumina/manta/releases/download/v${MANTA_VER}/manta-${MANTA_VER}.release_src.tar.bz2 -q --show-progress && \
    tar -jxf manta-${MANTA_VER}.release_src.tar.bz2 && \
    mkdir build && \
    cd build && \
    ../manta-${MANTA_VER}.release_src/configure --jobs=4 --prefix=/usr/local && \
    make -j4 install && \
    cd / && \
    rm -rf build manta-${MANTA_VER}.release_src manta-${MANTA_VER}.release_src.tar.bz2 manta-${MANTA_VER}.release_src manta-${MANTA_VER}.release_src.tar.bz2

CMD ["python", "/usr/local/bin/configManta.py", "--help"]
# ghcr.io/tgen/jetstream_containers/manta:1.6.0

