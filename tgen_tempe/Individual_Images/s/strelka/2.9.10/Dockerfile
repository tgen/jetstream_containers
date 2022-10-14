FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

# Versions
ARG STRELKA_VER=2.9.10

RUN dnf upgrade -y && dnf install -y python2 wget make bzip2 libstdc++-static gcc cmake gcc-c++ zlib-devel boost-devel tar &&\
    dnf clean all

# Download and install strelka
RUN wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_VER}/strelka-${STRELKA_VER}.release_src.tar.bz2 -q --show-progress && \
    tar -xjf strelka-${STRELKA_VER}.release_src.tar.bz2 && \
    mkdir build && \
    cd build && \
    ../strelka-${STRELKA_VER}.release_src/configure --jobs=4 --prefix=/usr/local && \
    make -j4 install && \
    cd / &&\
    rm -rf build strelka-${STRELKA_VER}.release_src strelka-${STRELKA_VER}.release_src.tar.bz2

# Version Checks
RUN echo "STRELKA_VER=${STRELKA_VER}"

# CMD ["bash", "usr/local/bin/runStrelkaSomaticWorkflowDemo.bash"]
# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/strelka:2.9.10