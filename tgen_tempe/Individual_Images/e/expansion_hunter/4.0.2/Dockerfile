FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ARG EXPANSION_HUNTER_VER=4.0.2
ENV SRC /usr/local/src


RUN dnf upgrade -y && dnf install -y \
    tar \
    wget \
    gcc \
    gcc-c++ \
    cmake \
    bzip2-devel \
    libcurl-devel \
    libstdc++-static \
    openssl-devel \
    xz-devel \
    zlib-devel &&\
    dnf clean all

WORKDIR $SRC
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/v${EXPANSION_HUNTER_VER}/ExpansionHunter-v${EXPANSION_HUNTER_VER}-linux_x86_64.tar.gz &&\
    tar zxf ExpansionHunter-v${EXPANSION_HUNTER_VER}-linux_x86_64.tar.gz &&\
    rm -rf ExpansionHunter-v${EXPANSION_HUNTER_VER}-linux_x86_64.tar.gz &&\
    mv ExpansionHunter-v${EXPANSION_HUNTER_VER}-linux_x86_64/bin/ExpansionHunter /usr/local/bin/ &&\
    mv ExpansionHunter-v4.0.2-linux_x86_64/variant_catalog/* /usr/local/lib/

CMD [ "ExpansionHunter" ]

# ghcr.io/tgen/jetstream_containers/expansion_hunter:4.0.2