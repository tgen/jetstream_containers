FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19 as build1

ARG HTSLIB_VER=1.15.1
ARG OCTOPUS_VER=0.7.4

ENV LD_LIBRARY_PATH=/usr/local/lib/

RUN dnf upgrade -y && \
    dnf install -y wget \
                   tar \
                   bzip2 \
                   bzip2-devel \
                   zlib-devel \
                   libcurl-devel \
                   xz \
                   xz-devel \
                   openssl-devel \
                   python39 \
                   curl \
                   gcc \
                   make \
                   python39-devel \
                   git \
                   cmake \
                   gcc-c++ \
                   boost-devel \
                   gmp-devel && \
    # Set python 3.9 as default python version. 
    alternatives --set python /usr/bin/python3.9 && \
    alternatives --install /usr/bin/pip pip /usr/bin/pip3.9 1 && \
    # Install HSTlib.  Required for Octopus and pysam python module.
    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
    tar jxf htslib-${HTSLIB_VER}.tar.bz2 && \
    cd htslib-${HTSLIB_VER} && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    # Install octopus
    pip install distro && \
    wget https://github.com/luntergroup/octopus/archive/v${OCTOPUS_VER}.tar.gz && \
    tar xf v${OCTOPUS_VER}.tar.gz && \
    octopus-${OCTOPUS_VER}/scripts/install.py --threads 4 && \
    cd / && \
    mv octopus-${OCTOPUS_VER} octopus

FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ENV LD_LIBRARY_PATH=/usr/local/lib/
ENV PATH=/octopus/bin:$PATH

COPY --from=build1 /usr/local/ /usr/local/
COPY --from=build1 /octopus /octopus

RUN dnf upgrade -y && \
    dnf install -y bzip2-libs \
                   zlib \
                   libcurl \
                   xz \
                   xz-libs \
                   openssl-libs \
                   python39 \
                   boost \
                   gmp

# ghcr.io/tgen/jetstream_containers/octopus:0.7.4