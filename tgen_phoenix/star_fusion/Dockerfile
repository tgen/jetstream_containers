ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_perl-5.32.0:${IMAGE_ARCH} AS opt_perl
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_db-6.2.32:${IMAGE_ARCH} AS opt_db
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_star-2.7.3a:${IMAGE_ARCH} AS opt_star
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_openjdk-8u:${IMAGE_ARCH} AS opt_openjdk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-3.8.6.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bowtie2-2.4.2:${IMAGE_ARCH} AS opt_bowtie2
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_samtools-1.10:${IMAGE_ARCH} AS opt_samtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_grep-3.4:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bzip2-1.0.8:${IMAGE_ARCH} AS opt_bzip2
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gmap-2017.11.15:${IMAGE_ARCH} AS opt_gmap

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build_1

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_perl /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/STAR-Fusion-v1.8.1.FULL.tar.gz || \
                       wget -t 3 https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v1.8.1/STAR-Fusion-v1.8.1.FULL.tar.gz && \
                       echo "daef3238904c9aafe197559b09fd5fc37034b2e38ae993c403c5c557b9c5623e  STAR-Fusion-v1.8.1.FULL.tar.gz" | sha256sum -c && \
    mkdir -p /opt && cd /opt && tar -xvf /tmp/scratch/STAR-Fusion-v1.8.1.FULL.tar.gz && \
    cd / && rm -rf /tmp/scratch

COPY --from=opt_perl /opt/ /opt/
COPY --from=opt_db /opt/ /opt/

RUN mkdir -p /usr/local/BerkeleyDB/ && ln -s /opt/include /usr/local/BerkeleyDB/ && \
                                       ln -s /opt/lib /usr/local/BerkeleyDB/ && \
    cpan App::cpanminus && \
    cpanm DB_File && \
    cpanm Set::IntervalTree && \
    cpanm JSON::XS && \
    cpanm Carp::Assert && \
    cpanm URI::Escape

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/jellyfish-2.3.0.tar.gz || \
                       wget -t 3 https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && \
                       echo "e195b7cf7ba42a90e5e112c0ed27894cd7ac864476dc5fb45ab169f5b930ea5a  jellyfish-2.3.0.tar.gz" | sha256sum -c && \
    mkdir -p /tmp/scratch/ && tar -xf jellyfish-2.3.0.tar.gz && \
    cd /tmp/scratch/jellyfish-2.3.0/ && CFLAGS="-O3 -march=${IMAGE_ARCH}" \
                                        CXXFLAGS="-O3 -march=${IMAGE_ARCH}" \
                                        ./configure --prefix=/opt && \
                                        make -j$(nproc) && \
                                        make install && \
    cd / && rm -rf /tmp/scratch

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/salmon_0.14.2_linux_x86_64.tar.gz || \
                       wget -t 3 https://github.com/COMBINE-lab/salmon/releases/download/v0.14.2/salmon_0.14.2_linux_x86_64.tar.gz && \
                       echo "d6973db43318375da6bb59ad86711c2edb874733f3e9b875225412afc20ce2ad  salmon_0.14.2_linux_x86_64.tar.gz" | sha256sum -c && \
    mkdir -p /opt && cd /opt/ && tar -xf /tmp/scratch/salmon_0.14.2_linux_x86_64.tar.gz && \
    cd / && rm -rf /tmp/scratch

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build_2

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

RUN mkdir -p /opt && \
    cd /opt/ && git clone --recursive https://github.com/trinityrnaseq/trinityrnaseq.git && \
    cd /opt/trinityrnaseq && git checkout v2.11.0 && \
    cd /opt/trinityrnaseq && \
        sed -i 's|#include "htslib/bgzf.h"|#include "htslib/bgzf.h"\n#include <string>|g' trinity-plugins/bamsifter/sift_bam_max_cov.cpp && \
    cd /opt/trinityrnaseq && CFLAGS="-O3 -march=${IMAGE_ARCH}" \
                             CXXFLAGS="-O3 -march=${IMAGE_ARCH}" \
                             make -j$(nproc)

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build_3

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_python /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install numpy
RUN python -m pip install igv_reports

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_star /opt/ /opt/
COPY --from=opt_openjdk /opt/ /opt/
COPY --from=opt_bowtie2 /opt/ /opt/
COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_bzip2 /opt/ /opt/
COPY --from=opt_gmap /opt/ /opt/
COPY --from=build_1 /opt/ /opt/
COPY --from=build_2 /opt/trinityrnaseq /opt/trinityrnaseq
COPY --from=build_3 /opt/ /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

ENV PATH /opt/java/bin/:/opt/STAR-Fusion-v1.8.1/:/opt/trinityrnaseq/:/opt/salmon-latest_linux_x86_64/bin/:${PATH}
ENV JAVA_HOME /opt/java
ENV TRINITY_HOME /opt/trinityrnaseq

RUN mkdir -p /tools/ && ln -s /opt/bin /tools/
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
