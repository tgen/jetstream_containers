ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_zlib-1.2.11:${IMAGE_ARCH} AS opt_zlib

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/subread-2.0.1-Linux-x86_64.tar.gz || \
                       wget -t 3 https://downloads.sourceforge.net/project/subread/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz && \
                       echo "bd4d3834ea4690d27c8b06095662343ebdabd3c24c7e6c09611bfeb61ef041fa  subread-2.0.1-Linux-x86_64.tar.gz" | sha256sum -c && \
    mkdir -p /opt && cd /opt && tar -xvf /tmp/scratch/subread-2.0.1-Linux-x86_64.tar.gz && \
    cd / && rm -rf /tmp/scratch

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/
COPY --from=opt_zlib /opt /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

ENV PATH /opt/subread-2.0.1-Linux-x86_64/bin/:${PATH}

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
