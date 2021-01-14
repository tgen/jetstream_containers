ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/salmon_0.14.2_linux_x86_64.tar.gz || \
                       wget -t 3 https://github.com/COMBINE-lab/salmon/releases/download/v0.14.2/salmon_0.14.2_linux_x86_64.tar.gz && \
                       echo "d6973db43318375da6bb59ad86711c2edb874733f3e9b875225412afc20ce2ad  salmon_0.14.2_linux_x86_64.tar.gz" | sha256sum -c && \
    mkdir -p /opt && cd /opt/ && tar -xf /tmp/scratch/salmon_0.14.2_linux_x86_64.tar.gz && \
    cd / && rm -rf /tmp/scratch

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/

ENV PATH /opt/salmon-latest_linux_x86_64/bin/:${PATH}
ENV LD_LIBRARY_PATH /opt/salmon-latest_linux_x86_64/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module


LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
