ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_zlib-1.2.11:${IMAGE_ARCH} AS opt_zlib
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-2.7.18.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bcftools-1.10.2:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_samtools-1.10:${IMAGE_ARCH} AS opt_samtools

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

ARG LOCAL_SOURCES_SERVER=http://localhost:8000

COPY --from=opt_zlib /opt /opt
COPY --from=opt_python /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch/ && wget -t 1 http://localhost:8000/strelka-2.9.10.release_src.tar.bz2 || \
                        wget -t 3 https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.release_src.tar.bz2 && \
    cd /tmp/scratch/ && tar -xvf strelka-2.9.10.release_src.tar.bz2 && \
    cd /tmp/scratch/build/ && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt ../strelka-2.9.10.release_src/ && make -j$(nproc) && make install && \
    cd / && rm -rf /tmp/scratch

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt /opt
COPY --from=opt_bcftools /opt /opt
COPY --from=opt_samtools /opt /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
