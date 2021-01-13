ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2.headless.genomic:${IMAGE_ARCH} AS opt_r

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build

COPY --from=opt_r /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /tmp && \
    mkdir -p ~/.R/ && \
    echo "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")" | R --no-save && \
    cd / && rm -rf /tmp/* && rm -rf ~/.R

RUN mkdir -p /tmp && \
    mkdir -p ~/.R/ && \
    echo "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")" | R --no-save && \
    cd / && rm -rf /tmp/* && rm -rf ~/.R

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/
COPY --from=build /tools/bin/which /tools/bin/which

ENV PATH_TO_PHOENIX_REPO /opt/phoenix/
COPY runIchorCNA_47ce8db.R /opt/phoenix/required_scripts/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
