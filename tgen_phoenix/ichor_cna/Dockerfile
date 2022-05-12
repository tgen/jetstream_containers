ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_r-4.1.2.genomic:${IMAGE_ARCH} AS opt_r
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_poppler-21.01.0:${IMAGE_ARCH} AS opt_poppler
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

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

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/
COPY --from=build /tools/bin/which /tools/bin/which
COPY --from=opt_poppler /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/

COPY runIchorCNA_47ce8db.R /opt/phoenix/required_scripts/
ENV JS_PIPELINE_PATH /opt/phoenix/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

