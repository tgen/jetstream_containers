ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2.genomic:${IMAGE_ARCH} as opt_r
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_samtools-1.10:${IMAGE_ARCH} as opt_samtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bedtools-2.29.2:${IMAGE_ARCH} as opt_bedtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_bootstrap as bootstrap

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/
COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY --from=bootstrap /usr/share/zoneinfo /usr/share/zoneinfo

RUN mkdir -p /opt/required_scripts/

COPY summarize_samstats_b3ddd38.R /opt/required_scripts/

ENV JS_PIPELINE_PATH /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

RUN mkdir -p /tools && ln -s /opt/bin /tools/
