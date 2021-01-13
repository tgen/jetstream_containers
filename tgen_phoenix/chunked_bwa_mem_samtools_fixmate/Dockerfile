ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_samtools-1.10:${IMAGE_ARCH} as opt_samtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bwa-0.7.17:${IMAGE_ARCH} as opt_bwa
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gawk-5.0.1:${IMAGE_ARCH} as opt_gawk

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_bwa /opt/bin/bwa /opt/bin/bwa
COPY --from=opt_gawk /opt/ /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
