ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_vardict-1.8.1:${IMAGE_ARCH} AS opt_vardict
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bedtools-2.29.2:${IMAGE_ARCH} AS opt_bedtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bcftools-1.10.2:${IMAGE_ARCH} as opt_bcftools

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_vardict /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/

ENV JAVA_HOME=/opt/java/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
