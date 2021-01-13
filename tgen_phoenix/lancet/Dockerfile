ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bcftools-1.10.2:${IMAGE_ARCH} as opt_bcftools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_lancet-1.1.0:${IMAGE_ARCH} as opt_lancet

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_bcftools /opt/ /opt/
COPY --from=opt_lancet /opt/ /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
