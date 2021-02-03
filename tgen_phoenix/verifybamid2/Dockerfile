ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2.genomic:${IMAGE_ARCH} AS opt_r
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_verifybamid2-1.0.6:${IMAGE_ARCH} AS opt_verifybamid2

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/
COPY --from=opt_verifybamid2 /opt/ /opt/
COPY --from=build /tools/bin/which /tools/bin/which

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
