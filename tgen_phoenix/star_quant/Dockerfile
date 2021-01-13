ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=haswell
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_star-2.7.3a:${IMAGE_ARCH} AS opt_star

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_star /opt/ /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module


LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
