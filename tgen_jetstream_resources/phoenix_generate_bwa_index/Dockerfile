ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_grep-3.4:${IMAGE_ARCH} as opt_grep
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_sed-4.8:${IMAGE_ARCH} as opt_sed
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bwa-0.7.17:${IMAGE_ARCH} as opt_bwa

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_sed /opt/ /opt/
COPY --from=opt_bwa /opt/ /opt/

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

