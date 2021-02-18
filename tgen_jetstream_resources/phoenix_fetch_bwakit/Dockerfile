ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_wget-1.21.1:${IMAGE_ARCH} as opt_wget
#FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_tar-1.34:${IMAGE_ARCH} as opt_tar
#FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bzip2-1.0.8:${IMAGE_ARCH} as opt_bzip2

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_wget /opt/ /opt/
#COPY --from=opt_tar /opt/ /opt/
#COPY --from=opt_bzip2 /opt/ /opt/

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
