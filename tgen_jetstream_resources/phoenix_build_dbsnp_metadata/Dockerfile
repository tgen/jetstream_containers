ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_edirect-20210219:${IMAGE_ARCH} as opt_edirect
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_grep-3.4:${IMAGE_ARCH} as opt_grep
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gawk-4.2.1:${IMAGE_ARCH} as opt_gawk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_sed-4.8:${IMAGE_ARCH} as opt_sed
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_wget-1.21.1:${IMAGE_ARCH} as opt_wget

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build

RUN mkdir -p /opt && \
    cd /opt && \
    wget https://github.com/step-/JSON.awk/archive/1.3.tar.gz && \
    tar -xvf 1.3.tar.gz && \
    rm 1.3.tar.gz
   
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_edirect /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=opt_sed /opt/ /opt/
COPY --from=opt_wget /opt/ /opt/
COPY --from=build /opt/ /opt/

ENV PATH ${PATH}:/opt/edirect/

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

