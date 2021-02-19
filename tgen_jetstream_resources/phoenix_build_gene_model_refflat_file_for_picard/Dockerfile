ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bzip2-1.0.8:${IMAGE_ARCH} as opt_bzip2
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_curl-7.70.0:${IMAGE_ARCH} as opt_curl
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gawk-5.0.1:${IMAGE_ARCH} as opt_gawk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_krb5-1.19:${IMAGE_ARCH} as opt_krb5

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build

RUN mkdir -p /opt/bin/ && cd /opt/bin/ && \
    ( wget http://localhost:8000/gtfToGenePred || \
      wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred ) && \
    chmod a+x gtfToGenePred

COPY --from=opt_bzip2 /opt/ /opt/

RUN ln -s libbz2.so.1.0 /opt/lib/libbz2.so.1

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/
COPY --from=opt_curl /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=opt_krb5 /opt/ /opt/

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

