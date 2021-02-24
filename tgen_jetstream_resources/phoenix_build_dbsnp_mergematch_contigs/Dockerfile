ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2.genomic:${IMAGE_ARCH} as opt_r

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/

RUN mkdir -p /tools && ln -s /opt/bin /tools/

COPY MergeMatch_dbSNP_GRCh38_Contigs.R /opt/utility_scripts/

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

