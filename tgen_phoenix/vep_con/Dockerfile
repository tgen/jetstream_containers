ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_vep-106:${IMAGE_ARCH} AS opt_vep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gzip-1.11:${IMAGE_ARCH} AS opt_gzip
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_tar-1.34:${IMAGE_ARCH} AS opt_tar

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

COPY --from=opt_vep /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=opt_gzip /opt/ /opt/
COPY --from=opt_tar /opt/ /opt/

ENV PATH /opt/bin/:/usr/bin/:/bin/:/usr/sbin/:/sbin/
ENV LD_LIBRARY_PATH /opt/lib/
ENV PERL5LIB /opt/ensembl-vep/modules/:/opt/ensembl-vep/

