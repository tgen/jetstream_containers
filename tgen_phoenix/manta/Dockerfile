ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_manta-1.6.0:${IMAGE_ARCH} AS opt_manta
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_samtools-1.13:${IMAGE_ARCH} AS opt_samtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_htslib-1.13:${IMAGE_ARCH} AS opt_htslib

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_manta /opt /opt
COPY --from=opt_bcftools /opt /opt
COPY --from=opt_samtools /opt /opt
COPY --from=opt_htslib /opt /opt

ENV PATH /opt/miniconda/bin/:${PATH}

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
