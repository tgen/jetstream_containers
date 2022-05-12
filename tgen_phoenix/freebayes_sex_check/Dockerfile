ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_freebayes-1.3.2:${IMAGE_ARCH} as opt_freebayes
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_samtools-1.13:${IMAGE_ARCH} as opt_samtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} as opt_bcftools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} as opt_gawk
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} as opt_grep

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_freebayes /opt/ /opt/
COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

