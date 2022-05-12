ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_r-4.1.2.genomic:${IMAGE_ARCH} AS opt_r
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_samtools-1.13:${IMAGE_ARCH} AS opt_samtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bedtools-2.29.2:${IMAGE_ARCH} AS opt_bedtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_htslib-1.13:${IMAGE_ARCH} AS opt_htslib
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/
COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY --from=opt_htslib /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/

COPY tgen_mutation_burden.sh mutation_burden_per_arm.R /opt/bin/
RUN chmod a+x /opt/bin/tgen_mutation_burden.sh

ENV R_HOME /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

