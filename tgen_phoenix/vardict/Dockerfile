ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_vardict-1.8.1:${IMAGE_ARCH} AS opt_vardict
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bedtools-2.29.2:${IMAGE_ARCH} AS opt_bedtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} as opt_bcftools

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_vardict /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=build /tools/bin/which /tools/bin/which

ENV JAVA_HOME=/opt/java/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
