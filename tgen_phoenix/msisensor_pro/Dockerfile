ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_msisensorpro-1.2.0:${IMAGE_ARCH} AS opt_msisensorpro
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_samtools-1.13:${IMAGE_ARCH} AS opt_samtools

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_msisensorpro /opt /opt
COPY --from=opt_samtools /opt /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
