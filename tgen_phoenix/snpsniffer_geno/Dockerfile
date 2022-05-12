ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gatk-4.2.2.0:${IMAGE_ARCH} AS opt_gatk
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} AS opt_bcftools

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_gatk /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/

ENV PATH /opt/java/bin/:/opt/miniconda/bin/:${PATH}
ENV JAVA_HOME /opt/java
ENV R_HOME /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

ENTRYPOINT ["/opt/gatk/run.sh"]

