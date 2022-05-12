ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_snpsniffer-7.0.0:${IMAGE_ARCH} AS opt_snpsniffer
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_rsync-3.2.3:${IMAGE_ARCH} AS opt_rsync
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_snpsniffer /opt/ /opt/
COPY --from=opt_rsync /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/

ENV PATH /opt/java/bin/:/opt/snpSniffer/:${PATH}
ENV JAVA_HOME /opt/java
ENV R_HOME /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

