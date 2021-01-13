ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gatk-4.1.4.1:${IMAGE_ARCH} AS opt_gatk

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_gatk /opt/ /opt/

ENV PATH /opt/java/bin/:/opt/miniconda/bin/:${PATH}
#ENV LD_LIBRARY_PATH /opt/lib/
ENV JAVA_HOME /opt/java
#ENV R_HOME /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

RUN mkdir -p /tools/ && ln -s /opt/bin /tools/

ENTRYPOINT ["/opt/gatk/run.sh"]

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
