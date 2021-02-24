ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gatk-4.1.4.1:${IMAGE_ARCH} as opt_gatk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bedtools-2.29.2:${IMAGE_ARCH} as opt_bedtools

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_gatk /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY capture_kits.csv /opt/phoenix/
COPY make_exome_refpack.py /opt/utility_scripts/
COPY formats /opt/utility_scripts/pylib/formats

ENV PATH /opt/java/bin/:/opt/miniconda/bin/:${PATH}
ENV JAVA_HOME /opt/java
ENV PYTHONPATH /opt/utility_scripts/pylib/

RUN mkdir -p /tools/ && ln -s /opt/bin /tools/

ENTRYPOINT ["/opt/gatk/run.sh"]

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

#RUN which python && ls /opt/utility_scripts/pylib/
#RUN python /opt/utility_scripts/make_exome_refpack.py && false