ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_r-4.1.2.genomic:${IMAGE_ARCH} AS opt_r
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/
COPY --from=opt_python /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/

RUN mkdir -p /tools/ && ln -s /opt/bin /tools/

ENV PATH /opt/java/bin/:/opt/miniconda/bin/:${PATH}
#ENV LD_LIBRARY_PATH /opt/lib/
ENV JAVA_HOME /opt/java
ENV R_HOME /opt

COPY seg_extend_229b8c7.py /opt/phoenix/required_scripts/
COPY plotCNVplus_4d89cb4.R /opt/phoenix/required_scripts/
ENV JS_PIPELINE_PATH /opt/phoenix/

RUN sed -i 's|<reset-dirs />||g' /opt/share/fontconfig/conf.avail/05-reset-dirs-sample.conf

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

ENTRYPOINT ["/opt/gatk/run.sh"]


