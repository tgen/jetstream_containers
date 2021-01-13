ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_openjdk-8u:${IMAGE_ARCH} AS opt_openjdk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bcftools-1.10.2:${IMAGE_ARCH} AS opt_bcftools

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

ARG LOCAL_SOURCES_SERVER

RUN mkdir -p /tmp/scratch/build && \
    cd /tmp/scratch && wget -t 1 ${LOCAL_SOURCES_SERVER}/snpEff_v4_3t_core.zip || \
                       wget -t 3 http://downloads.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip && \
                       echo "d55a7389a78312947c1e7dadf5e6897b42d3c6e942e7c1b8ec68bb35d2ae2244  snpEff_v4_3t_core.zip" | sha256sum -c && \
    cd /tmp/scratch && unzip snpEff_v4_3t_core.zip && mv snpEff /opt/ && \
    cd / && rm -rf /tmp/scratch

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_openjdk /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=build /opt /opt

RUN mkdir -p /opt/bin/ && echo -e '#!/bin/bash\njava -Xmx8G -jar /opt/snpEff/snpEff.jar $@' > /opt/bin/snpEff && chmod +x /opt/bin/snpEff

ENV PATH /opt/java/bin/:/opt/miniconda/bin/:${PATH}
ENV JAVA_HOME /opt/java

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
