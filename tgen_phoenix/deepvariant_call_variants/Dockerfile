ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/ubuntu-16.04_opt_deepvariant-0.10.0:${IMAGE_ARCH} 

ENV PATH /opt/deepvariant/bin/:/opt/bin/:${PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /usr/bin/ && echo "#!/bin/bash" > /usr/bin/module && chmod a+x /usr/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
