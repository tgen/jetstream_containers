ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-3.8.6.headless:${IMAGE_ARCH} AS opt_python

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_python /opt/ /opt/

ENV PATH_TO_PHOENIX_REPO /opt/phoenix/
COPY samStats2json.py /opt/phoenix/required_scripts/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module



LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
