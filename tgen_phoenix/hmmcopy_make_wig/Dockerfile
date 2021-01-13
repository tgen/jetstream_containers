ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_hmmcopy_utils-20170718:${IMAGE_ARCH} AS opt_hmmcopy_utils

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_hmmcopy_utils /opt /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

#RUN cat /etc/group && false
#RUN id && false
#RUN mkdir /tmp && chmod a+rwx /tmp

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
