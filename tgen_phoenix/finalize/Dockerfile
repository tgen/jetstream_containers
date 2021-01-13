ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-3.8.6.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_freetype-2.10.2:${IMAGE_ARCH} AS opt_freetype
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_libpng-1.6.37:${IMAGE_ARCH} AS opt_libpng

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

COPY --from=opt_python /opt /opt
COPY --from=opt_freetype /opt /opt
COPY --from=opt_libpng /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install numpy
RUN python -m pip install multiqc

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=build /opt /opt

RUN mkdir -p /usr/bin/ && echo "#!/bin/bash" > /usr/bin/module && chmod a+x /usr/bin/module

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
