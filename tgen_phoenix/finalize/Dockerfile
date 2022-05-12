ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_freetype-2.11.0:${IMAGE_ARCH} AS opt_freetype
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_libpng-1.6.37:${IMAGE_ARCH} AS opt_libpng

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} AS build

COPY --from=opt_python /opt /opt
COPY --from=opt_freetype /opt /opt
COPY --from=opt_libpng /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install numpy
RUN python -m pip install multiqc

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=build /opt /opt

RUN mkdir -p /usr/bin/ && echo "#!/bin/bash" > /usr/bin/module && chmod a+x /usr/bin/module

