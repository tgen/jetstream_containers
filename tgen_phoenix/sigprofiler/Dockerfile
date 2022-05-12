ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_rsync-3.2.3:${IMAGE_ARCH} AS opt_rsync
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_tar-1.34:${IMAGE_ARCH} AS opt_tar

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_python /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install SigProfilerExtractor

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_rsync /opt/ /opt/
COPY --from=opt_tar /opt/ /opt/
COPY --from=build /opt/ /opt/

ENV JS_PIPELINE_PATH /opt/phoenix/
COPY sigprofiler_7595614.py /opt/phoenix/required_scripts/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
