ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_pairoscope-0.4.2:${IMAGE_ARCH} AS opt_pairoscope
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_python /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install pandas

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_pairoscope /opt/ /opt/
COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=build /opt/ /opt/

ENV JS_PIPELINE_PATH /opt/phoenix/
COPY mm_igtx_pairoscope_calling_b38_356362b.py /opt/phoenix/required_scripts/

RUN sed -i 's|<reset-dirs />||g' /opt/share/fontconfig/conf.avail/05-reset-dirs-sample.conf

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
