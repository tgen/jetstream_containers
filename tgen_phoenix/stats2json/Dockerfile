ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_python /opt/ /opt/

COPY samStats2json_3a90a2f.py /opt/phoenix/required_scripts/
ENV JS_PIPELINE_PATH /opt/phoenix/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

