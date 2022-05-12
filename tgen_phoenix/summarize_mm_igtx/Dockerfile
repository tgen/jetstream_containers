ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_r-4.1.2.genomic:${IMAGE_ARCH} AS opt_r

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_r /opt/ /opt/

RUN mkdir -p /tools && ln -s /opt/bin /tools/

ENV JS_PIPELINE_PATH /opt/phoenix/
COPY summarize_Ig_875a823.R /opt/phoenix/required_scripts/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

