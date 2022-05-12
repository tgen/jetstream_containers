ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_perl-5.34.0:${IMAGE_ARCH} as opt_perl

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_perl /opt/ /opt/
COPY annotSeg_7102f1c.pl /opt/required_scripts/
RUN chmod a+x /opt/required_scripts/annotSeg_7102f1c.pl
ENV JS_PIPELINE_PATH /opt

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module