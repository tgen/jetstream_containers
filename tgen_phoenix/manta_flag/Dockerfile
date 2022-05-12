ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} as opt_python
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} as opt_bcftools

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build

COPY --from=opt_python /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN python -m pip install pysam

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=build /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/

COPY manta_prepare_sv_vcf_f94bcc1.py /opt/phoenix/required_scripts/
ENV JS_PIPELINE_PATH /opt/phoenix/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

