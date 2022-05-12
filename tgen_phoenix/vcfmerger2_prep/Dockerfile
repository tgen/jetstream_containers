ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=skylake

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-3.9.7.headless:${IMAGE_ARCH} AS opt_python3
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_python-2.7.18.headless:${IMAGE_ARCH} AS opt_python2
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_curl-7.78.0:${IMAGE_ARCH} AS opt_curl
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_grep-3.7:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gawk-5.1.0:${IMAGE_ARCH} AS opt_gawk
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_gzip-1.11:${IMAGE_ARCH} AS opt_gzip
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bcftools-1.13:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_util-linux-2.37.2:${IMAGE_ARCH} AS opt_util-linux
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_vt-0.57721:${IMAGE_ARCH} AS opt_vt
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_samtools-1.13:${IMAGE_ARCH} AS opt_samtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_bedtools-2.29.2:${IMAGE_ARCH} AS opt_bedtools
FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_htslib-1.13:${IMAGE_ARCH} AS opt_htslib

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build_1

COPY --from=opt_python3 /opt/ /opt/
COPY --from=opt_curl /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /opt && \
    cd /opt && git clone https://github.com/tgen/vcfMerger2.git && \
    cd /opt/vcfMerger2 && git checkout v0.8.7

RUN python -m pip install cyvcf2 intervene Pillow natsort

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build_2

COPY --from=opt_python2 /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN /opt/bin/python2 -m pip install cython

RUN mkdir -p /tmp/scratch && \
    cd /tmp/scratch && git clone https://github.com/secastel/phaser.git && \
    #cd /tmp/scratch/phaser/phaser/ && /opt/bin/python2 -m pip install ./ && \
    cd /tmp/scratch/phaser/phaser/ && python2 setup.py build_ext --inplace && \
    #mkdir -p /opt/lib/python2.7/site-packages/phaser && touch /opt/lib/python2.7/site-packages/phaser/__init__.py && \
    #cd /tmp/scratch/phaser/phaser/ && cp read_variant_map.so /opt/lib/python2.7/site-packages/phaser/ && \
    cd /tmp/scratch/phaser/phaser/ && cp phaser.py call_read_variant_map.py read_variant_map.py read_variant_map.so /opt/bin/ && \
    cd /opt/bin && chmod a+x phaser.py call_read_variant_map.py read_variant_map.py && \
    sed -i 's|/home/bin/|/opt/bin/|g' /opt/bin/phaser.py && \
    rm /opt/bin/python && \
    cd / && rm -rf /tmp/scratch

RUN /opt/bin/python2 -m pip install scipy pysam

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_tools:${IMAGE_ARCH} as build_3

ARG IMAGE_ARCH

RUN mkdir -p /tmp/scratch 
COPY xxd.c /tmp/scratch
RUN cd /tmp/scratch && mkdir -p /opt/bin &&  gcc -O3 -march=${IMAGE_ARCH} xxd.c -o /opt/bin/xxd

FROM ${SOURCE_CONTAINER_REGISTRY}/leanix-1.1_opt_base:${IMAGE_ARCH}

COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=opt_gzip /opt/ /opt/
COPY --from=opt_bcftools /opt/ /opt/
COPY --from=opt_util-linux /opt/ /opt/
COPY --from=opt_vt /opt/ /opt/
COPY --from=opt_samtools /opt/ /opt/
COPY --from=opt_bedtools /opt/ /opt/
COPY --from=build_1 /opt/ /opt/
COPY --from=build_2 /opt/ /opt/
COPY --from=build_3 /opt/ /opt/
COPY --from=opt_htslib /opt/ /opt/

ENV PATH /opt/vcfMerger2/bin:${PATH}

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module

