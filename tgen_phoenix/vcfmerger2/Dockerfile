ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge
ARG LOCAL_SOURCES_SERVER=http://localhost:8000

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-3.8.6.headless:${IMAGE_ARCH} AS opt_python3
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-2.7.18.headless:${IMAGE_ARCH} AS opt_python2
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_curl-7.70.0:${IMAGE_ARCH} AS opt_curl
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_grep-3.4:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gawk-5.0.1:${IMAGE_ARCH} AS opt_gawk
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gzip-1.10:${IMAGE_ARCH} AS opt_gzip
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bcftools-1.10.2:${IMAGE_ARCH} AS opt_bcftools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_util-linux-2.35.1:${IMAGE_ARCH} AS opt_util-linux
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_vt-0.57721:${IMAGE_ARCH} AS opt_vt
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_samtools-1.10:${IMAGE_ARCH} AS opt_samtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_bedtools-2.29.2:${IMAGE_ARCH} AS opt_bedtools
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2:${IMAGE_ARCH} AS opt_r

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build_1

COPY --from=opt_python3 /opt/ /opt/
COPY --from=opt_curl /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /opt && \
    cd /opt && git clone https://github.com/tgen/vcfMerger2.git && \
    cd /opt/vcfMerger2 && git checkout v0.8.7

RUN python -m pip install cyvcf2 intervene Pillow natsort

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build_2

COPY --from=opt_python2 /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN /opt/bin/python2 -m pip install cython

RUN mkdir -p /tmp/scratch && \
    cd /tmp/scratch && git clone https://github.com/secastel/phaser.git && \
    cd /tmp/scratch/phaser/phaser/ && /opt/bin/python -m pip install ./ && \
    cd /tmp/scratch/phaser/phaser/ && cp phaser.py call_read_variant_map.py read_variant_map.py /opt/bin/ && \
    cd /opt/bin && chmod a+x phaser.py call_read_variant_map.py read_variant_map.py && \
    sed -i 's|/home/bin/|/opt/bin/|g' /opt/bin/phaser.py && \
    rm /opt/bin/python && \
    cd / && rm -rf /tmp/scratch

RUN /opt/bin/python2 -m pip install scipy pysam

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build_3

ARG IMAGE_ARCH

RUN mkdir -p /tmp/scratch 
COPY xxd.c /tmp/scratch
RUN cd /tmp/scratch && mkdir -p /opt/bin &&  gcc -O3 -march=${IMAGE_ARCH} xxd.c -o /opt/bin/xxd

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} as build_4

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_r /opt/ /opt/

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/share/pkgconfig/:/opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /tmp && \
    mkdir -p ~/.R/ && \
    echo "CFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "CXXFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "options(Ncpus = 8)" > ~/.Rprofile && \
    echo "install.packages(\"Cairo\", repos=\"http://cran.rstudio.com\", dependencies = TRUE)" | R --no-save

RUN mkdir -p /tmp && \
    mkdir -p ~/.R/ && \
    echo "CFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "CXXFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "options(Ncpus = 8)" > ~/.Rprofile && \
    echo "install.packages(\"UpSetR\", repos=\"http://cran.rstudio.com\", dependencies = TRUE)" | R --no-save

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

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
COPY --from=build_4 /opt/ /opt/

ENV PATH /opt/vcfMerger2/bin:${PATH}
ENV R_HOME /opt/

COPY FreeSansBold.ttf /usr/share/fonts/gnu-free/FreeSansBold.ttf

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
