ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_r-4.0.2.headless:${IMAGE_ARCH} AS opt_r
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_grep-3.4:${IMAGE_ARCH} AS opt_grep
FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_gawk-5.0.1:${IMAGE_ARCH} AS opt_gawk

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_tools:${IMAGE_ARCH} AS build

ARG IMAGE_ARCH
ARG LOCAL_SOURCES_SERVER

COPY --from=opt_r /opt /opt

ENV PATH /opt/bin/:${PATH}
ENV PKG_CONFIG_PATH /opt/pkgconfig/:/opt/lib64/:/opt/lib64/pkgconfig/:/opt/lib/:/opt/lib/pkgconfig/:${PKG_CONFIG_PATH}
ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

RUN mkdir -p /tmp && \
    mkdir -p ~/.R/ && \
    echo "CFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "CXXFLAGS=-O3 -march=${IMAGE_ARCH}" >> ~/.R/Makevars && \
    echo "options(Ncpus = 8)" > ~/.Rprofile && \
    echo "install.packages(\"optparse\", repos=\"http://cran.rstudio.com\")" | R --no-save && \
    echo "install.packages(\"data.table\", repos=\"http://cran.rstudio.com\")" | R --no-save && \
    cd / && rm -rf /tmp/*

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_grep /opt/ /opt/
COPY --from=opt_gawk /opt/ /opt/
COPY --from=build /opt/ /opt/

ENV R_HOME /opt

ENV PATH_TO_PHOENIX_REPO /opt/phoenix/
COPY plotCNVplus_06b34ff.R /opt/phoenix/required_scripts/

RUN mkdir -p /opt/bin/ && echo "#!/bin/bash" > /opt/bin/module && chmod a+x /opt/bin/module
LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
