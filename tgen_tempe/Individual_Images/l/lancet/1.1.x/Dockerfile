FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

# Versions
ARG LANCET_VER=1.1.x

# Install requirements for lancet
RUN dnf upgrade -y && dnf install -y wget unzip cmake gcc-c++ make xz-devel libcurl-devel zlib-devel bzip2-devel openssl-devel &&\
    dnf clean all

# Download and install lancet
RUN wget https://github.com/nygenome/lancet/archive/v${LANCET_VER}.zip -q --show-progress && \
    unzip v${LANCET_VER}.zip && \
    rm v${LANCET_VER}.zip && \
    cd lancet-${LANCET_VER} && \
    make && \
    mv lancet /usr/local/bin && \
    cd / && \
    rm -rf lancet-${LANCET_VER}

# Version Checks
RUN echo "LANCET_VER=${LANCET_VER}"

CMD [ "lancet", "-h" ]
#pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/lancet:1.1.x