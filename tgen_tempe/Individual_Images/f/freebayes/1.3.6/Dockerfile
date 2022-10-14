FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

# Versions
ARG FREEBAYES_VER=1.3.6

# Installing requirements for freebayes
RUN dnf upgrade -y && dnf install -y \
    wget \
    tar &&\
    dnf clean all 

# Download and install freebayes
RUN wget https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES_VER}/freebayes-${FREEBAYES_VER}-linux-amd64-static.gz -q --show-progress && \
    gzip -d freebayes-${FREEBAYES_VER}-linux-amd64-static.gz &&\
    mv freebayes-${FREEBAYES_VER}-linux-amd64-static /usr/local/bin/freebayes &&\
    chmod +x /usr/local/bin/freebayes &&\
    rm -rf freebayes-${FREEBAYES_VER}-linux-amd64-static.gz

# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/freebayes:1.3.6