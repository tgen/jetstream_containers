FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/bwa_mem2_samtools:2.2.1
# Versions
ARG FGBIO_VER=2.0.2
ENV SRC=/usr/local/src

# Installing requirements for Fgbio
RUN dnf upgrade -y && dnf install -y \
    wget \
    java-1.8.0-openjdk \
    java-1.8.0-openjdk-devel

# Download and install Fgbio
RUN cd $SRC && \
    wget "https://github.com/fulcrumgenomics/fgbio/releases/download/${FGBIO_VER}/fgbio-${FGBIO_VER}.jar"

# Create Wrapper
RUN cd /usr/bin && \
    printf '#!/bin/bash\nexec java -jar /usr/local/src/fgbio-2.0.2.jar "$@"' >> fgbio && \
    chmod +x fgbio

WORKDIR /

# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/fgbio:2.0.2