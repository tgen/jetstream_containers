FROM ghcr.io/tgen/jetstream_containers/bcftools:1.10.2

# Environment Variables
ENV PATH=$PATH:/snpEff:/clinEff:/snpEff/scripts

# Versions
ARG SNPEFF_VER=4_3t

# Install additional requirements for snpEFF
RUN dnf upgrade -y && dnf install -y \
    wget \
    gzip \
    unzip \
    java-1.8.0-openjdk \ 
    java-1.8.0-openjdk-devel &&\
    dnf clean all 

# Download and install snpEff
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v${SNPEFF_VER}_core.zip/download -q --show-progress &&\
    mv download snpEff_v${SNPEFF_VER}_core.zip &&\
    unzip snpEff_v${SNPEFF_VER}_core.zip &&\
    rm snpEff_v${SNPEFF_VER}_core.zip

WORKDIR /snpEff

CMD [ "java","-jar","snpEff.jar","-h" ]

# ghcr.io/tgen/jetstream_containers/snpeff:4_3t
