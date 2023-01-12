FROM ghcr.io/tgen/jetstream_containers/r-base:3.6.1

ARG BCFTOOLS_VER=1.10.2
ARG SAMTOOLS_VER=1.10
ARG BAMTOOLS_VER=2.5.1
ARG BEDTOOLS_VER=2.29.0
ARG VT_VER=0.57721
ARG SNPEFF_VER=4_3t
ARG VCFMERGER_VER=0.9.1

# Environment Variables
ENV LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH
ENV PATH=/vcfMerger2/bin:/phaser/phaser:/snpEff:/clinEff:$PATH

RUN dnf upgrade -y && \

### BCFTools ###
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2 -q --show-progress &&\
    tar jxf bcftools-${BCFTOOLS_VER}.tar.bz2 &&\
    cd bcftools-${BCFTOOLS_VER} && ./configure --enable-libgsl && make && make install && cd / &&\
    rm bcftools-${BCFTOOLS_VER}.tar.bz2 &&\
    rm -rf bcftools-${BCFTOOLS_VER} && \
#
### SAMTOOLS ###
# Download and install SAMtools
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && make && make install && \
    cd / && rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
    rm -rf samtools-${SAMTOOLS_VER} && \

### BEDTOOLS ###
#Install dependencies BamTools

dnf install -y json-glib-devel \
#Installing dependencies for BEDTools
               python38 && \
    alternatives --set python /usr/bin/python3.8 && \
    alternatives --install /usr/bin/pip pip /usr/bin/pip3.8 1 && \

# Download and install BamTools
    wget https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VER}.tar.gz -q --show-progress &&\
    tar -zxf v${BAMTOOLS_VER}.tar.gz && \
    cd bamtools-${BAMTOOLS_VER} && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && make && make install && \
    cd / && rm v${BAMTOOLS_VER}.tar.gz && \
    rm -rf bamtools-${BAMTOOLS_VER} && \

# Download and install BedTools
    wget https://github.com/arq5x/bedtools2/archive/v${BEDTOOLS_VER}.tar.gz -q --show-progress &&\
    tar -zxf v${BEDTOOLS_VER}.tar.gz && \
    cd bedtools2-${BEDTOOLS_VER} && \
    make && \
    cp /bedtools2-${BEDTOOLS_VER}/bin/* /usr/local/bin && \
    cd / && rm v${BEDTOOLS_VER}.tar.gz && \
    rm -rf bedtools2-${BEDTOOLS_VER} && \

###  VT ###

# Download and install vt
    wget https://github.com/atks/vt/archive/refs/tags/${VT_VER}.tar.gz -q --show-progress && \
    tar zxf ${VT_VER}.tar.gz && \
    cd vt-${VT_VER} && \
    make && \
    mv vt /usr/local/bin && \
    cd / && rm ${VT_VER}.tar.gz && \
    rm -rf vt-${VT_VER} && \

### SNPEFF ###

# Install Dep snpEFF

# Install additional requirements for snpEFF
    dnf install -y unzip && \

# Download and install snpEff
    wget https://sourceforge.net/projects/snpeff/files/snpEff_v${SNPEFF_VER}_core.zip/download -q --show-progress && \
    mv download snpEff_v${SNPEFF_VER}_core.zip && \
    unzip snpEff_v${SNPEFF_VER}_core.zip && \
    rm snpEff_v${SNPEFF_VER}_core.zip && \
# java -jar /snpEff/snpEff.jar

### PHASER ###
# Install additonal requirements for phASER
    dnf install -y git \
                   python2-Cython \
                   python2-devel && \
    pip2 install scipy pysam && \

# Download and install phASER
    git clone https://github.com/secastel/phaser.git && \
    cd /phaser/phaser && \
    python2 setup.py build_ext --inplace && \
    chmod +x phaser.py && \
    sed -i 's/\#\!\/home\/bin\/python2/\#\!\/usr\/bin\/env python2/' phaser.py && \
    cd / && \
# python2.7 /phaser/phaser/phaser.py

### VCFMERGER ###
# Install Dependencies vcfMerger2
    dnf install -y vim-common \
                   python3-Cython \
                   gnu-free-sans-fonts \
                   python38-devel && \

# Install latest VER of pip, required for cyvcf2 install
    pip install --upgrade pip && \
# Pip install additional python requirements for vcfMerger2
    pip install cyvcf2 argparse==1.1 natsort==5.0.3 intervene==0.6.4  && \
# Download and install vcfMerger2
    wget https://github.com/tgen/vcfMerger2/archive/refs/tags/v${VCFMERGER_VER}.tar.gz -q --show-progress && \
    tar zxf v${VCFMERGER_VER}.tar.gz && \
    mv vcfMerger2-${VCFMERGER_VER} vcfMerger2 && \
    rm v${VCFMERGER_VER}.tar.gz && \

# Clean up dnf cache.
    dnf clean all    

CMD [ "vcfMerger2.py", "--help" ]
# ghcr.io/tgen/jetstream_containers/vcfmerger2:0.9.0
