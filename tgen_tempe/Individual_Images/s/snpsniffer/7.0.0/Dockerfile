FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/r:3.6.1

# Environment Variables
ENV PATH=$PATH:/snpSniffer

# Versions
ARG SNPSNIFFER_VER=7.0.0

# Copy Bash Script and CSV Module List
COPY module_install.sh /
COPY req_modules.csv /

# Installing pre requisits for snpSniffer
RUN dnf upgrade -y && dnf install -y wget tar gzip java-11-openjdk &&\
    dnf clean all

# Download and install snpSniffer
RUN wget https://github.com/tgen/snpSniffer/archive/v${SNPSNIFFER_VER}.tar.gz -q --show-progress && \
    tar zxf v${SNPSNIFFER_VER}.tar.gz && \
    rm v${SNPSNIFFER_VER}.tar.gz && \
    ln -s snpSniffer-${SNPSNIFFER_VER} snpSniffer &&\
    chmod +x /snpSniffer/snpSniffer.jar &&\
    chmod +x /snpSniffer/snpSniffer_Summarize.*

# Run Install R Modules Script
RUN chmod +x module_install.sh && \
    ./module_install.sh && \
    rm module_install.sh && \
    rm req_modules.csv

CMD [ "snpSniffer" , "--help"]
# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/snpsniffer:7.0.0
