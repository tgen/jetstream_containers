FROM starfusion:1.11.0

ARG BWA_MEM2_VER=2.2.1

# # Install BWA-Mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VER}/bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 &&\
    tar jxf bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 &&\
    rm -rf bwa-mem2-${BWA_MEM2_VER}_x64-linux.tar.bz2 &&\
    mv bwa-mem2-${BWA_MEM2_VER}_x64-linux bwa-mem2 &&\
    mv bwa-mem2/bwa* /usr/local/bin/

# ghcr.io/tgen/jetstream_containers/trinity-bwa-samtools:2.2.1