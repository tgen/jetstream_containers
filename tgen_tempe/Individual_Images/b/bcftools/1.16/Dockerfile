FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19 as build

ENV LD_LIBRARY_PATH /usr/local/lib/:$LD_LIBRARY_PATH

ARG bcftoosVer=1.16

RUN dnf upgrade -y
RUN dnf install -y wget tar bzip2 gcc make zlib-devel gsl-devel perl-libs bzip2-devel xz-devel libcurl-devel openssl-devel

RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoosVer}/bcftools-${bcftoosVer}.tar.bz2
RUN tar jxf bcftools-${bcftoosVer}.tar.bz2
RUN cd bcftools-${bcftoosVer} && ./configure --enable-libgsl && make && make install && cd /
RUN rm bcftools-${bcftoosVer}.tar.bz2
RUN rm -rf bcftools-${bcftoosVer}

FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ENV LD_LIBRARY_PATH /usr/local/lib/:$LD_LIBRARY_PATH

# copy from build stage
COPY --from=build /usr/local/ /usr/local

RUN dnf upgrade -y
RUN dnf install -y zlib gsl perl-libs bzip2-libs xz-libs libcurl openssl-libs

# ghcr.io/tgen/jetstream_containers/bcftools:1.16 
