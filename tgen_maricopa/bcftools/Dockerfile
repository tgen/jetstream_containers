FROM rockybase as build

ENV LD_LIBRARY_PATH /usr/local/lib/:$LD_LIBRARY_PATH

RUN dnf upgrade -y
RUN dnf install -y wget tar bzip2 gcc make zlib-devel gsl-devel perl-libs bzip2-devel xz-devel libcurl-devel openssl-devel

RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
RUN tar jxf bcftools-1.10.2.tar.bz2
RUN cd bcftools-1.10.2 && ./configure --enable-libgsl && make && make install && cd /
RUN rm bcftools-1.10.2.tar.bz2
RUN rm -rf bcftools-1.10.2

FROM rockybase

# copy from build stage
COPY --from=build /usr/local/* /usr/local

RUN dnf upgrade -y
RUN dnf install -y zlib gsl perl-libs bzip2-libs xz-libs libcurl openssl-libs

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
