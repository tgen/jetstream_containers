FROM debian10base

RUN apt update && apt full-upgrade -y
RUN apt -y install wget sed locales tar bzip2 g++ cmake libhts-dev libcurl4-gnutls-dev libbz2-dev libssl-dev curl x11-common pkg-config libgl1-mesa-dev libglu1-mesa-dev libcairo2-dev libreadline-dev libpcre++-dev libpng++-dev libjpeg-dev libtiff-dev openjdk-11-jdk tk8.6-dev libxml2-dev libgmp-dev libnlopt-dev libfftw3-dev libsndfile1-dev libhdf5-dev libudunits2-dev libgsl-dev libmagick++-dev gfortran libopenmpi-dev

# Uncomment en_US.UTF-8 for inclusion in generation to avoid perl warnings regarding locale generation.
RUN sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen

# Generate locale
RUN locale-gen

RUN wget https://cloud.r-project.org/src/base/R-3/R-3.6.1.tar.gz
RUN tar zxf R-3.6.1.tar.gz
RUN rm R-3.6.1.tar.gz
RUN cd R-3.6.1 && ./configure --enable-R-shlib && make && make install
RUN cd /
RUN rm -rf R-3.6.1

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
