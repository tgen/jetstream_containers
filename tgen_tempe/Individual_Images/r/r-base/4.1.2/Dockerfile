FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ARG HTSlib_Version=1.14
ARG NLopt_Version=2.6.1
ARG R_Base_Ver=4
ARG R_Minor_Ver=1.2

ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

RUN dnf upgrade -y && \
    dnf install -y wget \
                   tar \
                   bzip2 \
                   curl \
                   gcc-c++ \
                   make \
                   cmake \
                   zlib-devel \
                   libcurl-devel \
                   xz \
                   openssl-devel \
                   bzip2-devel \
                   xz-devel \
                   java-11-openjdk \
                   java-11-openjdk-devel \
                   gmp-devel \
                   hdf5-devel \
                   gsl-devel \
                   openmpi-devel \
                   xorg-x11-server-common \
                   pkgconf-pkg-config \
                   mesa-libGL-devel \
                   mesa-libGLU-devel \
                   cairo-devel \
                   readline-devel \
                   pcre-devel \
                   libpng-devel \
                   openjpeg2 \
                   libtiff-devel \
                   fftw-devel \
                   libsndfile \
                   udunits2-devel \
                   ImageMagick-c++-devel \
                   gcc-gfortran \
                   libxml2-devel \
                   tk-devel \
                   tcl-devel \
                   libXt-devel \
                   gdal-devel \
                   proj-devel \
                   sqlite-devel \
                   geos-devel \
                   libjpeg-turbo-devel && \
    dnf clean all && \
# Install HTSlib - Required for R.
    wget https://github.com/samtools/htslib/releases/download/${HTSlib_Version}/htslib-${HTSlib_Version}.tar.bz2 && \
    tar jxf htslib-${HTSlib_Version}.tar.bz2 && \
    cd htslib-${HTSlib_Version} && ./configure && make && make install && \
    cd / && \
    rm htslib-${HTSlib_Version}.tar.bz2 && \
    rm -rf htslib-${HTSlib_Version} && \

# Download and install NLopt - Required for R.
    wget https://github.com/stevengj/nlopt/archive/refs/tags/v${NLopt_Version}.tar.gz && \
    tar zxf v${NLopt_Version}.tar.gz && \
    cd nlopt-${NLopt_Version} && mkdir build && cd build && cmake .. && make && make install && \
    cd / && \
    rm v${NLopt_Version}.tar.gz && \
    rm -rf nlopt-${NLopt_Version} && \

# Download and install R
    wget https://cloud.r-project.org/src/base/R-${R_Base_Ver}/R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    tar zxf R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    cd R-${R_Base_Ver}.${R_Minor_Ver} && ./configure --enable-R-shlib --with-tcltk && make && make install && \
    cd / && \
    rm R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    rm -rf R-${R_Base_Ver}.${R_Minor_Ver} && \
    R CMD javareconf
