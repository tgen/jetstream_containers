FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ARG HTSlib_Version=1.14
ARG NLopt_Version=2.6.1
ARG R_Base_Ver=3
ARG R_Minor_Ver=6.1

ENV LD_LIBRARY_PATH="/usr/local/lib64":$LD_LIBRARY_PATH
ENV PKG_CONFIG_PATH="/usr/local/lib64/"/pkgconfig/

RUN dnf upgrade -y && \
    dnf install -y wget \
                   tar \
                   bzip2 \
                   curl \
                   gcc-c++ \
                   make \
                   cmake \
                   python38 \
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
                   xorg-x11-server-devel \
                   libX11-devel \
                   libXt-devel \
                   tcl-devel \
                   tk-devel \
                   libjpeg-turbo-devel && \
    # set Python 3.8 as the default version
    alternatives --set python /usr/bin/python3.8 && \
    alternatives --install /usr/bin/pip pip /usr/bin/pip3.8 1 && \
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

# Add /usr/local/lib64 to ld.so.conf
    echo "/usr/local/lib64" >> /etc/ld.so.conf.d/NLopt.conf && \
    ldconfig && \

# Download and install R
    wget https://cloud.r-project.org/src/base/R-${R_Base_Ver}/R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    tar zxf R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    cd R-${R_Base_Ver}.${R_Minor_Ver} && ./configure --enable-R-shlib && make && make install && \
    cd / && \
    sed -i "s/options(show.error.messages = TRUE)/options(bitmapType='cairo')/g" /usr/local/lib64/R/library/base/R/Rprofile && \
    rm R-${R_Base_Ver}.${R_Minor_Ver}.tar.gz && \
    rm -rf R-${R_Base_Ver}.${R_Minor_Ver}

# ghcr.io/tgen/jetstream_containers/r-base:3.6.1 