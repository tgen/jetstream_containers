FROM pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/rocky8_with_epel:2022.04.19

# Versions
ARG PYTHON_VER=3.7.2

# Install requirements for Python 3.7.2
RUN dnf upgrade -y && dnf install -y \
    wget \
    make \
    gcc \
    tar \
    openssl-devel \
    bzip2-devel \
    libffi-devel \
    zlib-devel \
    xz-devel \
    gcc-c++ &&\
    dnf clean all

# download and install Python 3.7.2
WORKDIR /usr/src 
RUN wget https://www.python.org/ftp/python/${PYTHON_VER}/Python-${PYTHON_VER}.tgz &&\
    tar xzf Python-${PYTHON_VER}.tgz &&\
    cd Python-${PYTHON_VER} &&\
    ./configure --enable-optimizations &&\
    make altinstall &&\
    cd .. &&\
    rm -rf Python-${PYTHON_VER}.tgz

RUN ln -sfn /usr/local/bin/python3.7 /usr/bin/python3.7 &&\
    ln -sfn /usr/local/bin/python3.7 /usr/bin/python3 &&\
    ln -sfn /usr/local/bin/python3.7 /usr/bin/python &&\
    ln -sfn /usr/local/bin/pip3.7 /usr/bin/pip3.7 &&\
    pip3.7 install --upgrade pip &&\
    python3.7 -m pip install pysam

# pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/python:3.7.2