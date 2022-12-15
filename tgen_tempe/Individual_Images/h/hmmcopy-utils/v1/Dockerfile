FROM ghcr.io/tgen/jetstream_containers/r:3.6.1 

# Installing requirements for HMMCopy_utils
RUN dnf upgrade -y && dnf install -y \
    wget  \
    unzip  \
    cmake  \
    gcc-c++  \
    zlib-devel &&\
    dnf clean all

# Download and install HMMCopy_utils
RUN wget https://github.com/shahcompbio/hmmcopy_utils/archive/master.zip &&\
    unzip master.zip &&\
    rm master.zip &&\
    cd hmmcopy_utils-master &&\
    cmake . &&\
    make &&\
    mv bin/* /usr/local/bin &&\
    mv util/* /usr/local/bin &&\
    mv lib/* /usr/local/lib &&\
    cd / &&\
    rm -rf hmmcopy_utils-master

# ghcr.io/tgen/jetstream_containers/hmmcopy-utils:v1

