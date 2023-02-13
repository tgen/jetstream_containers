FROM ghcr.io/tgen/jetstream_containers/python:3.7.2 

ARG POD5_VER=0.1.5
RUN dnf upgrade  -y &&\
    pip install --no-cache-dir --upgrade  pip && \
    pip install --no-cache-dir pod5==${POD5_VER}

# ghcr.io/tgen/jetstream_containers/pod5:0.1.5