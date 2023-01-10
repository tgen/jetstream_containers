FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

ENV PATH=$PATH:/tHReD-1.1.0

ARG THRED_VER=1.1.0
COPY tHReD-1.1.0.tar.gz /

RUN dnf upgrade -y && dnf install -y \
    tar \
    python38 \
    python38-devel &&\
    dnf clean all

# Change default Python to 3.8 to create /usr/bin/python
RUN alternatives --set python /usr/bin/python3 &&\
    alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

RUN python -m pip install natsort pandas matplotlib

RUN tar zxf tHReD-1.1.0.tar.gz &&\
    rm -rf tHReD-1.1.0.tar.gz &&\
    chmod a+x /tHReD-1.1.0/tHReD.py

# ghcr.io/tgen/jetstream_containers/thred:1.1.0 