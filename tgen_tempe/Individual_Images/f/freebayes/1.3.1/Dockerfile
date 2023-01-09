FROM ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19

# Versions
ARG FREEBAYES_VER=1.3.1

# Installing requirements for freebayes
RUN dnf upgrade -y && dnf install -y \
    wget &&\
    dnf clean all 

# Download and install freebayes
# `git clone --single-branch --branch v1.3.1 --recursive https://github.com/freebayes/freebayes.git`
RUN wget https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES_VER}/freebayes-v${FREEBAYES_VER} -q --show-progress && \
    mv freebayes-v${FREEBAYES_VER} /usr/local/bin/freebayes && \
    chmod +x /usr/local/bin/freebayes

CMD [ "freebayes", "--help" ]
# ghcr.io/tgen/jetstream_containers/freebayes:1.3.1 