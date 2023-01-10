FROM ghcr.io/tgen/jetstream_containers/r:3.6.1

# Versions
ARG VARDICT_VER=1.7.0

# Install requirements for BCFtool
RUN dnf upgrade -y && dnf install -y java-11-openjdk perl &&\
    dnf clean all

# Download and install VarDict
RUN wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/${VARDICT_VER}/VarDict-${VARDICT_VER}.tar -q --show-progress && \
    tar xf VarDict-${VARDICT_VER}.tar && \
    rm VarDict-${VARDICT_VER}.tar && \
    cd VarDict-${VARDICT_VER} && \
    cp -R bin/* /usr/local/bin && \
    cp -R lib/* /usr/local/lib && \
    cd / && \
    rm -rf VarDict-${VARDICT_VER}

CMD [ "VarDict", "-H" ]
# ghcr.io/tgen/jetstream_containers/vardict:1.7.9
