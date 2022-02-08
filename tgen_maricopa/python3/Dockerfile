FROM rockybase as build

RUN dnf upgrade -y
RUN dnf install -y python38 curl gcc libcurl-devel python38-devel zlib-devel

RUN alternatives --set python /usr/bin/python3.8
RUN alternatives --install /usr/bin/pip pip /usr/bin/pip3.8 1
RUN pip install requests multiqc requests deeptools numpy scipy py2bit pyBigWig pysam matplotlib MACS2

FROM rockybase

# copy build python modules from build stage
COPY --from=build /usr/local/* /usr/local

RUN dnf upgrade -y
RUN dnf install -y python38 zlib libcurl

RUN alternatives --set python /usr/bin/python3.8
RUN alternatives --install /usr/bin/pip pip /usr/bin/pip3.8 1

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
