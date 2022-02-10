FROM python3 as build

RUN dnf upgrade -y

RUN dnf install -y gcc python38-devel zlib-devel
# set matplotlib==3.4.3 to support SigProfiler 1.2.0
RUN pip install requests multiqc requests deeptools numpy scipy py2bit pyBigWig pysam matplotlib==3.4.3 MACS2

FROM python3

# copy from build stage
COPY --from=build /usr/local/* /usr/local

RUN dnf upgrade -y
RUN dnf install -y zlib

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
