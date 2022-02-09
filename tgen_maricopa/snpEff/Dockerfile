FROM bcftools as build

RUN dnf upgrade -y
RUN dnf install -y wget unzip java-1.8.0-openjdk

RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download
RUN mv download snpEff_v4_3t_core.zip
RUN unzip snpEff_v4_3t_core.zip

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

FROM bcftools

ENV PATH="/snpEff:/clinEff:/snpEff/scripts:$PATH"

RUN dnf upgrade -y
RUN dnf install -y java-1.8.0-openjdk


COPY --from=build /snpEff /snpEff
COPY --from=build /clinEff /clinEff

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
