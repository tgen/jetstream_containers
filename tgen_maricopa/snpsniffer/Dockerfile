# Uses R image from Phoenix Pipeline as base.
FROM r-3.6.1-modules as build

RUN apt update && apt full-upgrade
RUN apt install -y wget tar gzip

# Install snpSniffer
RUN wget https://github.com/tgen/snpSniffer/archive/v7.0.0.tar.gz
RUN tar zxf v7.0.0.tar.gz
RUN ln -s snpSniffer-7.0.0 snpSniffer

FROM r-3.6.1-modules

ENV PATH="/snpSniffer:$PATH"

# copy from build stage
COPY --from=build /snpSniffer-7.0.0/* /snpSniffer

# create snpSniffer script
RUN echo "#!/bin/bash" > /usr/local/bin/snpSniffer
RUN echo "" >> /usr/local/bin/snpSniffer
RUN echo "cd /snpSniffer" >> /usr/local/bin/snpSniffer
RUN echo 'java -jar /snpSniffer/snpSniffer.jar $@' >> /usr/local/bin/snpSniffer
RUN chmod +x /usr/local/bin/snpSniffer
RUN chmod +x /snpSniffer/snpSniffer.jar

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
