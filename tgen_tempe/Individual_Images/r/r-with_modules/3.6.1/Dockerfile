FROM ghcr.io/tgen/jetstream_containers/r:3.6.1

COPY module_install.sh /module_install.sh
COPY R_modules.csv /R_modules.csv
COPY Modules/* /

# Installing additional R Modules
# Created a bash script to crate this section from a csv file.
# Created the initial list with a bash script and extracting the list
# from the EasyBuild easyconfig file.
# Using script and CSV to install modules.  Manually insall  was to long to be practicle.

# Get the module CSV list.  AGain will need updated for external facing.
    wget --no-check-certificat https://gitlab01.tgen.org/container_scripts/phoenix-pipeline-2.0/-/raw/master/qc/bam_qc_constitutional_dna/R_modules.csv && \

RUN dnf update && \
    dnf install sed && \

RUN echo "/usr/lib64/openmpi/lib" > /etc/ld.so.conf.d/openmpi.conf && \
    ldconfig
RUN chmod +x module_install.sh && \
    ./module_install.sh && \
    rm module_install.sh && \
<<<<<<< HEAD
    rm R_modules.csv
RUN sed -i "s/options(show.error.messages = TRUE)/options(bitmapType='cairo')/g" /usr/local/lib/R/library/base/R/Rprofile
=======
    rm R_modules.csv
>>>>>>> JimDev
