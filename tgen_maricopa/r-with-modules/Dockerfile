FROM r-3.6.1

# Installing additional R Modules
# Created a bash script to crate this section from a csv file.
# Created the initial list with a bash script and extracting the list
# from the EasyBuild easyconfig file.
# Using script and CSV to install modules.  Adding list to %post was to long.  Build failed
# due to size of definition file.  This solution is also a lot easier to read and maintain.
#
# Copy csv file with list of modules and versions and script to install from csv file.
COPY R_modules.csv /R_modules.csv
COPY module_install.sh /module_install.sh

RUN apt install -y sed
RUN chmod +x module_install.sh
RUN ./module_install.sh
RUN rm module_install.sh
RUN rm R_modules.csv
RUN sed -i "s/options(show.error.messages = TRUE)/options(bitmapType='cairo')/g" /usr/local/lib/R/library/base/R/Rprofile

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
