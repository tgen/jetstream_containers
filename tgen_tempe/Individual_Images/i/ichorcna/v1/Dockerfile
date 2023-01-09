FROM ghcr.io/tgen/jetstream_containers/r-with_modules:3.6.1

# Installing requirements for HMMCopy_utils
RUN dnf upgrade -y && dnf install -y \
    wget  \
    unzip  \
    cmake  \
    gcc-c++  \
    harfbuzz-devel \
    fribidi-devel \
    zlib-devel &&\
    dnf clean all

# Download and install HMMCopy_utils
RUN wget https://github.com/shahcompbio/hmmcopy_utils/archive/master.zip &&\
    unzip master.zip &&\
    rm master.zip &&\
    cd hmmcopy_utils-master &&\
    cmake . &&\
    make &&\
    mv bin/* /usr/local/bin &&\
    mv util/* /usr/local/bin &&\
    mv lib/* /usr/local/lib &&\
    cd / &&\
    rm -rf hmmcopy_utils-master

# Configure R ichorCNA Modules 
RUN R -e 'install.packages("devtools",dependencies=TRUE,repos = "http://cran.us.r-project.org"); library(devtools);  install_github("broadinstitute/ichorCNA")' &&\
    echo "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")" | R --no-save

# ghcr.io/tgen/jetstream_containers/ichorcna:v1 