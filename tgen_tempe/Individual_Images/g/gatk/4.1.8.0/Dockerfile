FROM ghcr.io/tgen/jetstream_containers/r:3.6.1

# Environment Variables
ENV PATH=$PATH:/gatk
ENV PATH=$PATH:/miniconda3/bin/
ENV CONDA_DIR /miniconda3

# Versions
ARG GATK_VER=4.1.8.0

# Install requirements for GATK
RUN dnf upgrade -y && dnf install -y \
    wget \
    unzip \
# Java 8
    java-1.8.0-openjdk \ 
    java-1.8.0-openjdk-devel \
    python38 \
# Install Mini-Conda Dependencies
    glibc &&\
    dnf clean all 

# Change default Python to 3.8
RUN alternatives --set python /usr/bin/python3.8

# Download and install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip -q --show-progress &&\
    unzip gatk-${GATK_VER}.zip &&\
    rm -rf gatk-${GATK_VER}.zip &&\
    mv gatk-${GATK_VER} gatk

# Install Mini-Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh &&\
    chmod -R a+rw /gatk &&\
    conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -afy && \
    find /miniconda3/ -follow -type f -name '*.a' -delete && \
    find /miniconda3/ -follow -type f -name '*.pyc' -delete && \
    rm -rf /root/.cache/pip

# Downgraded h5py to version 2.10.0 to resolve issue
#https://gatk.broadinstitute.org/hc/en-us/articles/360052489832-Known-Issue-with-CNNScoreVariants-version-4-1-9-0
RUN conda run -n gatk conda install -y h5py=2.10.0

# Version Check
RUN echo "GATK_VER=${GATK_VER}"

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]
# ghcr.io/tgen/jetstream_containers/gatk:4.1.8.0