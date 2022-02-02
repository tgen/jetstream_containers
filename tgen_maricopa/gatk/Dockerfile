FROM r-3.6.1-modules as build

RUN apt install -y wget unzip python3 apt-transport-https gnupg
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.0/gatk-4.1.8.0.zip
RUN unzip gatk-4.1.8.0.zip

FROM r-3.6.1-modules

RUN apt install -y python3
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1
RUN mkdir gatk

ENV CONDA_DIR /miniconda3

ENV PATH="/miniconda3/envs/gatk/bin:/miniconda3/bin:/gatk:$PATH"
# ENTRYPOINT ["conda", "run", "-n", "gatk"]
# CMD [".", "/miniconda3/etc/profile.d/conda.sh"]
# CMD ["conda", "run", "-n", "gatk"]

# copy from build stage
COPY --from=build /gatk-4.1.8.0/* /gatk

# Install OpenJDK 8 from adoptopenjdk.net

# Add adoptopenjdk.net repository
RUN wget  https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public
RUN apt-key add public
RUN rm public
RUN echo "deb https://adoptopenjdk.jfrog.io/adoptopenjdk/deb buster main" | tee /etc/apt/sources.list.d/adoptopenjdk.list
RUN apt update
RUN apt install -y adoptopenjdk-8-hotspot


# Install Miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# install in batch (silent) mode, does not edit PATH or .bashrc or .bash_profile
# -p path
# -f force
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3

# Start GATK Python environment

WORKDIR /gatk
RUN chmod -R a+rw /gatk
RUN conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -afy && \
    find /miniconda3/ -follow -type f -name '*.a' -delete && \
    find /miniconda3/ -follow -type f -name '*.pyc' -delete && \
    rm -rf /root/.cache/pip

# Downgraded h5py to version 2.10.0 to resolve issue
#https://gatk.broadinstitute.org/hc/en-us/articles/360052489832-Known-Issue-with-CNNScoreVariants-version-4-1-9-0
RUN conda run -n gatk conda install -y h5py=2.10.0

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers
