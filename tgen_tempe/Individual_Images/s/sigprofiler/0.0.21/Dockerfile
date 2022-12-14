FROM python:3.8.16-slim

ENV PATH=$PATH:/sigprofiler

# wget is needed to pull signatures, genInstall.py will use rsync in specific environments
RUN apt update && \
    apt install -y wget rsync

# Install SigProfilerAssignment and dependencies with locked versions
COPY genInstall.py .
COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install -r requirements.txt && rm requirements.txt && \
    python genInstall.py && rm genInstall.py

# Copy in sigprofiler_assignment script
WORKDIR /sigprofiler
COPY sigprofiler_assignment.py .

# Pull GRCh38 signatures for sigprofiler
WORKDIR /sigprofiler_resources
RUN wget https://cancer.sanger.ac.uk/signatures/documents/2047/COSMIC_v3.3.1_SBS_GRCh38.txt && \
    wget https://cancer.sanger.ac.uk/signatures/documents/1903/COSMIC_v3.3_DBS_GRCh38.txt

WORKDIR /

# ghcr.io/tgen/jetstream_containers/sigprofiler:0.0.21