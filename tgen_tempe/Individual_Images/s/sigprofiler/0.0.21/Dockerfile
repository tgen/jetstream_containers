FROM python:3.8.16-slim

ENV PATH=$PATH:/sigprofiler

# wget is needed to pull signatures, genInstall.py will use rsync in specific environments
RUN sed -i 's/main/main contrib/g' /etc/apt/sources.list && \
    apt update && \
    apt install -y wget rsync ttf-mscorefonts-installer

# Install SigProfilerAssignment and dependencies with locked versions
COPY genInstall.py .
RUN pip install --upgrade pip && \
    pip install about-time==3.1.1 \
    alive-progress==2.4.1 \
    cycler==0.11.0 \
    grapheme==0.6.0 \
    joblib==1.2.0 \
    kiwisolver==1.4.4 \
    matplotlib==3.4.3 \
    numpy==1.23.5 \
    packaging==22.0 \
    pandas==1.5.2 \
    patsy==0.5.3 \
    pdf2image==1.16.0 \
    Pillow==9.3.0 \
    psutil==5.9.4 \
    PyMuPDF==1.21.0 \
    pyparsing==3.0.9 \
    PyPDF2==2.11.2 \
    python-dateutil==2.8.2 \
    pytz==2022.6 \
    reportlab==3.6.12 \
    scikit-learn==1.2.0 \
    scipy==1.9.3 \
    SigProfilerAssignment==0.0.21 \
    SigProfilerMatrixGenerator==1.2.13 \
    sigProfilerPlotting==1.3.6 \
    six==1.16.0 \
    statsmodels==0.13.5 \
    threadpoolctl==3.1.0 \
    typing_extensions==4.4.0 && \
    python genInstall.py && rm genInstall.py

# Copy in sigprofiler_assignment script
WORKDIR /sigprofiler
COPY sigprofiler_assignment.py .
RUN chmod +x sigprofiler_assignment.py

# Pull GRCh38 signatures for sigprofiler
WORKDIR /sigprofiler_resources
RUN wget https://cancer.sanger.ac.uk/signatures/documents/2047/COSMIC_v3.3.1_SBS_GRCh38.txt && \
    wget https://cancer.sanger.ac.uk/signatures/documents/1903/COSMIC_v3.3_DBS_GRCh38.txt

WORKDIR /

# ghcr.io/tgen/jetstream_containers/sigprofiler:0.0.21 