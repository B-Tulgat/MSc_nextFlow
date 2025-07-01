FROM mambaorg/micromamba:1.5.8

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
        wget \
        curl \
        git \
        build-essential \
        && rm -rf /var/lib/apt/lists/*

RUN apt install -y openjdk-17-jdk
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x nextflow
RUN mv nextflow /usr/local/bin/

USER root

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.2-linux-x64.tar.gz -O /tmp/dorado.tar.gz && \
mkdir -p /opt && \
tar -xzf /tmp/dorado.tar.gz -C /opt && \
mv /opt/dorado-1.0.2-linux-x64 /opt/dorado && \
ln -s /opt/dorado/bin/dorado /usr/local/bin/dorado && \
rm /tmp/dorado.tar.gz


USER $MAMBA_USER

# Create environment and install all tools in one step
RUN micromamba create -n nanopore -c conda-forge -c bioconda -y \
python=3.9 \
samtools=1.20 \
minimap2=2.17 \
pandas \
numpy \
scipy \
pysam \
ont-modkit && \
micromamba clean --all --yes

# Activate environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /workspace

RUN micromamba activate nanopore

# Test installations
RUN dorado --version && \
modkit --version && \
samtools --version && \
minimap2 --version
