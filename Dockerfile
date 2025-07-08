FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        curl \
        git \
        build-essential \
        openjdk-17-jdk \
        && rm -rf /var/lib/apt/lists/*

# Install Nextflow in a single layer, ensuring it's executable and in PATH
RUN curl -sL https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mv nextflow /usr/local/bin/

# Install Dorado in a single layer, cleaning up tarball
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.2-linux-x64.tar.gz -O /tmp/dorado.tar.gz && \
    mkdir -p /opt/dorado && \
    tar -xzf /tmp/dorado.tar.gz -C /opt/dorado --strip-components=1 && \
    ln -s /opt/dorado/bin/dorado /usr/local/bin/dorado && \
    rm /tmp/dorado.tar.gz

USER $MAMBA_USER

# Now, micromamba create will use /opt/micromamba as its root
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

ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Explicitly tell micromamba which env to activate by default
ENV ENV_NAME=nanopore

# Test installations as the non-root user
RUN micromamba run -n nanopore dorado --version && \
    micromamba run -n nanopore modkit --version && \
    micromamba run -n nanopore samtools --version && \
    micromamba run -n nanopore minimap2 --version
