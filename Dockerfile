# Use an official Ubuntu base image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/bioinformatics_tools/FastQC:/bioinformatics_tools/Trimmomatic-0.39:/bioinformatics_tools/STAR-2.7.11b/bin/Linux_x86_64_static:/bioinformatics_tools/samtools-1.19.2:${PATH}"

# Update the package list and install required dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    build-essential \
    openjdk-8-jdk \
    python3.11 \
    python3-pip \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a directory for bioinformatics tools
RUN mkdir -p /bioinformatics_tools
WORKDIR /bioinformatics_tools

# Install FastQC v0.12.1
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc

# Install Trimmomatic v0.39
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip

# Install STAR v2.7.11b
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz && \
    tar -xvzf 2.7.11b.tar.gz && \
    cd STAR-2.7.11b/source && make STAR && cd ../../

# Install Samtools v1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.2.tar.bz2 && \
    tar -xjf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && make && cd ..

# Install Picard v3.1.1
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar

# Install featureCounts (subread) v2.0.6
RUN wget https://downloads.sourceforge.net/project/subread/subread-2.0.6-Linux-x86_64.tar.gz && \
    tar -xvzf subread-2.0.6-Linux-x86_64.tar.gz

# Install Python dependencies
RUN pip3 install numpy==1.26.4 matplotlib==3.8.3 pysam==0.22.0 deeptools==3.5.3

# Install R dependencies
RUN Rscript -e 'install.packages(c("BiocManager", "pheatmap", "ggplot2", "dplyr", "igraph", "biomaRt"), repos="http://cran.us.r-project.org")' && \
    Rscript -e 'BiocManager::install(c("DESeq2", "vsn", "BiocGenerics"))'

# Define an entrypoint for the container
ENTRYPOINT ["/bin/bash"]

