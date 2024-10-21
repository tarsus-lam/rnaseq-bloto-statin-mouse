#!/bin/bash

# Create a directory for tools
mkdir -p bioinformatics_tools
cd bioinformatics_tools

# Download and install FastQC v0.12.1
echo "Downloading FastQC..."
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod +x FastQC/fastqc

# Download and install Trimmomatic v0.39
echo "Downloading Trimmomatic..."
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

# Download and install STAR v2.7.11b
echo "Downloading STAR..."
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
tar -xvzf 2.7.11b.tar.gz
cd STAR-2.7.11b/source
make STAR
cd ../../

# Download and install Samtools v1.19.2
echo "Downloading Samtools..."
wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.2.tar.bz2
tar -xjf samtools-1.19.2.tar.bz2
cd samtools-1.19.2
make
cd ..

# Download and install Picard v3.1.1
echo "Downloading Picard..."
wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar

# Download and install featureCounts (subread) v2.0.6
echo "Downloading featureCounts..."
wget https://downloads.sourceforge.net/project/subread/subread-2.0.6-Linux-x86_64.tar.gz
tar -xvzf subread-2.0.6-Linux-x86_64.tar.gz

# Install Java v1.8.0_401
echo "Installing Java..."
sudo apt-get install openjdk-8-jdk -y

# Install Python v3.11.5 and libraries
echo "Installing Python and dependencies..."
sudo apt-get install python3.11 python3-pip -y
pip3 install numpy==1.26.4 matplotlib==3.8.3 pysam==0.22.0 deeptools==3.5.3

# Install R v4.3.2 and libraries
echo "Installing R and required packages..."
sudo apt-get install r-base -y
Rscript -e 'install.packages(c("BiocManager", "pheatmap", "ggplot2", "dplyr", "igraph", "biomaRt"), repos="http://cran.us.r-project.org")'
Rscript -e 'BiocManager::install(c("DESeq2", "vsn", "BiocGenerics"))'

echo "All tools installed successfully!"

# Additional commands to ensure executables are available system-wide
echo 'export PATH=$PATH:$(pwd)/FastQC:$(pwd)/Trimmomatic-0.39:$(pwd)/STAR-2.7.11b/bin/Linux_x86_64_static:$(pwd)/samtools-1.19.2' >> ~/.bashrc
source ~/.bashrc

echo "Installation complete."

