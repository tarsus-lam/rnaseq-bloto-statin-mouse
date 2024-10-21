# Role of ZNF335 in Cellular Cholesterol Metabolism and Statin Response

## Project Overview

This project investigates the role of **ZNF335** in regulating cellular cholesterol metabolism and statin response. Liver samples from CRISPRi/CRISPRa-modified mice were analyzed using RNA sequencing (RNA-seq) to identify differentially expressed genes and key pathways involved in lipid metabolism, focusing on the effects of Zfp335 and its potential impact on high cholesterol and statin efficacy.

## Objectives

1. Analyze the role of Zfp335 in cholesterol metabolism using CRISPRi/CRISPRa-modified mice.
2. Identify differentially expressed genes involved in cholesterol and lipid metabolism.
3. Conduct functional enrichment analysis to explore regulatory networks related to Zfp335.
4. Assess Zfp335’s potential impact on statin response.

## Installation Methods

You can install the required dependencies in two ways: using a bash script or a Docker container.

### Method 1: Install Dependencies via Bash Script

Download and run the install_tools.sh script:

```bash
bash install_tools.sh
```

### Method 2: Install Dependencies via Docker

Alternatively, build and run a Docker container with all dependencies:

1. Build the Docker Image:
```bash
docker build -t rnaseq_tools .
```
2. Run the Docker Container:
```bash
docker run -it --rm -v /path/to/your/data:/data rnaseq_tools
```

## Workflow and Key Scripts

### 1. Quality Control (QC)
The function of this step is to assess the quality of raw sequencing reads to identify any potential issues with the data before further analysis.

- **Input**: Raw FASTQ files are used as input to check for quality.

- **Output**: FastQC quality reports are generated, detailing per-base quality scores and other key metrics.

- **Script**: scriptFastQC

### 2. Trimming and Filtering
This step removes low-quality bases and adapter sequences from the raw reads, ensuring high-quality data for downstream analysis.

- **Input**: Raw FASTQ files are provided as input to be trimmed and filtered for quality.
  
- **Output**: The output consists of trimmed FASTQ files, with adapters and low-quality sequences removed.
  
- **Script**: scriptTrimmomatic

### 3. Genome Indexing
This step generates a genome index to prepare the reference genome for efficient read alignment in later steps.

- **Input**: The input includes the mouse genome FASTA file and GTF annotation file.

- **Output**: STAR genome index files are created, which are used for aligning RNA-seq reads.

- **Script**: scriptSTARIndex

### 4. First Pass Alignment
This step aligns the RNA-seq reads to the mouse genome to generate preliminary BAM files and identify splice junctions.

- **Input**: Trimmed FASTQ files are used as input to perform the first pass of read alignment.

- **Output**: The output is BAM files containing aligned reads, along with identified splice junctions.

- **Script**: scriptSTARPass1

### 5. Second Pass Alignment
This step refines the alignment process by using the splice junctions identified in the first pass to improve the accuracy of the alignments.

- **Input**: The input includes BAM files from the first pass, along with splice junction data.

- **Output**: Final aligned BAM files are produced, with improved read alignment around splice sites.

- **Script**: scriptSTARPass2

### 6. Mappability Statistics
This step evaluates the mapping quality of the RNA-seq data by calculating percentages of mapped, unmapped, and multi-mapped reads.

- **Input**: STAR log files from the alignment process are provided as input to analyze mappability statistics.

- **Output**: Mappability statistics are generated, along with a pie chart showing the proportion of mapped, unmapped, and multi-mapped reads.

- **Script**: alignment_analysis.py

### 7. BAM File Processing (Indexing, Marking Duplicates, Merging)
This step involves processing BAM files by indexing them, marking duplicate reads, and merging multiple BAM files into combined datasets.

- **Input**: The input includes the BAM files generated from the second pass alignment.

- **Output**: Indexed, duplicate-marked, and merged BAM files are produced, along with duplicate metrics plots for quality assessment.

- **Script**s: scriptIndexBam, scriptMarkDuplicates, scriptMergeBam, duplicate_metrics_plot.py

### 8. Read Coverage
This step calculates and visualizes the read coverage across the genome to assess the evenness and depth of sequencing coverage.

- **Input**: BAM files from the alignment process are used as input to calculate read coverage.

- **Output**: A read coverage plot is generated, and raw counts of coverage across the genome are output for further analysis.

- **Script**s: scriptPlotCoverage, calculate_coverage.py

### 9. Feature Counting
This step quantifies the number of reads mapped to each gene, providing a gene-level count matrix for differential expression analysis.

- **Input**: BAM files are provided as input to count the number of reads mapped to genes.

- **Output**: A gene count matrix is produced, detailing the number of reads aligned to each gene.

- **Script**: scriptFeatureCounts

### 10. Differential Expression Analysis
This step identifies genes that are differentially expressed between experimental conditions, such as wild-type vs. CRISPR-modified samples.

- **Input**: The input includes the gene count matrix and sample information.

- **Output**: Differentially expressed genes are identified, and visualizations such as heatmaps and PCA plots are generated to explore the data.

- **Script**s: scriptDESeq2, scriptDESeq2.R

### 11. Functional Enrichment
This step performs functional enrichment analysis to identify biological pathways and functions that are significantly associated with the differentially expressed genes.

- **Input**: A list of differentially expressed genes is used as input for functional enrichment.

- **Output**: The output includes a list of enriched biological functions or pathways based on the gene list.

- **Script**s: scriptEnrichR, scriptEnrichR.R

## Running the Workflow

### Option 1: Running the Entire Workflow with Nextflow
Once you have Nextflow installed and configured to use Docker, you can execute the entire workflow with:

```bash
nextflow run main.nf -with-docker
```

This command will run all steps of the pipeline automatically, leveraging the Docker container for a consistent environment across all processes.

### Option 2: Running Individual Scripts on an SGE Cluster
You can also manually submit individual scripts to an SGE-managed cluster. For example, to submit the **FastQC** step, use the following:

```bash
qsub scriptFastQC.sh
```

Repeat this for each step in the pipeline as required by your cluster environment.

## Conclusion

By following the steps outlined in this README, you will be able to replicate the RNA-seq analysis workflow to study the role of ZNF335 in cholesterol metabolism and statin response. The workflow includes quality control, alignment, differential expression analysis, and functional enrichment.
