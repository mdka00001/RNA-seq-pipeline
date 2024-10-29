# RNA-seq Data Processing Pipeline

This repository contains a Snakemake pipeline for processing RNA-seq data from FASTQ download to read alignment and feature quantification. The workflow is divided into three main stages: **FASTQ Download and Quality Control**, **Processing**, and **Alignment**. Each stage has its own Snakemake file for modular processing, enhancing readability and flexibility.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Workflow](#pipeline-workflow)
  - [1. FASTQ Download and QC](#1-fastq-download-and-qc)
  - [2. Processing](#2-processing)
  - [3. Alignment](#3-alignment)
- [Output Structure](#output-structure)
- [License](#license)

## Overview

This RNA-seq processing pipeline covers:
1. **FASTQ Download and QC** - Downloads FASTQ files, performs quality control, and generates a MultiQC report.
2. **Processing** - Adapts trimming, quality checks the trimmed reads, and compiles QC reports.
3. **Alignment** - Aligns reads to a reference genome and performs feature counting.

## Installation

The pipeline relies on Snakemake, Conda, and specific bioinformatics tools. To set up the environment, you can use the provided `environment.yaml` file. Ensure you have Conda installed.

### Steps
1. Clone this repository:
    ```bash
    git clone https://github.com/yourusername/your-repository-name.git
    cd your-repository-name
    ```
2. Install the required conda environment:
    ```bash
    conda env create -f environment.yaml
    conda activate rna-seq-pipeline
    ```

## Configuration

The pipeline uses a `config.yaml` file to specify parameters like output directory, adapter sequences, accessions, reference genome, and other necessary inputs. Update `config.yaml` with paths and variables specific to your dataset and analysis requirements.

### Example `config.yaml`
```
accessions:
  - ERR10949268
  - ERR10949270
  - ERR10949271
  - ERR10949272
  - ERR10949274
  - ERR10949273
  - ERR10949269
  - ERR10949275
adapter: /home/sysgen/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
refseq: /data/karim_thesis/NGS_test/grch38_tran/genome_tran
index: /data/karim_thesis/AD_organoids/Homo_sapiens.GRCh38.112.gtf
output_dir: /data/karim_thesis/AD_organoids/forebrain_organoids_2/
```


