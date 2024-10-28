import os
import glob
from snakemake.io import glob_wildcards

configfile : "config.yaml"

# Define the directory containing the raw sequence files

# Use glob_wildcards to match SRA files within subdirectories and extract the SRA identifiers
# Define the final output of the workflow
rule all:
    input:
        expand(config["output_dir"]+"rawQC/"+"{sra}_{frr}_fastqc.{extensions}", sra=config["accessions"], frr=[1,2],extensions=["zip","html"]),
        expand(config["output_dir"]+"rawQC/multiqc_report.html")

# Define the rule to convert SRA files to FASTQ format
rule fastqdump:
    output:
        output1 = config["output_dir"]+"rawSeq/"+"{sra}_1.fastq",
        output2 = config["output_dir"]+"rawSeq/"+"{sra}_2.fastq"
    params:
        acc = "{sra}"
    shell:
        """
        fasterq-dump --split-files --outdir {config[output_dir]}/rawSeq/ {params.acc}
        """
rule fastqc:
    input:
        input1 = config["output_dir"]+"rawSeq/"+"{sra}_{frr}.fastq"
    output:
        zip = config["output_dir"]+"rawQC/"+"{sra}_{frr}_fastqc.zip",
        html = config["output_dir"]+"rawQC/"+"{sra}_{frr}_fastqc.html",
    shell:
        """
        fastqc {input.input1} -o {config[output_dir]}/rawQC/
        """

rule multiQC:
    input:
        expand(config["output_dir"]+"rawQC/"+"{sra}_{frr}_fastqc.{extensions}", sra=config["accessions"], frr=[1,2],extensions=["zip","html"])
    output:
        config["output_dir"]+"rawQC/multiqc_report.html"
    shell:
        """
        multiqc {config[output_dir]}rawQC/ -o {config[output_dir]}rawQC/
        """
