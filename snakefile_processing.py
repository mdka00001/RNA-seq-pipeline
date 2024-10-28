import os
import glob
from snakemake.io import glob_wildcards

configfile : "config.yaml"

# Define the directory containing the raw sequence files

# Use glob_wildcards to match SRA files within subdirectories and extract the SRA identifiers
# Define the final output of the workflow
rule all:
    input:
        expand(config["output_dir"]+"trimmedQC/"+"{sra}_trimmed_{frr}_fastqc.{extensions}",sra=config["accessions"], frr=[1,2],extensions=["zip","html"]),
        expand(config["output_dir"]+"trimmedQC/multiqc_report.html")

rule trim_reads:
    input:
        r1 = config["output_dir"]+"rawSeq/{sra}_1.fastq",
        r2 = config["output_dir"]+"rawSeq/{sra}_2.fastq",
        adapter = config["adapter"]
    output:
        r1 = config["output_dir"]+"trimmedSeq/{sra}_trimmed_1.fq.gz",
        r2 = config["output_dir"]+"trimmedSeq/{sra}_trimmed_2.fq.gz",
        r1_unpaired = config["output_dir"]+"trimmedSeq/{sra}_trimmed_1.unpaired.fq.gz",
        r2_unpaired = config["output_dir"]+"trimmedSeq/{sra}_trimmed_2.unpaired.fq.gz"
    log:
        config["output_dir"]+"trimmedSeq/{sra}.log"
    shell:
        "trimmomatic-0.39.jar PE "
        "-threads 2 "
        "{input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "ILLUMINACLIP:{input.adapter}:2:60:7 SLIDINGWINDOW:2:30"

rule fastqc:
    input:
        input1 = config["output_dir"]+"trimmedSeq/"+"{sra}_trimmed_{frr}.fq.gz"
    output:
        zip = config["output_dir"]+"trimmedQC/"+"{sra}_trimmed_{frr}_fastqc.zip",
        html = config["output_dir"]+"trimmedQC/"+"{sra}_trimmed_{frr}_fastqc.html",
    shell:
        """
        fastqc {input.input1} -o {config[output_dir]}/trimmedQC/
        """

rule multiQC:
    input:
        expand(config["output_dir"]+"trimmedQC/"+"{sra}_trimmed_{frr}_fastqc.html", sra=config["accessions"], frr=[1,2]),
        expand(config["output_dir"]+"trimmedQC/"+"{sra}_trimmed_{frr}_fastqc.zip", sra=config["accessions"], frr=[1,2])
    output:
        config["output_dir"]+"trimmedQC/multiqc_report.html"
    shell:
        """
        multiqc {config[output_dir]}/trimmedQC/ -o {config[output_dir]}/trimmedQC/
        """