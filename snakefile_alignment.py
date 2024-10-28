import os
import glob
from snakemake.io import glob_wildcards

configfile : "config.yaml"

# Define the directory containing the raw sequence files

# Use glob_wildcards to match SRA files within subdirectories and extract the SRA identifiers
# Define the final output of the workflow
rule all:
    input:
        expand(config["output_dir"]+"aligned/{sra}_paired.bam", sra=config["accessions"]),
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}.txt", sra=config["accessions"]),
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}_filtered.txt", sra=config["accessions"]),
        expand(config["output_dir"]+"featureCounts/genes.txt"),
        expand( config["output_dir"]+"featureCounts/output.txt")


rule alignment:
    input:
        rawread_1 = config["output_dir"]+"trimmedSeq/{sra}_trimmed_1.fq.gz",
        rawread_2 = config["output_dir"]+"trimmedSeq/{sra}_trimmed_2.fq.gz"
    output:
        aligned = config["output_dir"]+"aligned/{sra}_paired.bam"
        
    params:
        basename=config["output_dir"]+"aligned/{sra}_paired.bam"
    shell:
        """
        hisat2  -p 8 -q --rna-strandness RF -x {config[refseq]} -1 {input.rawread_1} -2 {input.rawread_2} | samtools sort -o {params.basename}
        """
rule features:
    input:
        expand(config["output_dir"]+"aligned/{sra}_paired.bam", sra=config["accessions"]),
        rawread_1 = config["output_dir"]+"aligned/{sra}_paired.bam",
        index = config["index"]
    output:
        trimmed_1 = config["output_dir"]+"featureCounts/featureCounts_{sra}.txt"
    threads:
        4
    shell:
        """
        featureCounts -p -O -T 8 -a {input.index} -o {output.trimmed_1} {input.rawread_1} 
        """

rule filter:
    input:
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}.txt", sra=config["accessions"]),
        rawread = config["output_dir"]+"featureCounts/featureCounts_{sra}.txt"
    output:
        filtered = config["output_dir"]+"featureCounts/featureCounts_{sra}_filtered.txt"
    shell:
        """
        cut -f 7 {input.rawread} > {output.filtered}
        """
rule gene:
    input:
        rawread = expand(config["output_dir"]+"featureCounts/featureCounts_{sra}.txt", sra=config["accessions"])
    output:
        gene_ = config["output_dir"]+"featureCounts/genes.txt"
    shell:
        """
        ls -1  {input.rawread} | head -1 | xargs cut -f1 > {output.gene_}
        """

rule merged:
    input:
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}_filtered.txt", sra=config["accessions"])
    output:
        output_ = config["output_dir"]+"featureCounts/output.txt"
    shell:
        """
        paste {config[output_dir]}featureCounts/genes.txt {config[output_dir]}featureCounts/*_filtered.txt | sed "1d" > {output.output_} |
        awk 'BEGIN{{OFS="\\t"}} NR==1 {{
            for(i=2; i<=NF; i++) {{
                split($i, a, "/")
                split(a[length(a)], b, "_")
                $i = b[1]
            }}
        }} {{print}}' {output.output_} > {output.output_}
        """


