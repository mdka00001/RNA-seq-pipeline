import os
import glob
from snakemake.io import glob_wildcards

configfile : "config.yaml"

rule all:
    input:
        expand(config["output_dir"]+"aligned/{sra}_paired.bam", sra=config["accession"]),
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}.txt", sra=config["accession"]),
        expand(config["output_dir"]+"featureCounts/featureCounts_{sra}_filtered.txt", sra=config["accession"]),
        config["output_dir"]+"featureCounts/genes.txt",
        config["output_dir"]+"featureCounts/output.txt"

rule hisat2_align:
    input:
        rawread_1 = config["input_dir"]+"trimmedSeq/{sra}_trimmed_1.fq.gz",
        rawread_2 = config["input_dir"]+"trimmedSeq/{sra}_trimmed_2.fq.gz"
    output:
        sam = temp(config["output_dir"]+"aligned/{sra}.sam")
    singularity:
        "docker://quay.io/biocontainers/hisat2:2.2.1--h503566f_8"
    shell:
        """
        hisat2 -p 8 -q --rna-strandness RF -x {config[refseq]} \
            -1 {input.rawread_1} -2 {input.rawread_2} -S {output.sam}
        """

rule sam_to_bam:
    input:
        sam = config["output_dir"]+"aligned/{sra}.sam"
    output:
        bam = config["output_dir"]+"aligned/{sra}_paired.bam"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.6--h5fe306e_12"
    shell:
        """
        samtools sort -o {output.bam} {input.sam}
        """

rule features:
    input:
        rawread_1 = config["output_dir"]+"aligned/{sra}_paired.bam",
        index = config["index"]
    output:
        trimmed_1 = config["output_dir"]+"featureCounts/featureCounts_{sra}.txt"
    threads: 4
    singularity:
        "docker://quay.io/biocontainers/subread:2.1.1--h577a1d6_0"
    shell:
        """
        featureCounts -p -O -T {threads} -a {input.index} -o {output.trimmed_1} {input.rawread_1}
        """

rule filter:
    input:
        rawread = config["output_dir"]+"featureCounts/featureCounts_{sra}.txt"
    output:
        filtered = config["output_dir"]+"featureCounts/featureCounts_{sra}_filtered.txt"
    singularity:
        "docker://quay.io/biocontainers/coreutils:9.5"
    shell:
        """
        cut -f 7 {input.rawread} > {output.filtered}
        """

rule gene:
    input:
        rawread = expand(config["output_dir"]+"featureCounts/featureCounts_{sra}.txt", sra=config["accession"])
    output:
        gene_ = config["output_dir"]+"featureCounts/genes.txt"
    singularity:
        "docker://quay.io/biocontainers/coreutils:9.5"
    shell:
        """
        ls -1 {input.rawread} | head -1 | xargs cut -f1 > {output.gene_}
        """

rule merged:
    input:
        filtered = expand(config["output_dir"] + "featureCounts/featureCounts_{sra}_filtered.txt", sra=config["accession"]),
        genes = config["output_dir"] + "featureCounts/genes.txt"
    output:
        output_ = config["output_dir"] + "featureCounts/output.txt"
    run:
        header = "Geneid\t" + "\t".join(config["accession"])
        filtered_files = " ".join(input.filtered)

        shell(f"""
            paste {input.genes} {filtered_files} > {output.output_}
            echo -e "{header}" | cat - {output.output_} > tmp && mv tmp {output.output_}
        """)
