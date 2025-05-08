configfile: "config.yaml"

rule all:
    input:
        expand(
            config["output_dir"] + "trimmedQC/{sra}_trimmed_{frr}_fastqc.{ext}",
            sra=config["accession"], frr=[1, 2], ext=["zip", "html"]
        ),
        config["output_dir"] + "trimmedQC/multiqc_report.html"


rule trim_reads:
    input:
        r1 = config["input_dir"] + "{sra}/{sra}_1.fq.gz",
        r2 = config["input_dir"] + "{sra}/{sra}_2.fq.gz",
        adapter = config["adapter"]
    output:
        r1 = config["output_dir"] + "trimmedSeq/{sra}_trimmed_1.fq.gz",
        r2 = config["output_dir"] + "trimmedSeq/{sra}_trimmed_2.fq.gz",
        r1_unpaired = config["output_dir"] + "trimmedSeq/{sra}_trimmed_1.unpaired.fq.gz",
        r2_unpaired = config["output_dir"] + "trimmedSeq/{sra}_trimmed_2.unpaired.fq.gz"
    log:
        config["output_dir"] + "trimmedSeq/{sra}.log"
    singularity:
        "docker://quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    shell:
        """
        trimmomatic PE \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:{input.adapter}:2:60:7 \
            SLIDINGWINDOW:4:5 \
            MINLEN:36 \
            AVGQUAL:6
        """

rule fastqc:
    input:
        fq = config["output_dir"] + "trimmedSeq/{sra}_trimmed_{frr}.fq.gz"
    output:
        zip = config["output_dir"] + "trimmedQC/{sra}_trimmed_{frr}_fastqc.zip",
        html = config["output_dir"] + "trimmedQC/{sra}_trimmed_{frr}_fastqc.html"
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    shell:
        """
        fastqc {input.fq} -o {config[output_dir]}trimmedQC/
        """

rule multiQC:
    input:
        expand(config["output_dir"] + "trimmedQC/{sra}_trimmed_{frr}_fastqc.html",
               sra=config["accession"], frr=[1, 2]),
        expand(config["output_dir"] + "trimmedQC/{sra}_trimmed_{frr}_fastqc.zip",
               sra=config["accession"], frr=[1, 2])
    output:
        config["output_dir"] + "trimmedQC/multiqc_report.html"
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0"
    shell:
        """
        multiqc {config[output_dir]}trimmedQC/ -o {config[output_dir]}trimmedQC/
        """
