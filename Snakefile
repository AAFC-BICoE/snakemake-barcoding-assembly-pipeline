# Snakemake file to process CO1 Barcoding Genes Sequenced on a MiSeq
# Author: Jackson Eyres jackson.eyres@canada.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

import glob
import os
from shutil import copyfile

# Configuration Settings
SAMPLES = set([os.path.basename(f).replace("_L001_R1_001.fastq.gz","").replace("_L001_R2_001.fastq.gz","") for f in glob.glob('fastq/*.fastq.gz')])

# Location of adaptor.fa for trimming
adaptors = "pipeline_files/adapters.fa"
primers = "pipeline_files/diptera_primers.fasta"

rule all:
    input:
        r1_trimmed = expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2_trimmed = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES),

        merged = expand("merged/{sample}_merged.fq", sample=SAMPLES),
        unmerged = expand("unmerged/{sample}_unmerged.fq", sample=SAMPLES),

        sliced = expand("sliced/{sample}_sliced.fq", sample=SAMPLES),
        vsearch = expand("vsearch/{sample}.fas", sample=SAMPLES),
        consensus = expand("consensus/{sample}/{sample}.fasta", sample=SAMPLES),

        summary_output = "Summary_Output.csv",
        multifasta = "Curated_Barcodes.fasta"

rule bbduk:
    # Sequencing Adaptor and quality trimming
    input:
        r1 = 'fastq/{sample}_L001_R1_001.fastq.gz',
        r2 = 'fastq/{sample}_L001_R2_001.fastq.gz'
    output:
        out1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        out2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz",
    log: "logs/bbduk.{sample}.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref={adaptors} qtrim=rl trimq=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &>{log}; touch {output.out1} {output.out2}"


rule bbmerge:
    # Merges paired end reads with overlapping regions into a single long fragement.
    # Useful for amplicon based sequencing
    input:
        r1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        r2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz"
    output:
        merged = "merged/{sample}_merged.fq",
        unmerged = "unmerged/{sample}_unmerged.fq"
    log: "logs/bbmerge.{sample}.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "bbmerge.sh in={input.r1} in2={input.r2} outm={output.merged} outu={output.unmerged} &>{log}; touch {output.merged} {output.unmerged}"


rule remove_primers:
    # Custom script to clean up reads of any degenerate primers and spurious sequencing bases
    input:
        "merged/{sample}_merged.fq"
    output:
        "sliced/{sample}_sliced.fq"
    log: "logs/slice.{sample}.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "python pipeline_files/trim_primers.py -f {input} -p {primers} -o {output} &>{log}"


rule vsearch:
    # Vsearch removes duplicate reads
    input:
        "sliced/{sample}_sliced.fq"
    output:
        "vsearch/{sample}.fas"
    log: "logs/vsearch.{sample}.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "vsearch --derep_fulllength {input} --sizein --fasta_width 0 --sizeout --output {output} &>{log} || true; touch {output}"


rule evaluate_vsearch:
    # Vsearch removes duplicate reads
    input:
        "vsearch/{sample}.fas"
    output:
        consensus="consensus/{sample}/{sample}.fasta",
        contam="consensus/{sample}/{sample}_contamination.fasta"
    log: "logs/consensus.{sample}.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "python pipeline_files/evaluate_vsearch.py -i {input} -o consensus/{wildcards.sample} &>{log}; touch {output.consensus} {output.contam}"


rule evaluate_consensus:
    input:
        consensus = expand("consensus/{sample}/{sample}.fasta", sample=SAMPLES),
    output:
        summary="Summary_Output.csv",
        multifasta="Curated_Barcodes.fasta"
    log: "logs/evaluate_consensus.log"
    conda: "pipeline_files/vsearch_env.yml"
    shell: "python pipeline_files/evaluate_consensus.py -d consensus &>{log}"
