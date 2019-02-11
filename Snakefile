 Snakemake file to process CO1 Barcoding Genes Sequenced on a MiSeq
# Author: Jackson Eyres jackson.eyres@canada.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

import glob
import os
from shutil import copyfile

# Configuration Settings

# Location of fastq folder, default "fastq". Phyluce requires files to not mix "-" and "_", so fastq files renamed
for f in glob.glob('fastq/*.fastq.gz'):
    basename = os.path.basename(f)
    new_basename = ""

SAMPLES = set([os.path.basename(f).replace("_L001_R1_001.fastq.gz","").replace("_L001_R2_001.fastq.gz","") for f in glob.glob('fastq/*.fastq.gz')])

# Location of adaptor.fa for trimming
adaptors = "pipeline_files/adapters.fa"


rule all:
    input:
        r1_trimmed = expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2_trimmed = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES),

        # fastqc_trimmed_dir = directory("fastqc_trimmed"),

        spades_assemblies = expand("spades_assemblies/{sample}/contigs.fasta", sample=SAMPLES),

        all_spades_assemblies = directory("all_spades_assemblies"),
        final_good_contigs = "final_good_contigs.fasta",
        final_medium_contigs = "final_medium_contigs.fasta",
        final_good_contigs_aligned = "final_good_contigs_aligned.fasta",
        final_medium_contigs_aligned = "final_medium_contigs_aligned.fasta",

        bold_good = "final_good_contigs_aligned.fasta_output.csv",
        nohits_good = "final_good_contigs_aligned_nohits.txt"

rule fastqc:
    # Quality Control check on raw data before adaptor trimming
    input:
        directory("fastq")
    output:
         fastqc_dir = directory("fastqc")
    log: "logs/fastqc.log"
    conda: "pipeline_files/barcoding.yml"
    shell:
        "mkdir fastqc; fastqc -o fastqc fastq/*.fastq.gz 2>{log} 2>&1"


rule bbduk:
    # Sequencing Adaptor trimming
    input:
        r1 = 'fastq/{sample}_L001_R1_001.fastq.gz',
        r2 = 'fastq/{sample}_L001_R2_001.fastq.gz'
    output:
        out1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        out2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz",
    log: "logs/bbduk.{sample}.log"
    conda: "pipeline_files/barcoding.yml"
    shell: "bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref={adaptors} qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>{log} 2>&1; touch {output.out1} {output.out2}"


rule fastqc_trimmed:
    # Quality Control check after adaptor trimming
    input: r1=expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2 = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES)
    output:
        fastqc_trimmed_dir = directory("fastqc_trimmed")
    log: "logs/fastqc_trimmed.log"
    conda: "pipeline_files/barcoding.yml"
    shell:
        "mkdir fastqc_trimmed; fastqc -o fastqc_trimmed trimmed/*.fastq.gz 2>{log} 2>&1"


rule spades:
    # Assembles fastq files using default settings
    input:
        r1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        r2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz"
    output:
        "spades_assemblies/{sample}/contigs.fasta"
    log: "logs/spades.{sample}.log"
    conda: "pipeline_files/barcoding.yml"
    threads: 2
    shell:
        "spades.py -t 2 -1 {input.r1} -2 {input.r2} -o spades_assemblies/{wildcards.sample} 2>{log} 2>&1"


rule gather_assemblies:
    # Rename all spades assemblies and copy to a folder for further analysis
    input:
        assemblies = expand("spades_assemblies/{sample}/contigs.fasta", sample=SAMPLES)
    output:
        all_spades_assemblies = directory("all_spades_assemblies")
    run:
        for assembly in input.assemblies:
            if os.path.exists(assembly):
                if os.path.exists("all_spades_assemblies"):
                    pass
                else:
                    os.makedirs("all_spades_assemblies")
            newname = assembly.split("/")[1] + "_S.fasta"
            copyfile(assembly, os.path.join(output.all_spades_assemblies, newname))


rule fastq_quality_metrics:
    # BBMap's Stats.sh assembly metrics for fastq files
    input: directory("fastq")
    output: "fastq_metrics.tsv"
    conda: "pipeline_files/barcoding.yml"
    shell: "statswrapper.sh {input}/*.fastq.gz > {output}"


rule gather_contigs:
    input: directory("all_spades_assemblies")
    output: "final_good_contigs.fasta", "final_medium_contigs.fasta", directory("problem_fastas")
    conda: "pipeline_files/barcoding.yml"
    shell: "python pipeline_files/gather_contigs.py -s {input}"


rule align_to_reference:
    input: "final_good_contigs.fasta"
    output: "final_good_contigs_aligned.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell: "cp pipeline_files/co1.fasta temp.fasta && cat final_good_contigs.fasta >> temp.fasta && mafft --adjustdirection temp.fasta > final_good_contigs_aligned.fasta && rm temp.fasta"


rule align_medium_to_reference:
    input: "final_medium_contigs.fasta"
    output: "final_medium_contigs_aligned.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell: "cp pipeline_files/co1.fasta temp.fasta && cat final_medium_contigs.fasta >> temp.fasta && mafft --adjustdirection temp.fasta > final_medium_contigs_aligned.fasta && rm temp.fasta"


rule align_poor_to_reference:
    input: directory("problem_fastas")
    output: directory("problem_fastas")
    conda: "pipeline_files/barcoding.yml"
    shell: "for f in problem_fastas/*; do cp pipeline_files/co1.fasta temp.fasta && cat $f >> temp.fasta && mafft --adjustdirection temp.fasta > problem_fastas_aligned/${f//+(*\/|.*)} && rm temp.fasta; done"


rule bold_retriever:
    input: "final_good_contigs_aligned.fasta"
    output: "final_good_contigs_aligned.fasta_output.csv"
    conda: "pipeline_files/barcoding.yml"
    shell: "python pipeline_files/bold_retriever-master/bold_retriever.py -f final_good_contigs_aligned.fasta -db COX1_SPECIES"


rule bold_parser:
    input: "final_good_contigs_aligned.fasta_output.csv"
    output: "final_good_contigs_aligned_nohits.txt"
    conda: "pipeline_files/barcoding.yml"
    shell: "python pipeline_files/bold_retriever_parser.py -f final_good_contigs_aligned.fasta_output.csv"