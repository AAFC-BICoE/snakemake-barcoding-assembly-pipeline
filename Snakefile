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
        # When SPAdes fails, it wont create contigs.fasta, but should create input_dataset.yaml
        spades_datasets = expand("spades_assemblies/{sample}/input_dataset.yaml", sample=SAMPLES),

        spades_assemblies_temp = expand("spades_assemblies/{sample}/contigs_temp.fasta", sample=SAMPLES),
        spades_assemblies = expand("spades_assemblies/{sample}/contigs.fasta", sample=SAMPLES),

        assemblies_renamed = expand("all_spades_assemblies/{sample}_S.fasta", sample=SAMPLES),

        final_good_contigs = "final_good_contigs.fasta",
        final_medium_contigs = "final_medium_contigs.fasta",
        final_good_contigs_aligned = "final_good_contigs_aligned.fasta",
        final_medium_contigs_aligned = "final_medium_contigs_aligned.fasta",
#
        problem_contigs = "problem_fastas.txt",
#        problem_fasta_directory = directory("problem_fastas_aligned")


rule bbduk:
    # Sequencing Adaptor and quality trimming
    input:
        r1 = 'fastq/{sample}_L001_R1_001.fastq.gz',
        r2 = 'fastq/{sample}_L001_R2_001.fastq.gz'
    output:
        out1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        out2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz",
    log: "logs/bbduk.{sample}.log"
    conda: "pipeline_files/barcoding.yml"
    shell: "bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref={adaptors} qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo &>{log}; touch {output.out1} {output.out2}"


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
    conda: "pipeline_files/barcoding.yml"
    shell: "bbmerge.sh in={input.r1} in2={input.r2} outm={output.merged} outu={output.unmerged} &>{log}; touch {output.merged} {output.unmerged}"


rule remove_primers:
    # Custom script to clean up reads of any degenerate primers and spurious sequencing bases
    input:
        "merged/{sample}_merged.fq"
    output:
        "sliced/{sample}_sliced.fq"
    log: "logs/slice.{sample}.log"
    conda: "pipeline_files/barcoding.yml"
    shell: "python pipeline_files/trim_primers.py -f {input} -p {primers} -o {output} &>{log}"


rule spades:
    # Assembles spliced fastq files using default settings
    input:
        s = "sliced/{sample}_sliced.fq",
    output:
        "spades_assemblies/{sample}/input_dataset.yaml"
    log: "logs/spades.{sample}.log"
    conda: "pipeline_files/barcoding.yml"
    # On certain samples, SPAdes fails and produces a non zero exit code which causes snakemake to end prematurely.
    # This code is subsequently ignored in the shell command and converted to a zero exit code.
    # Failed SPAdes runs are excluded from further analysis
    threads: 4
    shell:
        "spades.py -t 4 -s {input.s} -o spades_assemblies/{wildcards.sample}  &>{log}|| exit 0"


rule spades_touch:
    # Just to workaround snakemake restrictions with failed spades runs
    input:
        "spades_assemblies/{sample}/input_dataset.yaml"
    output:
        "spades_assemblies/{sample}/contigs_temp.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell:
        "touch spades_assemblies/{wildcards.sample}/contigs.fasta; cp spades_assemblies/{wildcards.sample}/contigs.fasta {output}"


rule spades_touch_2:
    # Just to workaround snakemake restrictions with failed spades runs
    input:
        "spades_assemblies/{sample}/contigs_temp.fasta"
    output:
        "spades_assemblies/{sample}/contigs.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell:
        "cp {input} {output}"


rule gather_assemblies:
    # Rename all spades assemblies and copy to a folder for further analysis
    input:
        assemblies = expand("spades_assemblies/{sample}/contigs.fasta", sample=SAMPLES)
    output:
        assemblies_renamed = expand("all_spades_assemblies/{sample}_S.fasta", sample=SAMPLES)
    run:
        for assembly in input.assemblies:
            if os.path.exists(assembly):
                if os.path.exists("all_spades_assemblies"):
                    pass
                else:
                    os.makedirs("all_spades_assemblies")
            newname = assembly.split("/")[1] + "_S.fasta"
            copyfile(assembly, os.path.join("all_spades_assemblies", newname))


rule gather_contigs:
    # Custom script to seperate out contigs into high, medium and low quality
    input: expand("all_spades_assemblies/{sample}_S.fasta", sample=SAMPLES)
    output: "final_good_contigs.fasta", "final_medium_contigs.fasta", "problem_fastas.txt"
    conda: "pipeline_files/barcoding.yml"
    shell: "python pipeline_files/gather_contigs.py -s all_spades_assemblies"


rule align_to_reference:
    # Orients contigs to a COI reference for 5' to 3' for BOLD submission
    input: "final_good_contigs.fasta"
    output: "final_good_contigs_aligned.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell: "cp pipeline_files/co1.fasta temp.fasta && cat final_good_contigs.fasta >> temp.fasta && mafft --adjustdirection temp.fasta > final_good_contigs_aligned.fasta && rm temp.fasta"


rule align_medium_to_reference:
    # Orients contigs to a COI reference for 5' to 3' for BOLD submission
    input: "final_medium_contigs.fasta",
    output: "final_medium_contigs_aligned.fasta"
    conda: "pipeline_files/barcoding.yml"
    shell: "cp pipeline_files/co1.fasta temp1.fasta && cat final_medium_contigs.fasta >> temp1.fasta && mafft --adjustdirection temp1.fasta > final_medium_contigs_aligned.fasta && rm temp1.fasta"


#rule align_poor_to_reference:
#    # Orients contigs to a COI reference for 5' to 3' for BOLD submission
#    input: "problem_fastas.txt"
#    output: directory("problem_fastas_aligned")
#    conda: "pipeline_files/barcoding.yml"
#    shell: "for f in problem_fastas/*; do cp pipeline_files/co1.fasta temp.fasta && cat $f >> temp.fasta && mafft --adjustdirection temp.fasta > problem_fastas_aligned/$f//+(*\/|.*)} && rm temp.fasta; done"

