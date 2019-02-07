with open ("problem_fastas.txt") as f:
    lines = f.readlines()
    for line in lines:
        specimen = line.replace("_S.fasta\n", "")
        print("cp {}_L001_R1_001.fastq.gz problem_specimens/{}_L001_R1_001.fastq.gz".format(specimen, specimen))
        print("cp {}_L001_R2_001.fastq.gz problem_specimens/{}_L001_R2_001.fastq.gz".format(specimen, specimen))