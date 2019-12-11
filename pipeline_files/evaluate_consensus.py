"""
Examines consensus files for length, read count, contamination
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import glob
import argparse


def main():
    parser = argparse.ArgumentParser(description='Examines all the consensus sequences and '
                                                 'generates a multi-fasta and csv with results')
    parser.add_argument('-d', type=str,
                        help='Consensus Parent Directory', required=True)
    args = parser.parse_args()

    evaluate(args.d)


def evaluate(input_directory):
    consensus_directories = sorted(glob.glob(os.path.join(input_directory, "*")))

    summary_lines = []
    curated_fastas = []
    for consensus_dir in consensus_directories:
        files = sorted(glob.glob(os.path.join(consensus_dir, "*")))
        sample = os.path.basename(consensus_dir)
        if len(files) == 5:

            fasta = files[0]
            merger = files[1]
            contamination = files[4]

            with open(fasta) as f:
                sequence = SeqIO.read(f, "fasta")
                sequence.id = sample
                sequence.description = ""

            possible_contamination = False
            with open(contamination) as g:
                data = g.read()
                if data:
                    possible_contamination = True

            correct_length = True
            if len(sequence.seq) != 646:
                correct_length = False

            with open(merger) as h:
                lines = h.readlines()
                size_a = int(lines[15].rstrip().split("=")[1])
                size_b = int(lines[16].rstrip().split("=")[1])

            enough_reads = False
            if size_a + size_b > 40:
                enough_reads = True

            if correct_length and enough_reads and not possible_contamination:
                curated_fastas.append(sequence)

            s = "{},{},{},{},{}\n".format(sample, correct_length, enough_reads, possible_contamination, sequence.seq)
            summary_lines.append(s)

        else:
            s = "{},{},{},{},{}\n".format(sample, "NA", "NA", "NA", "")
            summary_lines.append(s)

    with open("Summary_Output.csv", "w") as f:
        f.write("Sample Name, Correct Length (T/F), >40 Reads? (T/F), Possible Contamination (T/F), CO1 Sequence\n")
        for line in summary_lines:
            f.write(line)

    with open("Curated_Barcodes.fasta", "w") as f:
        SeqIO.write(curated_fastas, f, "fasta")


if __name__ == "__main__":
    main()
