"""
Examines vsearch deprelicated reads and grabs the top two fragments and merges with emboss.
Also reports if contamination likely by checking size of contamination file
by looking at abundance of other reads
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Extracts and merges fragements from vsearch')
    parser.add_argument('-i', type=str,
                        help='Vsearch Fas', required=True)
    parser.add_argument('-o', type=str,
                        help='Emboss Output Dir', required=True)
    args = parser.parse_args()

    extract_contig(args.i, args.o)


def extract_contig(input_file, output_directory):

    # Verify folders exist
    if os.path.isfile(input_file):
        pass
    else:
        print("Missing either {}".format(input_file))
        return

    seqs = []
    for seq in SeqIO.parse(input_file, 'fasta'):
        seqs.append(seq)

    if len(seqs) > 1:
        max_size = 0
        actual_seqs = []
        potential_contamination_seqs = []

        index = 0
        for seq in seqs:
            index += 1
            size = int(seq.id.split("=")[1])

            if index == 1:
                max_size = int(seq.id.split("=")[1])
                actual_seqs.append(seq)
                continue

            if index == 2:
                actual_seqs.append(seq)
                continue

            if size > float(max_size)*0.1:
                potential_contamination_seqs.append(seq)

        file_prefix = os.path.split(output_directory)[1]

        a_sequence = os.path.join(output_directory, file_prefix + "_a.fasta")
        with open(a_sequence, "w") as f:
            SeqIO.write(actual_seqs[0], f, "fasta")

        b_sequence = os.path.join(output_directory, file_prefix + "_b.fasta")
        with open(b_sequence, "w") as f:
            SeqIO.write(actual_seqs[1], f, "fasta")

        output_file_1 = os.path.join(output_directory, file_prefix + ".merger")
        output_file_2 = os.path.join(output_directory, file_prefix + ".fasta")

        contamination_file = os.path.join(output_directory, file_prefix + "_contamination.fasta")
        with open(contamination_file, "w") as f:
            if len(potential_contamination_seqs) > 0:
                SeqIO.write(potential_contamination_seqs, f, "fasta")

        shell_call = "merger -asequence {} -bsequence {} -outfile {} -outseq {}".format(a_sequence,
                                                                                        b_sequence,
                                                                                        output_file_1,
                                                                                        output_file_2)
        print(shell_call)
        subprocess.call(shell_call, shell=True)


if __name__ == "__main__":
    main()
