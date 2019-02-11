"""
Simple script to parse a Bold_Retriever output file along with its source fasta file to create a multifasta of no hits
suitable for BOLD manual submission
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""

import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description='Parses Bold Retriever Script')
    parser.add_argument('-f', type=str,
                        help='Bold Retriever Output File', required=True)
    parser.add_argument('-i', type=str,
                        help='Multifasta of COI contigs', required=True)
    args = parser.parse_args()

    no_hits = parse_output(args.f)
    contig_extractor(no_hits, args.i)


def parse_output(file):
    # Parses a Bold_Retriever output for file lines with nohits and adds the contig names to a text file
    no_hits = []
    with open(file) as f:

        lines = f. readlines()
        for line in lines:
            split = line.split(",")
            if line.startswith("nohit"):
                no_hits.append(split[1])

    new_file = file.replace(".fasta_output.csv", "_nohits.txt")
    with open(new_file, "w") as g:
        for item in no_hits:
            g.write("{}\n".format(item))

    return no_hits


def contig_extractor(no_hits, contigs):
    # Takes the contigs that did not return hits from Bold_Retriever and adds them to a multifasta
    matching_seq = []
    with open(contigs) as f:
        for seq in SeqIO.parse(f, 'fasta'):
            if seq.id in no_hits:
                matching_seq.append(seq)

    new_name = contigs.replace(".fasta", "_nohits.fasta")
    with open(new_name, "w") as g:
        SeqIO.write(matching_seq, g, "fasta")


if __name__ == "__main__":
    main()
