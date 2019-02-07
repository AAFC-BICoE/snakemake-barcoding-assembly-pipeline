"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description='Parses Bold Retriever Script')
    parser.add_argument('-f', type=str,
                        help='Bold Retriever Output File', required=True)
    args = parser.parse_args()

    parse_output(args.f)


def parse_output(file):

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


if __name__ == "__main__":
    main()
