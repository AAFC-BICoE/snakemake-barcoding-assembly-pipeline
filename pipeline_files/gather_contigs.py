"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import glob
import argparse
import shutil


def main():
    parser = argparse.ArgumentParser(description='Creates a multifasta from the best SPAdes contigs')
    parser.add_argument('-s', type=str,
                        help='Folder containing SPAdes assemblies', required=True)
    args = parser.parse_args()

    combine_uces(args.s)


def combine_uces(spades_directory):
    """
    Takes the UCES from an rnaSPAdes and SPAdes run and creates a seperate file taking only the best sequence per UCE
    :return:
    """

    # Verify folders exist
    if os.path.isdir(spades_directory):
        pass
    else:
        print("Missing either {}".format(spades_directory))
        return

    spades_fastas = glob.glob(os.path.join(spades_directory, "*.fasta"))

    # Put all the contigs into a single dictionary
    final_good_contigs = []
    final_medium_contigs = []
    problematic_fastas = []
    for fasta in spades_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace(".fasta", "")
        contigs = []
        for seq in SeqIO.parse(fasta, 'fasta'):
            contigs.append(seq)

        temp_final_good_contigs = []
        temp_final_medium_contigs = []
        for seq in contigs:
            length = len(seq.seq)
            coverage = float(seq.description.split("_")[-1])
            seq.id = specimen_name

            if coverage > 50 and length > 600:
                temp_final_good_contigs.append(seq)
            elif coverage > 50 and 400 < length < 601:
                temp_final_medium_contigs.append(seq)

        if len(temp_final_good_contigs) == 1:
            final_good_contigs.append(temp_final_good_contigs[0])

        elif len(temp_final_medium_contigs) == 1:
            final_medium_contigs.append(temp_final_medium_contigs[0])
        else:
            problematic_fastas.append(fasta)

    print(len(spades_fastas))
    print("Good Quality Contigs: ", len(final_good_contigs))
    print("Medium Quality Contigs: ", len(final_medium_contigs))
    print("Problem Fastas: ", len(problematic_fastas))

    problematic_path = "problem_fastas.txt"
    os.makedirs("problem_fastas")
    with open(problematic_path, "w") as e:
        for item in problematic_fastas:
            e.write(os.path.split(item)[1] + "\n")
            new_path = os.path.join("problem_fastas", os.path.split(item)[1])
            shutil.copy(item, new_path)

    # Write Final Contigs to files
    final_good_path = "final_good_contigs.fasta"
    with open(final_good_path, "w") as f:
        for seq in final_good_contigs:
            SeqIO.write(seq, handle=f, format="fasta")

    final_medium_path = "final_medium_contigs.fasta"
    with open(final_medium_path, "w") as g:
        for seq in final_medium_contigs:
            SeqIO.write(seq, handle=g, format="fasta")


if __name__ == "__main__":
    main()
