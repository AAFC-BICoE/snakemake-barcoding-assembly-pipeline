"""
Script that takes in dual fragment merged paired end reads of Diptera CO1 and removed degenerate primers
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""

from Bio import SeqIO
import subprocess
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='Trims Diptera CO1 merged reads of degenerate primers')

    parser.add_argument('-f', type=str,
                        help='Merged Read file in .fq format', required=True)
    parser.add_argument('-p', type=str,
                        help='Primer File in .fasta format', required=True)
    parser.add_argument('-o', type=str,
                        help='Output file', required=True)
    parser.add_argument('-dg', type=str,
                        help='Path to degenerate oligo Script', default="pipeline_files/dg")
    args = parser.parse_args()

    find_primers(args.f, args.p, args.o, args.dg)


def get_primers(file):
    # Adds primers to list, reverse compliments the reverse primers
    primer_sequences = []
    with open(file) as f:
        index = 0
        for seq in SeqIO.parse(f, "fasta"):
            index += 1
            if index % 2 == 0:
                primer_sequences.append(seq.reverse_complement())
            else:
                primer_sequences.append(seq)

    return primer_sequences


def find_primers(fastq_file, primer_file, output_file, dg):
    """
    Scans reads for forward or reverse primers. If fragment A, removes everything before the forward primer, and
    removes the degenerate reverse primer. If fragement B, removed the degenerate forward primer, and everything beyond
    the reverse primer. Also eliminates any reads that are too short, or have no primers detected. The remaining reads
    can properly assembly without any internal degeneracy.

    -----------------------
    Forward A       --------------------
                               Reverse B
    :param fastq_file:
    :param primer_file:
    :param output_file:
    :param dg: Script to generate oligos of degenerate primers
    :return:
    """
    with open(fastq_file) as f:
        primers = get_primers(primer_file)
        primer_a_f = str(primers[0].seq) # Forward primer Fragment A
        primer_b_r = str(primers[3].seq) # Reverse Primer Fragement B

        # Degenerate primers must be converted into their respective oligos for detection in the reads
        dg_result = subprocess.run(['./{}'.format(dg), primer_b_r], stdout=subprocess.PIPE).stdout.decode('utf-8')
        reverse_primers = dg_result.split('\n')[:-1]

        curated_fragment_a_reads = []
        curated_fragment_b_reads = []
        total_reads = 0
        too_short_a = 0
        too_short_b = 0
        no_primers = 0

        for seq in SeqIO.parse(f, "fastq"):
            total_reads += 1
            """
            Look for fragement A primer in read. If found, trim the read to the minimum_fragement length, 
            which excludes the degenerate reverse primer.
            """
            index = seq.seq.find(primer_a_f)
            if index > 0:
                min_size = 454
                if len(seq.seq) >= index + min_size:  # Fragment is minimum size without including degenerate primer
                    cut_read = seq[index:index+min_size+1]
                    curated_fragment_a_reads.append(cut_read)
                else:
                    too_short_a += 1
            else:
                """
                Since fragement a forward primer wasn't found, likely dealing with fragement B read
                Look for the reverse primer, since its degenerate look at all possible oligos. 
                Then trim the read to just the reverse primer at tail, and trim the first 20 bases to 
                capture the degenerate forward primer
                """

                for primer in reverse_primers:
                    temp_index = seq.seq.find(primer)
                    if temp_index > index:
                        index = temp_index
                        if len(seq.seq) > 450:
                            cut_read = seq[20:index+23]
                            curated_fragment_b_reads.append(cut_read)
                        else:
                            too_short_b += 1
                        break
                else:
                    no_primers += 1
        log_string = "File {}, Total Reads: {}, Missing Primers: {}, Reads Too Short: {}, Fragment A Reads: {}, " \
                     "Fragement B Reads: {}".format(os.path.basename(fastq_file), total_reads, no_primers,
                                                    (too_short_a+too_short_b),
                                                    len(curated_fragment_a_reads), len(curated_fragment_b_reads))
        print(log_string)

        "Write sliced reads to output file"
        curated_reads = curated_fragment_a_reads + curated_fragment_b_reads
        with open(output_file, "w") as g:
            SeqIO.write(curated_reads, g, "fastq")


if __name__ == "__main__":
    main()
