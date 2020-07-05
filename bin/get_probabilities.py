#!/usr/bin/env python3
import argparse
import re
import sys
from os import makedirs, path
from typing import List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    args = parse_arguments()

    if args.output_dir:
        makedirs(args.outdir, exist_ok=True)
        output_dir = args.output_dir
    else:
        output_dir = ""

    references = {x.id: x for x in SeqIO.parse(args.references, "fasta")}
    seqs = [x for x in SeqIO.parse(args.sequences, "fasta")]

    if args.target_name:
        name = args.target_name
    else:
        # get name directly from filename
        name = re.sub(r"^.*/|[.].*$", "", args.sequences)
    ref = references[name]

    pos_counts, base_counts, mutations_by_seq = get_mutation_frequencies(seqs, ref)

    # format all information in a dataframe
    df = pd.DataFrame(
        {
            "position": range(1, len(ref) + 1),
            "mutation_count": pos_counts.astype(int),
            "mutation_probability": pos_counts / len(seqs),
        }
    )
    df_bases = pd.DataFrame(
        np.divide(
            base_counts,
            pos_counts,
            out=np.zeros_like(base_counts),
            where=pos_counts != 0,
        ).T,
        columns=["A", "C", "G", "T"],
    )
    df = pd.concat([df, df_bases], axis=1, sort=False)

    # save files
    df.to_csv(path.join(output_dir, f"{name}.tsv"), sep="\t", index=False)

    with open(path.join(output_dir, args.mutations_by_seq), "w") as f:
        f.write(f"{name}\t{','.join(str(x) for x in mutations_by_seq)}")


def get_mutation_frequencies(
    sequences: List[Seq], reference: Seq
) -> Tuple[np.array, np.array, List[int]]:
    """
    Gets mutations frequencies from a list of sequences comparing againt a reference

    :param sequences: list of sequences to find mutations
    :param reference: reference sequence
    :return: mutational frequencies vectors and list of mutations by sequence
    """
    b_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    # mutational frequency of each position
    pos_counts = np.zeros(len(reference))
    # mutational frequency of each base on each position
    base_counts = np.zeros((4, len(reference)))
    # count of mutations for each sequence
    seq_mutation_counts = []

    for s in sequences:
        c = 0
        for i, (x, y) in enumerate(zip(reference, s)):
            if x != y and y in "ACGT":
                pos_counts[i] += 1
                c += 1
                base_counts[b_idx[y], i] += 1
        seq_mutation_counts.append(c)

    return pos_counts, base_counts, seq_mutation_counts


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Gets frequecies of mutation in each position counting mutations from a FASTA file"
    )

    parser.add_argument(
        "sequences",
        metavar="sequences_file",
        help="FASTA file with sequences to get mutational frequencies from",
    )

    parser.add_argument(
        "references",
        metavar="references_file",
        help="FASTA file with references to detect mutations",
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        metavar="output_dir",
        action="store",
        help="directory to save output files",
    )

    parser.add_argument(
        "--target-name",
        "-n",
        metavar="target_name",
        action="store",
        help="use this name to get reference and in output files instead of sequences filename",
    )

    parser.add_argument(
        "--mutations-by-seq",
        "-m",
        metavar="mutations_by_seq",
        action="store",
        default="mutations_by_seq.txt",
        help="mutations by sequence output filename"
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
