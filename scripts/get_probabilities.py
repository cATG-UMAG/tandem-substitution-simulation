#!/usr/bin/env python3
# Gets proabilites of mutation in each position counting mutations in fasta files
import re
import sys
from os import makedirs, path

import numpy as np
import pandas as pd

from Bio import SeqIO


def main():
    if len(sys.argv) < 3:
        print(
            "./get_probabilities.py <references_fasta> <out_dir> <seq_file_1> [seq_file_2 ...]"
        )
        exit(-1)

    references_fn, outdir = sys.argv[1:3]
    fasta_files = sys.argv[3:]

    makedirs(outdir, exist_ok=True)
    references = {x.id: x for x in SeqIO.parse(references_fn, "fasta")}

    mutations_by_seq = {}

    for f in fasta_files:
        name = re.sub(r"^.*/|[.].*$", "", f)

        seqs = [x for x in SeqIO.parse(f, "fasta")]
        ref = references[name]

        pos_counts = np.zeros(len(ref))
        seq_counts = []

        # to store mutation probability to each base
        base_counts = np.zeros((4, len(ref)))
        b_idx = {"A": 0, "C": 1, "G": 2, "T": 3}

        for s in seqs:
            c = 0
            for i, (x, y) in enumerate(zip(ref, s)):
                if x != y and y in "ACGT":
                    pos_counts[i] += 1
                    c += 1
                    base_counts[b_idx[y], i] += 1
            seq_counts.append(c)

        mutations_by_seq[name] = seq_counts
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

        df.to_csv(path.join(outdir, f"{name}.tsv"), sep="\t", index=False)

    with open(path.join(outdir, "mutations_by_seq.txt"), "w") as f:
        f.write(
            "\n".join(
                "{}\t{}".format(k, ",".join(str(x) for x in mutations_by_seq[k]))
                for k in mutations_by_seq
            )
        )


if __name__ == "__main__":
    main()
