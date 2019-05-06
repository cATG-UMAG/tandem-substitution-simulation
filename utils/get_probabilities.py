#!/usr/bin/env python3
# Gets proabilites of mutation in each position counting mutations in fasta files
import re
from glob import glob

import numpy as np
import pandas as pd

from Bio import SeqIO

# pattern to match all the fasta files to find tandems
FASTA_FILES = "fasta/IG*.fasta"
# references fasta file: the fasta filename without extension must match with the sequence name in the references
REFERENCES = "fasta/references.fasta"


def main():
    references = {x.id: x for x in SeqIO.parse(REFERENCES, "fasta")}
    mutations_by_seq = {}

    for f in glob(FASTA_FILES):
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

        df.to_csv("mutation_info/{}.tsv".format(name), sep="\t", index=False)

    with open("mutation_info/mutations_by_seq.txt", "w") as f:
        f.write(
            "\n".join(
                "{}\t{}".format(k, ",".join(str(x) for x in mutations_by_seq[k]))
                for k in mutations_by_seq
            )
        )


if __name__ == "__main__":
    main()
