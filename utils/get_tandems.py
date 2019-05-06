#!/usr/bin/env python3
# Finds tandems in sequences from fasta files comparing them to their reference
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

    for f in glob(FASTA_FILES):
        # no path, no extension
        name = re.sub(r"^.*/|[.].*$", "", f)
        ref = references[name]
        seqs = [x for x in SeqIO.parse(f, "fasta")]
        # save all mutations in a 2 level dict: seq_id in 1st level, position in 2nd level
        seq_mutations = {
            x.id: {
                i: dict(ref=ref[i], alt=x[i]) for i in range(len(x)) if ref[i] != x[i]
            }
            for x in seqs
        }

        mutations = []
        for s in seq_mutations:
            for v in get_1d_clusters(sorted(seq_mutations[s])):
                ref = "".join(seq_mutations[s][i]["ref"] for i in sorted(v))
                alt = "".join(seq_mutations[s][i]["alt"] for i in sorted(v))

                mutations.append([s, min(v), ref, alt, len(ref)])

        df = pd.DataFrame(mutations, columns=["seq_id", "pos", "ref", "alt", "size"])
        df.to_csv("tandem_info/real/{}.tsv".format(name), index=False, sep="\t")


def get_1d_clusters(data, stepsize=1):
    """
    Finds clusters of elements in an array. Returns only clusters with more than 1 element.
    :param data: input array
    :param stepsize: minimun distance between elements to cluster them
    :return: a list of arrays with the clusters found
    """
    consecutive = np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)
    return [x for x in consecutive if len(x) > 1]


if __name__ == "__main__":
    main()
