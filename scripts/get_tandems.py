#!/usr/bin/env python3
# Finds tandems in sequences from fasta files comparing them to their reference
import re
import sys

import numpy as np
import pandas as pd

from Bio import SeqIO


def main():
    if len(sys.argv) < 4:
        print("./get_tandems.py <fasta_file> <references_fasta> <output>")
        exit(-1)

    fasta_fn, references_fn, output_fn = sys.argv[1:]

    references = {x.id: x for x in SeqIO.parse(references_fn, "fasta")}

    # no path/no extension, this name will be the one to search inside reference files
    name = re.sub(r"^.*/|[.].*$", "", fasta_fn)
    ref = references[name]
    seqs = [x for x in SeqIO.parse(fasta_fn, "fasta")]
    # save all mutations in a 2 level dict: seq_id in 1st level, position in 2nd level
    seq_mutations = {
        x.id: {i: dict(ref=ref[i], alt=x[i]) for i in range(len(x)) if ref[i] != x[i]}
        for x in seqs
    }

    mutations = []
    for s in seq_mutations:
        for v in get_1d_clusters(sorted(seq_mutations[s])):
            ref = "".join(seq_mutations[s][i]["ref"] for i in sorted(v))
            alt = "".join(seq_mutations[s][i]["alt"] for i in sorted(v))

            mutations.append([s, min(v), ref, alt, len(ref)])

    df = pd.DataFrame(mutations, columns=["seq_id", "pos", "ref", "alt", "size"])
    df.to_csv(output_fn, index=False, sep="\t")


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
