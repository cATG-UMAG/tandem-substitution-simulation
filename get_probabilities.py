#!/usr/bin/env python3
import re
from glob import glob
from Bio import SeqIO
import pandas as pd
import numpy as np


def main():
    ref_file = "fasta/references.fasta"

    references = {x.id: x for x in SeqIO.parse(ref_file, 'fasta')}
    mutations_by_seq = {}
    for f in glob("fasta/IG*.fasta"):
        name = re.sub('^.*/|[.].*$', '', f)

        seqs = [x for x in SeqIO.parse(f, 'fasta')]
        ref = references[name]

        pos_counts = np.zeros(len(ref))
        seq_counts = []

        for s in seqs:
            c = 0
            for i, (x, y) in enumerate(zip(ref, s)):
                if x != y:
                    pos_counts[i] += 1
                    c += 1
            seq_counts.append(c)

        mutations_by_seq[name] = seq_counts
        df = pd.DataFrame({
            'position': range(1, len(ref) + 1),
            'mutation_count': pos_counts.astype(int),
            'mutation_probability': pos_counts / len(seqs)})

        df.to_csv("mutation_info/{}.tsv".format(name), index=False)

    with open('mutation_info/mutations_by_seq.txt', 'w') as f:
        f.write('\n'.join("{}\t{}".format(k, ','.join(str(x) for x in mutations_by_seq[k])) for k in mutations_by_seq))


if __name__ == '__main__':
    main()
