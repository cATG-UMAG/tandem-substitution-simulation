#!/usr/bin/env python3
import re
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np


STOP_CODONS = ['TAA', 'TAG', 'TGA']


def main():
    if len(sys.argv) < 6:
        print("./tandem_generation.py <mutation_info_file> <reference_seq_file> <mutations_byseq_file> <N> <outfile>")
        sys.exit(-1)

    # read arguments
    mutation_info, ref_file, mutations_by_seq, n, outfile = sys.argv[1:7]
    n = int(n)
    target_name = re.sub('^.*/|[.].*$', '', mutation_info)

    # load data
    references = {x.id: x for x in SeqIO.parse(ref_file, 'fasta')}
    with open(mutations_by_seq) as f:
        max_mutations_byseq = {x.split('\t')[0]: max([int(y) for y in x.split('\t')[1].split(',')]) for x in f.read().splitlines()}

    m_info = pd.read_csv(mutation_info)
    ref = references[target_name]
    max_mutations = max_mutations_byseq[target_name]
    n_size = len(str(n))

    tandem_info = []
    # generate sequences with a fixed number of mutations, from 2 to max_mutations
    for i in range(2, max_mutations + 1):
        for j in range(n):
            mutations = generate_mutations(ref, m_info.mutation_probability, i)
            for v in get_1d_clusters(sorted(mutations)):
                tandem = ''.join(mutations[m] for m in v)
                t_size = len(v)
                # check if a stop codon is produced
                mutated_subseq = re.sub('[.]+', tandem, get_context(ref, v[0], 2, t_size))
                stop = any(s in mutated_subseq for s in STOP_CODONS)

                tandem_info.append([i, str(j+1).zfill(n_size), v[0], ref[v[0]: v[0]+t_size].seq, tandem, t_size, stop])

    tandem_df = pd.DataFrame(tandem_info, columns=['mutations_byseq', 'seq_id', 'pos', 'ref', 'alt', 'size', 'stop_codon'])
    tandem_df.to_csv(outfile + '.gz', compression='gzip', sep='\t', index=False)


# generate a set of mutations on a sequence based on a vector of probabilities
def generate_mutations(ref, mutation_probabilities, n):
    seq_len = len(ref)
    # working with numpy arrays to speed up the proccess
    positions = np.arange(0, seq_len)
    mutation_prob = np.array(mutation_probabilities)
    mutations = {}
    while len(mutations) < n:
        mutation_candidates = np.empty(0)
        # try generating different random arrays if no mutations candidates were selected
        while len(mutation_candidates) == 0:
            random = np.random.rand(seq_len)
            mutation_candidates = positions[random <= mutation_prob]
        # select randomly one of the mutation candidates and generate a mutation
        pos = int(np.random.choice(mutation_candidates))
        mutated_base = mutate(ref[pos])
        mutations[pos] = mutated_base

    return mutations


# choose randomly a different base than the one from the argument
def mutate(base):
    return np.random.choice([x for x in "ACGT" if x != base])


# gets context from a sequence at a given position
# n defines how many bases are needed by each side
# s defines the size of the central region
def get_context(seq, pos, n=1, s=1):
    return ''.join(['.' if pos <= n < pos+s else seq[n] if 0 <= n < len(seq) else '' for n in range(pos-n, pos+n+s)])


# given a array, finds clusters of elmeents on it.
# stepsize defines the minimun distance between elements
def get_1d_clusters(data, stepsize=1):
    consecutive = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    # only return cluters with size > 1
    return [x for x in consecutive if len(x) > 1]


if __name__ == '__main__':
    main()
