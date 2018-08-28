#!/usr/bin/env python3
import re
import sys
from multiprocessing import Pool
from itertools import chain
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
        for line in f.read().splitlines():
            line = line.split('\t')
            if line[0] == target_name:
                mutations_by_seq = [int(x) for x in line[1].split(',')]
                break

    m_info = pd.read_csv(mutation_info)
    ref = references[target_name]
    n_size = len(str(n))

    # runs in parallel n simulations with len(mutations_by_seq) sequences
    with Pool() as p:
        tandem_info = p.starmap(run_simulation, ((i, m_info, ref, mutations_by_seq, n_size) for i in range(n)))

    # flatten list of lists to a single list
    tandem_info = list(chain.from_iterable(tandem_info))

    tandem_df = pd.DataFrame(tandem_info, columns=['simulation', 'mutations_byseq', 'seq_id', 'pos', 'ref', 'alt', 'size', 'stop_codon'])
    tandem_df.to_csv(outfile + '.gz', compression='gzip', sep='\t', index=False)


# "main" operation wrapped in a function, to being able to paralelize it
def run_simulation(sim_id, m_info, ref, mutations_by_seq, n_size):
    n_mutations = random_fit_nonnegative(mutations_by_seq, len(mutations_by_seq))
    info_list = []
    for j, k in enumerate(n_mutations):
        mutations = generate_mutations(ref, m_info.mutation_probability, k)
        for v in get_1d_clusters(sorted(mutations)):
            tandem = ''.join(mutations[m] for m in v)
            t_size = len(v)

            mutated_subseq = re.sub('[.]+', tandem, get_context(ref, v[0], 2, t_size))
            # check if a stop codon is produced
            stop = any(s in mutated_subseq for s in STOP_CODONS)

            info_list.append([sim_id + 1, k, str(j + 1).zfill(n_size), v[0], ref[v[0]: v[0] + t_size].seq, tandem, t_size, stop])

    return info_list


def generate_mutations(ref, mutation_probabilities, n):
    """Generates a set of mutations on a sequence based on a vector of probabilities."""
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
        # randomly select one of the mutation candidates and generate a mutation
        pos = int(np.random.choice(mutation_candidates))
        mutated_base = mutate(ref[pos])
        mutations[pos] = mutated_base

    return mutations


def mutate(base):
    """Chooses randomly a different base than the one from the argument."""
    return np.random.choice([x for x in "ACGT" if x != base])


def get_context(seq, pos, n=1, s=1):
    """
    Gets context from a sequence at a given position.
    :param seq: the sequence
    :param pos: target position
    :param n: context size (each side)
    :param s: size of the target region, displaces the position of right context
    :return: a string matching [ACTG]*[.]+[ACTG]*
    """
    context = []
    for n in range(pos-n, pos+n+s):
        if pos <= n < pos+s:
            # target region, will be represented with dots
            context.append('.')
        elif 0 <= n < len(seq):
            # context region, if it is in sequence boundaries
            context.append(seq[n])
        # else (positions outside sequence) => nothing

    return ''.join(context)


def get_1d_clusters(data, stepsize=1):
    """
    Finds clusters of elements in an array. Returns only clusters with more than 1 element.
    :param data: input array
    :param stepsize: minimun distance between elements to cluster them
    :return: a list of arrays with the clusters found
    """
    consecutive = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    return [x for x in consecutive if len(x) > 1]


def random_fit_nonnegative(values, n):
    """
    Generates n random values using a normal distribution fitted from values used as argument.
    Returns only non-negative values.
    :param values: array/list to use as model fot the random data
    :param n: number of random elements to return
    :returns: an array of n random non-negative numbers
    """
    values = np.array(values)
    mean = np.mean(values)
    sd = np.std(values)

    random_values = np.empty(0)
    offset = 0.05  # 5% offset to compa values less than 0
    while len(random_values) < n:
        random_values = np.round(np.random.normal(mean, sd, round(n*(1+offset))))
        random_values = random_values[random_values >= 0]
        offset *= 2  # If the while loop check fail, next time will try with a larger offset

    # slice n first elements and shape the array to int
    return random_values[:n].astype('int')


if __name__ == '__main__':
    main()
