#!/usr/bin/env python3
import re
import argparse
from multiprocessing import Pool
from itertools import chain
from Bio import SeqIO
import pandas as pd
import numpy as np


def main():
    # read arguments
    args = parse_arguments()
    n = int(args.n)
    target_name = re.sub('^.*/|[.].*$', '', args.mutation_info)

    # load data
    references = {x.id: x for x in SeqIO.parse(args.ref_file, 'fasta')}
    with open(args.mutations_by_seq) as f:
        for line in f.read().splitlines():
            line = line.split('\t')
            if line[0] == target_name:
                mutations_by_seq = [int(x) for x in line[1].split(',')]
                break

    m_info = pd.read_table(args.mutation_info)
    ref = references[target_name]
    n_size = len(str(n))

    # runs in parallel n simulations with len(mutations_by_seq) sequences
    with Pool() as p:
        tandem_info = p.starmap(run_simulation, ((i, m_info, ref, mutations_by_seq, n_size, args.method, args.productiveonly, args.fitnmutations) for i in range(n)))

    # flatten list of lists to a single list
    tandem_info = list(chain.from_iterable(tandem_info))

    tandem_df = pd.DataFrame(tandem_info, columns=['simulation', 'mutations_byseq', 'seq_id', 'pos', 'ref', 'alt', 'size', 'stop_codon'])
    tandem_df.to_csv(args.outfile + '.gz', compression='gzip', sep='\t', index=False)


# "main" operation wrapped in a function, to being able to paralelize it
def run_simulation(sim_id, m_info, ref, mutations_by_seq, n_size, method, productive_only=False, fit_normal_nmutations=False):
    mutations_fn = generate_mutations_precandidating if method == 'precandidating' else generate_mutations_sampling
    n_mutations = random_fit_nonnegative(mutations_by_seq, len(mutations_by_seq)) if fit_normal_nmutations else mutations_by_seq
    info_list = []
    for j, k in enumerate(n_mutations):
        stop_codon = True
        while stop_codon:
            info_sublist = []
            mutations = mutations_fn(ref, m_info.mutation_probability, k)
            for v in get_1d_clusters(sorted(mutations)):
                tandem = ''.join(mutations[m] for m in v)
                t_size = len(v)

                mutated_subseq = re.sub('[.]+', tandem, get_context(ref, v[0], 2, t_size))
                # check if a stop codon is produced
                stop_codon = has_stop_codons(mutated_subseq, max(v[0] - 2, 0))
                if stop_codon and productive_only:
                    break

                info_sublist.append([sim_id + 1, k, str(j + 1).zfill(n_size), v[0], ref[v[0]: v[0] + t_size].seq, tandem, t_size, stop_codon])
            stop_codon = False
        info_list += info_sublist

    return info_list


def generate_mutations_precandidating(ref, mutation_probabilities, n):
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
        mutations[pos] = mutate(ref[pos])

    return mutations


def generate_mutations_sampling(ref, mutation_probabilities, n):
    """Generates a set of mutations on a sequence based on a vector of probabilities, using a sampling method"""
    seq_len = len(ref)
    positions = np.arange(0, seq_len)
    mutation_prob = np.array(mutation_probabilities)
    mutation_prob /= np.sum(mutation_prob)  # array sum must me 1
    mutations = {int(m): mutate(ref[int(m)]) for m in np.random.choice(positions, size=n, replace=False, p=mutation_prob)}

    return mutations


def mutate(base):
    """Chooses randomly a different base than the one from the argument."""
    return np.random.choice([x for x in "ACGT" if x != base])


def has_stop_codons(seq, pos):
    """Checks stop codons"""
    stop_codons = ('TAA', 'TAG', 'TGA')
    # get codons using reading frame based on pos
    codons = re.findall('.{3}', seq[pos % 3:])
    return any(s in codons for s in stop_codons)


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
    for n in range(pos - n, pos + n + s):
        if pos <= n < pos + s:
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
    consecutive = np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)
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
    offset = 0.05  # 5% offset to compensate values less than 0
    while len(random_values) < n:
        random_values = np.round(np.random.normal(mean, sd, round(n * (1 + offset))))
        random_values = random_values[random_values >= 0]
        offset *= 2  # If the while loop check fail, next time will try with a larger offset

    # slice n first elements and shape the array to int
    return random_values[:n].astype('int')


def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulates tandem generation based on dataset")

    parser.add_argument('mutation_info', metavar='mutation_info_file', help="file with mutation probabilities (.tsv)")
    parser.add_argument('ref_file', metavar='reference_seq_file', help="references file (.fasta)")
    parser.add_argument('mutations_by_seq', metavar='mutations_by_seq_file', help="file with mutations by each sequence in dataset")
    parser.add_argument('n', metavar='n', help="number of simulations")
    parser.add_argument('outfile', metavar='output_file', help="output file (.tsv)")
    parser.add_argument('--method', '-m', metavar='mutation_method', choices=['precandidating', 'sampling'], default='precandidating',
                        help="method used to generate mutations")
    parser.add_argument('--productiveonly', '-p', action='store_true', help="do not generate stop codons")
    parser.add_argument('--fitnmutations', '-f', action='store_true',
                        help="get number of mutations by sequence from a normal distribution fitted from real data")

    return parser.parse_args()


if __name__ == '__main__':
    main()
