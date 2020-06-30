#!/usr/bin/env python3
import argparse
import csv
import re
from multiprocessing import Pool
from random import choices

import numpy as np
import pandas as pd

from Bio import Seq, SeqIO


def main():
    # read arguments
    args = parse_arguments()
    target_name = re.sub(r"^.*/|[.].*$", "", args.mutation_info)

    # load data
    references = {x.id: x for x in SeqIO.parse(args.ref_file, "fasta")}
    with open(args.mutations_by_seq) as f:
        for line in f.read().splitlines():
            line = line.split("\t")
            if line[0] == target_name:
                mutations_by_seq = [int(x) for x in line[1].split(",")]
                break

    m_info = pd.read_csv(args.mutation_info, sep="\t")
    ref = (
        [references[target_name]]
        if not args.referencenames
        else [references[x] for x in args.referencenames]
    )
    n_size = len(str(args.n))

    # if references are multiple, then make a modified m_info for each one.
    if len(ref) > 1:
        m_info_by_ref = [correct_mutation_info(m_info, x.seq) for x in ref]

    # runs in parallel n simulations with len(mutations_by_seq) sequences
    with Pool(args.threads) as p:
        tandem_info = p.starmap(
            run_simulation,
            (
                (
                    i,
                    m_info_by_ref[r] if len(ref) > 1 else m_info,
                    ref[r],
                    mutations_by_seq,
                    n_size,
                    args.method,
                    args.productiveonly,
                    args.fitnmutations,
                )
                for i, r in zip(range(args.n), choices(range(len(ref)), k=args.n))
            ),
        )

    df_columns = [
        "simulation",
        "mutations_byseq",
        "seq_id",
        "pos",
        "ref",
        "alt",
        "size",
        "context",
        "frame_pos",
        "aa_change",
        "stop_codon",
    ]

    # writes csv directly. Sorting and compression goes outside.
    with open(args.outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(df_columns)  # header

        while tandem_info:
            writer.writerows(tandem_info.pop())


# "main" operation wrapped in a function, to being able to paralelize it
def run_simulation(
    sim_id,
    m_info,
    ref,
    mutations_by_seq,
    n_size,
    method,
    productive_only=False,
    fit_normal_nmutations=False,
):
    mutations_fn = (
        generate_mutations_precandidating
        if method == "precandidating"
        else generate_mutations_sampling
    )
    n_mutations = (
        random_fit_nonnegative(mutations_by_seq, len(mutations_by_seq))
        if fit_normal_nmutations
        else mutations_by_seq
    )

    # converting data before the loop to save some time
    mutation_probability = np.array(m_info.mutation_probability)
    nucl_probabilities = m_info[["A", "C", "G", "T"]].to_dict(orient="index")

    info_list = []
    for i, n in enumerate(n_mutations):
        stop_codon = True
        while stop_codon:
            info_sublist = []
            mutations = mutations_fn(ref, mutation_probability, n, nucl_probabilities)
            for v in get_1d_clusters(sorted(mutations)):
                tandem = "".join(mutations[m] for m in v)
                t_size = len(v)

                context = get_context(ref, v[0], 2, t_size)
                # create the "mutated subsequence" replacing the tandem alt into the context
                mutated_subseq = re.sub(r"[.]+", tandem, context)
                # check if a stop codon is produced
                stop_codon = has_stop_codons(mutated_subseq, max(v[0] - 2, 0))
                if stop_codon and productive_only:
                    break

                info_sublist.append(
                    (
                        sim_id + 1,
                        n,
                        str(i + 1).zfill(n_size),
                        v[0] + 1,
                        str(ref[v[0] : v[0] + t_size].seq),
                        tandem,
                        t_size,
                        context,
                        (v[0] % 3) + 1,
                        get_aa_change(ref, mutated_subseq, max(v[0] - 2, 0)),
                        stop_codon,
                    )
                )

            stop_codon = False
        info_list += info_sublist

    return info_list


def generate_mutations_precandidating(
    ref, mutation_probabilities, n, nucl_probabilities
):
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
        mutations[pos] = mutate(ref[pos], nucl_probabilities[pos])

    return mutations


def generate_mutations_sampling(ref, mutation_probabilities, n, nucl_probabilities):
    """Generates a set of mutations on a sequence based on a vector of probabilities, using a sampling method"""
    seq_len = len(ref)
    positions = np.arange(0, seq_len)
    mutation_prob = np.array(mutation_probabilities)
    mutation_prob /= np.sum(mutation_prob)  # array sum must me 1
    mutations = {
        int(m): mutate(ref[int(m)], nucl_probabilities[int(m)])
        for m in np.random.choice(positions, size=n, replace=False, p=mutation_prob)
    }

    return mutations


def mutate(base, nucl_probability=None):
    """Chooses randomly a different base than the one from the argument or one based on a vector of probabilities."""
    if nucl_probability:
        return np.random.choice(
            ("A", "C", "G", "T"), p=[nucl_probability[x] for x in "ACGT"]
        )
    else:
        return np.random.choice([x for x in "ACGT" if x != base])


def has_stop_codons(seq, pos):
    """
    Checks stop codons.
    :param seq: mutated subsequence
    :param pos: position of the subsequence in terms of the original full sequence
    :return: boolean value indicating if the sequence has stop codons or not.
    """
    stop_codons = ("TAA", "TAG", "TGA")
    # get codons using reading frame based on pos
    codons = re.findall(".{3}", seq[(3 - pos % 3) % 3 :])
    return any(s in codons for s in stop_codons)


def get_aa_change(ref, seq, pos):
    """
    Builds aa change using reference and mutated sequence.
    :param ref: full reference sequence (Bio.Seq)
    :param seq: mutated subsequence (mutation + 2nt context each side)
    :param pos: position of the subsequence in terms of the original full sequence
    :return: a string in the form "{}>{}" with the aa changes
    """
    # position of the first full frame in the subsequence
    framed_pos = int((pos + 2) / 3) * 3
    coding_alt = seq[framed_pos - pos :]
    coding_alt = Seq.Seq(coding_alt[: 3 * int(len(coding_alt) / 3)])
    coding_ref = ref[framed_pos : framed_pos + len(coding_alt)]

    return f"{str(coding_ref.seq.translate())}>{str(coding_alt.translate())}"


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
            context.append(".")
        elif 0 <= n < len(seq):
            # context region, if it is in sequence boundaries
            context.append(seq[n])
        # else (positions outside sequence) => nothing

    return "".join(context)


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
        # If the while loop check fail, next time will try with a larger offset
        offset *= 2

    # slice n first elements and shape the array to int
    return random_values[:n].astype("int")


def correct_mutation_info(m_info, reference):
    """
    Makes a modified version of the mutation_info dataframe using the reference to correct it.
    Used when a single mutation_info object is used against multiple references.
    """
    ref = str(reference)
    m_info = m_info.copy()

    # iterate over sequence positions
    for b, r in zip(ref, m_info.itertuples()):
        # if the probability of the reference base is not 0 a correction is needed
        if b in "ACGT" and getattr(r, b) != 0:
            m_info.loc[m_info.position == r.position, b] = 0
            if getattr(r, b) == 1:
                # If the probability was 1, is enough with this
                m_info.loc[m_info.position == r.position, "mutation_probability"] = 0
            else:
                # If not was 1, then the other probabilities have to be corrected as well
                bases = [x for x in "ACGT" if x != b]
                psum = sum(getattr(r, x) for x in bases)
                m_info.loc[m_info.position == r.position, bases] = [
                    getattr(r, x) / psum for x in bases
                ]

    return m_info


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Simulates tandem generation based on dataset"
    )

    parser.add_argument(
        "mutation_info",
        metavar="mutation_info_file",
        help="file with mutation probabilities (.tsv)",
    )
    parser.add_argument(
        "ref_file", metavar="reference_seq_file", help="references file (.fasta)"
    )
    parser.add_argument(
        "mutations_by_seq",
        metavar="mutations_by_seq_file",
        help="file with mutations by each sequence in dataset",
    )
    parser.add_argument("n", type=int, metavar="n", help="number of simulations")
    parser.add_argument("outfile", metavar="output_file", help="output file (.tsv)")
    parser.add_argument(
        "--method",
        "-m",
        metavar="mutation_method",
        choices=["precandidating", "sampling"],
        default="sampling",
        help="method used to generate mutations",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="number of threads to use in the simulation",
    )
    parser.add_argument(
        "--referencenames",
        "-r",
        action="store",
        nargs="+",
        metavar="REF",
        help="reference sequences to make the simulation, if not provided then mutation_info filename will be used",
    )
    parser.add_argument(
        "--productiveonly",
        "-p",
        action="store_true",
        help="do not generate stop codons",
    )
    parser.add_argument(
        "--fitnmutations",
        "-f",
        action="store_true",
        help="get number of mutations by sequence from a normal distribution fitted from real data",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
