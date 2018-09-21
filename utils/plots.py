# This module contains plots to show data
# Assummes the files are located in certains places
from itertools import product
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats


def tandem_distplot(family_name):
    """Plots the distribution of the tandems by simulation vs the real value"""
    # load data
    df_real = pd.read_table("tandem_info/real/{}.tsv".format(family_name))
    df_sim = pd.read_table("tandem_info_summarized/{}/count_sim_size.tsv".format(family_name))

    # filter size
    df_real = df_real[df_real['size'] == 2]
    df_sim = df_sim[df_sim['size'] == 2]

    real_tandems = len(df_real)
    sim_tandems = df_sim.groupby('simulation').agg({'n': 'sum'})

    pdf = stats.norm.pdf(real_tandems, np.mean(sim_tandems.n), np.std(sim_tandems.n))

    # plot
    bin_step = np.ceil(real_tandems / 500)  # guarantees good visualization in the pdf output
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.hist(sim_tandems.n, bins=np.arange(min(sim_tandems.n), max(sim_tandems.n), bin_step), align='left', lw=0.1, label='simulated')
    ax.axvline(x=real_tandems, ymin=0, ymax=1, linewidth=1.5, color='forestgreen', label='real')
    ax.text(real_tandems, ax.get_ylim()[1] * 0.85, "real = {}".format(real_tandems), ha='center', bbox=dict(fc="w", ec="0.5", alpha=0.7))
    ax.text(real_tandems, ax.get_ylim()[1] * 0.03, "pdf = {:.2e}".format(pdf), ha='right', bbox=dict(fc="w", ec="0.5", alpha=0.7))
    ax.set(title="Distribution of tandems by simulation size=2 ({})".format(family_name), xlabel="Tandems by simulation", ylabel="Simulations")
    ax.legend(loc='upper center')
    sns.despine()

    return fig


def mutations_byseq_distplot(family_name):
    """Plots the distribution of the mutations by each sequence"""
    with open("mutation_info/mutations_by_seq.txt") as f:
        for line in f.read().splitlines():
            line = line.split('\t')
            if line[0] == family_name:
                mutations_by_seq = [int(x) for x in line[1].split(',')]
                break

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(mutations_by_seq, bins=np.arange(min(mutations_by_seq), max(mutations_by_seq)), align='left')
    ax.set(title="Distribution of mutations by sequence ({})".format(family_name), xlabel="Mutations by sequence", ylabel="Sequences")
    sns.despine()

    return fig


def mutation_probability(family_name):
    """Plots the probability of mutation of each position"""
    df = pd.read_table("mutation_info/{}.tsv".format(family_name))

    # plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x=df.position, height=df.mutation_probability, lw=0)
    ax.set(title="Mutation probabilities by position ({})".format(family_name), xlabel="Position in sequence", ylabel="Mutation probability")
    sns.despine()

    return fig


def cluster_size_distribution(family_name):
    """Plots the distribution of the sizes of clusters found (simulated vs real)"""
    df_real = pd.read_table("tandem_info/real/{}.tsv".format(family_name))
    df_sim = pd.read_table("tandem_info_summarized/{}/count_sim_size.tsv".format(family_name))

    # group clusters by size
    df_real_bysize = df_real.groupby('size').agg({'seq_id': 'count'}).reset_index().rename(columns={'seq_id': 'n'})
    df_sim_bysize = df_sim.groupby('size').agg({'n': 'sum'}).reset_index()
    df_sim_bysize['n'] = np.round(df_sim_bysize['n'] / 100000)  # get the average of the 100000 simulations

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(df_real_bysize['size'], df_real_bysize['n'], label='real', alpha=0.5)
    ax.bar(df_sim_bysize['size'], df_sim_bysize['n'], label='simulation', alpha=0.5)
    ax.set(title="Distributions of mutation clusters by size ({})".format(family_name), xlabel="Cluster size", ylabel="Number of clusters")
    ax.legend()
    sns.despine()

    return fig


def tandem_heatmap(family_name, target="real"):
    """Plots a heatmap of all the possible tandem substitutions"""
    if target == "real":
        df = pd.read_table("tandem_info/real/{}.tsv".format(family_name))
    else:  # simulated
        df = pd.read_table("tandem_info_summarized/{}/count_sim_mut.tsv".format(family_name))
    df = df[df['size'] == 2]

    # get all combinations of dinucleotides to fill with 0 later
    dinucleotides = [''.join(x) for x in product("ACGT", repeat=2)]
    dinucleotide_combinations = [[x, y, 0] for x in dinucleotides for y in dinucleotides if (x[0] != y[0] and x[1] != y[1])]
    dc_df = pd.DataFrame(dinucleotide_combinations, columns=['ref', 'alt', 'n'])

    # group by substitution and count
    grouped = pd.DataFrame(df.groupby(['ref', 'alt']).size()).reset_index().rename(columns={0: 'n'})
    # add rows with 0 and make the table
    table = pd.concat([grouped, dc_df]).drop_duplicates(subset=['ref', 'alt'], keep='first').pivot('ref', 'alt', 'n')

    # if the dataset is from the simulation, get the average between the 100000 simulations
    if target == "simulated":
        table /= 100000

    # some variables
    annot_fmt = ".0f" if target == "real" else ".2f"
    cbar_label = "Number of tandems" + (" (average of simulations)" if target == "simulated" else "")

    # plot
    with sns.axes_style("darkgrid"):
        fig, ax = plt.subplots(figsize=(13, 11))
        sns.heatmap(table, ax=ax, cmap="YlGnBu", annot=True, fmt=annot_fmt, square=True, linewidths=1, cbar_kws={'label': cbar_label})
        ax.set(title="Tandem distribution from {} data ({})".format(target, family_name), xlabel="ALT", ylabel="REF")

    return fig
