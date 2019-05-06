# This module contains plots to show data
# Assummes the files are located in certains places
from itertools import product
from os import path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats

# Not the best idea, but for now...
TANDEM_REAL_DIR = "tandem_info/real"
TANDEM_SIM_SUMM_DIR = "tandem_info_summarized/sampling"


def tandem_distplot(family_name, show_pdf=True, split_axis=False):
    """Plots the distribution of the tandems by simulation vs the real value"""
    # load data
    df_real = pd.read_table(
        _find_filename("{}/{}.tsv".format(TANDEM_REAL_DIR, family_name))
    )
    df_sim = pd.read_table(
        _find_filename(
            "{}/{}/count_sim_size.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
        )
    )

    # filter size
    df_real = df_real[df_real["size"] == 2]
    df_sim = df_sim[df_sim["size"] == 2]

    real_tandems = len(df_real)
    sim_tandems = df_sim.groupby("simulation").agg({"n": "sum"})

    # this value determines bin size and also where to cut axis if split_axis is True
    adjust_value = 50 if max(sim_tandems.n) - min(sim_tandems.n) < 150 else 100

    bin_step = np.ceil(
        real_tandems / (adjust_value * 5)
    )  # guarantees good visualization in the pdf output
    if split_axis:
        # Make a figure with 2 subplots
        fig, (ax, ax2) = plt.subplots(
            1, 2, sharey=True, figsize=(12, 8), gridspec_kw={"width_ratios": [5, 1]}
        )
        ax.set_xlim(
            adjust_value * np.floor(min(sim_tandems.n) / adjust_value),
            adjust_value * np.ceil(max(sim_tandems.n) / adjust_value),
        )
        ax2_limits = [
            adjust_value * np.floor(real_tandems / adjust_value),
            adjust_value * np.ceil(real_tandems / adjust_value),
        ]
        ax2.set_xlim(*ax2_limits)
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
        ax2 = (
            ax
        )  # this way is not neccesary to use conditions to know where to plot the real value

    # plot simulated and real values
    ax.hist(
        sim_tandems.n,
        bins=np.arange(min(sim_tandems.n), max(sim_tandems.n), bin_step),
        align="left",
        lw=0.1,
        label="simulated",
    )
    ax2.axvline(
        x=real_tandems,
        ymin=0,
        ymax=1,
        linewidth=1.5,
        color="forestgreen",
        label="observed",
    )

    # annotations and labels
    ax2.text(
        real_tandems,
        ax.get_ylim()[1] * 0.9,
        "observed = {}".format(real_tandems),
        ha="center",
        bbox=dict(fc="w", ec="0.5", alpha=0.5),
    )
    ax.set(ylabel="Simulations")
    fig.suptitle(
        "Distribution of tandems by simulation ({})".format(family_name),
        fontsize=rcParams["axes.titlesize"],
    )
    fig.text(0.5, -0.01, "Tandems by simulations", ha="center")
    # ax.legend(loc='upper left')

    if show_pdf:
        pdf = stats.norm.pdf(
            real_tandems, np.mean(sim_tandems.n), np.std(sim_tandems.n)
        )
        ax2.text(
            real_tandems,
            ax.get_ylim()[1] * 0.03,
            "pdf = {:.2e}".format(pdf),
            ha="right",
            bbox=dict(fc="w", ec="0.5", alpha=0.7),
        )

    sns.despine()
    if split_axis:
        # options to make the two plots look like a single one
        ax.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax.tick_params(labelright=False)
        ax2.xaxis.set_ticks(ax2_limits)
        ax2.yaxis.set_ticks_position("none")

    return fig


def mutations_byseq_distplot(family_name):
    """Plots the distribution of the mutations by each sequence"""
    with open("mutation_info/mutations_by_seq.txt") as f:
        for line in f.read().splitlines():
            line = line.split("\t")
            if line[0] == family_name:
                mutations_by_seq = [int(x) for x in line[1].split(",")]
                break

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(
        mutations_by_seq,
        bins=np.arange(min(mutations_by_seq), max(mutations_by_seq)),
        align="left",
    )
    ax.set(
        title="Distribution of mutations by sequence ({})".format(family_name),
        xlabel="Mutations by sequence",
        ylabel="Sequences",
    )
    sns.despine()

    return fig


def mutation_probability(family_name):
    """Plots the probability of mutation of each position"""
    df = pd.read_table(_find_filename("mutation_info/{}.tsv".format(family_name)))

    # plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x=df.position, height=df.mutation_probability, lw=0)
    ax.set(
        title="Mutation probabilities by position ({})".format(family_name),
        xlabel="Position in sequence",
        ylabel="Mutation probability",
    )
    sns.despine()

    return fig


def cluster_size_distribution(family_name):
    """Plots the distribution of the sizes of clusters found (simulated vs real)"""
    df_real = pd.read_table(
        _find_filename("{}/{}.tsv".format(TANDEM_REAL_DIR, family_name))
    )
    df_sim = pd.read_table(
        _find_filename(
            "{}/{}/count_sim_size.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
        )
    )

    # group clusters by size
    df_real_bysize = (
        df_real.groupby("size")
        .agg({"seq_id": "count"})
        .reset_index()
        .rename(columns={"seq_id": "real"})
    )
    df_sim_bysize = (
        df_sim.groupby("size")
        .agg({"n": "sum"})
        .reset_index()
        .rename(columns={"n": "simulation"})
    )
    df_sim_bysize["simulation"] = np.round(
        df_sim_bysize["simulation"] / len(df_sim.simulation.unique())
    )  # get the average of all the simulations
    df_sim_bysize = df_sim_bysize[
        df_sim_bysize.simulation != 0
    ]  # drop the rows with 0 tandems

    # join data in a single dataframe
    df_by_size = pd.merge(df_real_bysize, df_sim_bysize, how="outer").fillna(0)
    df_by_size = df_by_size.melt(
        id_vars="size",
        value_vars=["real", "simulation"],
        var_name="target",
        value_name="n",
    )

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=df_by_size, x="size", y="n", hue="target", ax=ax)
    # annotations
    for p in ax.patches:
        ax.annotate(
            int(p.get_height()),
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="bottom",
            rotation=90,
            xytext=(0, 5),
            color="gray",
            textcoords="offset points",
        )
    ax.set_ylim([0, max(df_by_size.n) * 1.2])
    ax.set(
        title="Distribution of mutation clusters by size ({})".format(family_name),
        xlabel="Cluster size",
        ylabel="Number of clusters",
    )
    ax.legend()
    sns.despine()

    return fig


def stop_codon_distribution(family_name):
    """Plots the distribution of the clusters found (with/without stop codons)"""
    df_nostop = pd.read_table(
        _find_filename(
            "{}/{}/count_sim_size.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
        )
    )
    df_stop = pd.read_table(
        _find_filename(
            "{}/{}/count_sim_size_stop.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
        )
    )

    # group clusters by size
    df_nostop_bysize = (
        df_nostop.groupby("size")
        .agg({"n": "sum"})
        .reset_index()
        .rename(columns={"n": "no stop codon"})
    )
    df_nostop_bysize["no stop codon"] = np.round(
        df_nostop_bysize["no stop codon"] / len(df_nostop.simulation.unique())
    )
    df_nostop_bysize = df_nostop_bysize[df_nostop_bysize["no stop codon"] != 0]

    df_stop_bysize = (
        df_stop.groupby("size")
        .agg({"n": "sum"})
        .reset_index()
        .rename(columns={"n": "stop codon"})
    )
    df_stop_bysize["stop codon"] = np.round(
        df_stop_bysize["stop codon"] / len(df_stop.simulation.unique())
    )
    df_stop_bysize = df_stop_bysize[df_stop_bysize["stop codon"] != 0]

    # join data in a single dataframe
    df_by_size = pd.merge(df_nostop_bysize, df_stop_bysize, how="outer").fillna(0)
    df_by_size = df_by_size.melt(
        id_vars="size",
        value_vars=["no stop codon", "stop codon"],
        var_name="category",
        value_name="n",
    )

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=df_by_size, x="size", y="n", hue="category", ax=ax)
    # annotations
    for p in ax.patches:
        ax.annotate(
            int(p.get_height()),
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="bottom",
            rotation=90,
            xytext=(0, 5),
            color="gray",
            textcoords="offset points",
        )
    ax.set_ylim([0, max(df_by_size.n) * 1.2])
    ax.set(
        title="Distribution of mutation clusters with stop codons ({})".format(
            family_name
        ),
        xlabel="Cluster size",
        ylabel="Number of clusters",
    )
    ax.legend()
    sns.despine()

    return fig


def tandem_heatmap(family_name, target="real", stop_codons=False):
    """Plots a heatmap of all the possible tandem substitutions"""
    if target == "real":
        df = pd.read_table(
            _find_filename("{}/{}.tsv".format(TANDEM_REAL_DIR, family_name))
        )
    elif target == "simulated" and stop_codons is False:  # simulated
        df = pd.read_table(
            _find_filename(
                "{}/{}/count_sim_mut.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
            )
        )
    else:  # stop codons
        df = pd.read_table(
            _find_filename(
                "{}/{}/count_sim_mut_stop.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
            )
        )
    df = df[df["size"] == 2]

    # get all combinations of dinucleotides to fill with 0 later
    dinucleotides = ["".join(x) for x in product("ACGT", repeat=2)]
    dinucleotide_combinations = [
        [x, y, 0]
        for x in dinucleotides
        for y in dinucleotides
        if (x[0] != y[0] and x[1] != y[1])
    ]
    dc_df = pd.DataFrame(dinucleotide_combinations, columns=["ref", "alt", "n"])

    # group by substitution and count
    if target == "real":
        grouped = (
            df.groupby(["ref", "alt"])
            .agg({"seq_id": "count"})
            .reset_index()
            .rename(columns={"seq_id": "n"})
        )
    else:
        grouped = df.groupby(["ref", "alt"]).agg({"n": "sum"}).reset_index()
    # add rows with 0 and make the table
    table = (
        pd.concat([grouped, dc_df])
        .drop_duplicates(subset=["ref", "alt"], keep="first")
        .pivot("ref", "alt", "n")
    )

    # if the dataset is from the simulation, get the average between all the simulations
    if target == "simulated":
        table = np.round(table / len(df.simulation.unique()), 1)

    # some variables
    annot_fmt = ".0f" if target == "real" else ".1f"
    cbar_label = "Number of tandems" + (
        " (average of simulations)" if target == "simulated" else ""
    )
    stop_msg = "with stop codons" if stop_codons else ""

    # plot
    with sns.axes_style("darkgrid"):
        fig, ax = plt.subplots(figsize=(13, 11))
        sns.heatmap(
            table.iloc[::-1],
            ax=ax,
            cmap="YlGnBu",
            annot=True,
            fmt=annot_fmt,
            square=True,
            linewidths=1,
            cbar_kws={"label": cbar_label},
        )
        ax.set(
            title="Tandem distribution {} from {} data ({})".format(
                stop_msg, target, family_name
            ),
            xlabel="ALT",
            ylabel="REF",
        )

    return fig


def tandem_heatmap_percentage(family_name):
    """Plots a heatmap of all the possible tandem substitutions"""
    df_real = pd.read_table(
        _find_filename("{}/{}.tsv".format(TANDEM_REAL_DIR, family_name))
    )
    df_real = df_real[df_real["size"] == 2]

    df_sim = pd.read_table(
        _find_filename(
            "{}/{}/count_sim_mut.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
        )
    )
    df_sim = df_sim[df_sim["size"] == 2]

    # get all combinations of dinucleotides to fill with 0 later
    dinucleotides = ["".join(x) for x in product("ACGT", repeat=2)]
    dinucleotide_combinations = [
        [x, y, 0]
        for x in dinucleotides
        for y in dinucleotides
        if (x[0] != y[0] and x[1] != y[1])
    ]
    dc_df = pd.DataFrame(dinucleotide_combinations, columns=["ref", "alt", "n"])

    # group by substitution and count
    grouped_real = (
        df_real.groupby(["ref", "alt"])
        .agg({"seq_id": "count"})
        .reset_index()
        .rename(columns={"seq_id": "n"})
    )
    grouped_sim = df_sim.groupby(["ref", "alt"]).agg({"n": "sum"}).reset_index()

    # add rows with 0 and make the tables
    table_real = (
        pd.concat([grouped_real, dc_df])
        .drop_duplicates(subset=["ref", "alt"], keep="first")
        .pivot("ref", "alt", "n")
    )
    table_sim = (
        pd.concat([grouped_sim, dc_df])
        .drop_duplicates(subset=["ref", "alt"], keep="first")
        .pivot("ref", "alt", "n")
    )

    # get the average between all the simulations
    table_sim = table_sim / len(df_sim.simulation.unique())

    # get percentage table
    table = np.round(table_sim / table_real * 100, 0)
    table = table.replace([np.inf, -np.inf], np.nan)

    # plot
    with sns.axes_style("darkgrid"):
        fig, ax = plt.subplots(figsize=(13, 11))
        sns.heatmap(
            table.iloc[::-1],
            ax=ax,
            cmap="YlGnBu",
            annot=True,
            square=True,
            fmt=".0f",
            linewidths=1,
            cbar_kws={"label": "Percentage of false tandems (simulated/real)"},
        )
        ax.set(
            title="Percentage of false Tandems ({})".format(family_name),
            xlabel="ALT",
            ylabel="REF",
        )

    return fig


def tandem_mutation_probability_heatmap(family_name, target="real", base2=True):
    """"Groups tandems by position mutation probability and reference in each nucleotide of the tandem"""
    if target == "real":
        df = pd.read_table(
            _find_filename("{}/{}.tsv".format(TANDEM_REAL_DIR, family_name))
        )
        df = df[df["size"] == 2]
    else:
        df = pd.read_table(
            _find_filename(
                "{}/{}/count_sim_pos.tsv".format(TANDEM_SIM_SUMM_DIR, family_name)
            )
        )
        df = df[df["size"] == 2]
        n_simulations = len(df.simulation.unique())
        df = df.groupby(["pos", "ref", "size"]).agg({"n": "sum"}).reset_index()

    df_prob = pd.read_table("mutation_info/{}.tsv".format(family_name))

    df["2"] = df.pos + 1  # 2nd nucleotide
    df = df.reset_index().rename(columns={"pos": "1", "index": "id"})
    df[["ref1", "ref2"]] = (
        df["ref"].str.extractall("(.)")[0].unstack()
    )  # split reference into 2

    # prepare 2 long-format dataframes with sequence position and reference base
    dfl_pos = pd.melt(
        df,
        id_vars=["id"],
        value_vars=["1", "2"],
        var_name="tposition",
        value_name="position",
    )
    dfl_ref = pd.melt(
        df[["id", "ref1", "ref2"]].rename(columns={"ref1": "1", "ref2": "2"}),
        id_vars=["id"],
        value_vars=["1", "2"],
        var_name="tposition",
        value_name="ref",
    )

    # merge all
    dfl = pd.merge(
        pd.merge(dfl_pos, dfl_ref), df_prob[["position", "mutation_probability"]]
    )

    # transform dataframe to put ref and mutation probability into separated columns for each nucleotide in tandem
    dfw = dfl.pivot(
        index="id", columns="tposition", values=["mutation_probability", "ref"]
    ).reset_index()
    dfw.columns = ["id", "mp1", "mp2", "ref1", "ref2"]

    # add "n"
    if target == "real":
        dfw["n"] = 1
    else:
        dfw = pd.merge(dfw, df[["id", "n"]])

    # group data by reference nucleotides and mutation probability intervals
    intervals = (
        np.insert(np.logspace(-6, 0, 7, base=2), 0, 0)
        if base2
        else np.arange(0, 1, 0.02)
    )
    dfw_gr = (
        dfw.groupby(
            [
                "ref1",
                "ref2",
                pd.cut(dfw["mp1"], intervals),
                pd.cut(dfw["mp2"], intervals),
            ]
        )
        .agg({"n": "sum"})
        .reset_index()
    )

    # combine mutation probability and ref in a single field to group later
    fmt_fn = _base2_str if base2 else lambda x: x
    dfw_gr["nucl1"] = dfw_gr.apply(
        lambda x: r"$({},{}]$ / {}".format(
            fmt_fn(x.mp1.left), fmt_fn(x.mp1.right), x.ref1
        ),
        axis=1,
    )
    dfw_gr["nucl2"] = dfw_gr.apply(
        lambda x: r"$({},{}]$ / {}".format(
            fmt_fn(x.mp2.left), fmt_fn(x.mp2.right), x.ref2
        ),
        axis=1,
    )

    # make the table
    t = dfw_gr.pivot(index="nucl1", columns="nucl2", values="n")
    # sort it correctly
    names_sorted = sorted(
        t.columns, key=lambda x: np.exp2(int(x[5:7])) if x[2] != "0" else 0
    )
    t = t.reindex(names_sorted, axis=1).reindex(names_sorted, axis=0)

    if target == "simulated":
        t = t / n_simulations
    t = t.fillna(0)

    # plot
    with sns.axes_style("darkgrid"):
        fig, ax = plt.subplots(figsize=(13, 11))
        sns.heatmap(
            t.iloc[::-1],
            ax=ax,
            cmap="YlGnBu",
            square=True,
            linewidths=1,
            cbar_kws={"label": "Number of tandems"},
        )
        ax.set(
            title="Tandem distribution by mutation probability and reference of {} data ({})".format(
                target, family_name
            ),
            xlabel="Mutation probability / Reference nucleotide nucleotide 2",
            ylabel="Mutation probability / Reference nucleotide nucleotide 1",
        )

    return fig


def _base2_str(x):
    return r"2^{" + str(int(round(np.log2(x), 0))) + "}" if x > 0 else "0"


def _find_filename(filepath):
    for suffix in ["", ".gz", ".xz"]:
        if path.isfile(filepath + suffix):
            return filepath + suffix
    raise FileNotFoundError("File {} does not exist".format(filepath))
