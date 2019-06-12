# This file contains plots to show data
from itertools import product

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt, rcParams
from scipy import stats


def tandem_distplot(
    real_tandems_fn, count_sim_size_fn, family_name, show_pdf=True, split_axis_auto=True
):
    """Plots the distribution of the tandems by simulation vs the real value"""
    # load data
    df_real = pd.read_csv(real_tandems_fn, sep="\t")
    df_sim = pd.read_csv(count_sim_size_fn, sep="\t")

    # filter size
    df_real = df_real[df_real["size"] == 2]
    df_sim = df_sim[df_sim["size"] == 2]

    real_tandems = len(df_real)
    sim_tandems = df_sim.groupby("simulation").agg({"n": "sum"})

    # this value determines bin size and also where to cut axis if split_axis is True
    n_sim_range = max(sim_tandems.n) - min(sim_tandems.n)
    adjust_value = 50 if n_sim_range < 150 else 100

    # guarantees good visualization in the pdf output
    bin_step = np.ceil(real_tandems / (adjust_value * 5))

    # split the axis in two if the real value is too far away from the simulation distribution
    split_axis = split_axis_auto and (real_tandems - max(sim_tandems.n)) / n_sim_range > 1
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
        # this way is not neccesary to use conditions to know where to plot the real value
        ax2 = ax

    # plot simulated and real values
    ax.hist(
        sim_tandems.n,
        bins=np.arange(min(sim_tandems.n), max(sim_tandems.n), bin_step),
        align="left",
        lw=0.1,
        edgecolor="w",
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
        y=0.92,
    )
    fig.text(0.5, 0.04, "Tandems by simulations", ha="center")

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


def mutations_byseq_distplot(mutations_by_seq_fn, family_name):
    """Plots the distribution of the mutations by each sequence"""
    with open(mutations_by_seq_fn) as f:
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
        rwidth=0.8
    )
    ax.set(
        title="Distribution of mutations by sequence ({})".format(family_name),
        xlabel="Mutations by sequence",
        ylabel="Sequences",
    )
    sns.despine()

    return fig


def mutation_probability(mutation_info_fn, family_name):
    """Plots the probability of mutation of each position"""
    df = pd.read_csv(mutation_info_fn, sep="\t")

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


def cluster_size_distribution(
    real_tandems_fn, count_size_fn, family_name, n_simulations
):
    """Plots the distribution of the sizes of clusters found (simulated vs real)"""
    df_real = pd.read_csv(real_tandems_fn, sep="\t")
    df_sim = pd.read_csv(count_size_fn, sep="\t")

    # group clusters by size
    df_real_bysize = (
        df_real.groupby("size")
        .agg({"seq_id": "count"})
        .reset_index()
        .rename(columns={"seq_id": "real"})
    )
    df_sim_bysize = df_sim.rename(columns={"n": "simulation"})
    # get the average of all the simulations
    df_sim_bysize["simulation"] = np.round(df_sim_bysize["simulation"] / n_simulations)
    # drop the rows with 0 tandems
    df_sim_bysize = df_sim_bysize[df_sim_bysize.simulation != 0]

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


def stop_codon_distribution(
    count_size_fn, count_size_stop_fn, family_name, n_simulations
):
    """Plots the distribution of the clusters found (with/without stop codons)"""
    df_nostop = pd.read_csv(count_size_fn, sep="\t")
    df_stop = pd.read_csv(count_size_stop_fn, sep="\t")

    # group clusters by size
    df_nostop_bysize = df_nostop.rename(columns={"n": "no stop codon"})
    df_nostop_bysize["no stop codon"] = np.round(
        df_nostop_bysize["no stop codon"] / n_simulations
    )
    df_nostop_bysize = df_nostop_bysize[df_nostop_bysize["no stop codon"] != 0]

    df_stop_bysize = df_stop.rename(columns={"n": "stop codon"})
    df_stop_bysize["stop codon"] = np.round(
        df_stop_bysize["stop codon"] / n_simulations
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


def tandem_heatmap(filename, family_name, target="real", n_simulations=1):
    """Plots a heatmap of all the possible tandem substitutions"""
    df = pd.read_csv(filename, sep="\t")
    df = df[df["size"] == 2]

    stop_codons = "_stop" in filename

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
        table = np.round(table / n_simulations, 1)

    # some variables
    annot_fmt = ".1f" if stop_codons else ".0f"
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


def tandem_heatmap_percentage(
    real_tandems_fn, count_mut_fn, family_name, n_simulations, limit_perc=False
):
    """Plots a heatmap of all the possible tandem substitutions"""
    df_real = pd.read_csv(real_tandems_fn, sep="\t")
    df_real = df_real[df_real["size"] == 2]

    df_sim = pd.read_csv(count_mut_fn, sep="\t")
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
    table_sim = table_sim / n_simulations

    # get percentage table
    table = np.round(table_sim / table_real * 100, 0)
    table = table.replace([np.inf, -np.inf], np.nan)

    if limit_perc:
        table.clip(upper=100, inplace=True)

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


def tandem_mutation_probability_heatmap(
    tandems_fn, mutation_info_fn, family_name, n_simulations, target="real", base2=True
):
    """"Groups tandems by position mutation probability and reference in each nucleotide of the tandem"""
    df = pd.read_csv(tandems_fn, sep="\t")

    df_prob = pd.read_csv(mutation_info_fn, sep="\t")

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
