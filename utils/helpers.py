import re
from itertools import chain
from os import makedirs
from os.path import join, isfile

import pandas as pd


def get_name(filename):
    return re.sub(r"^.*/|[.].*$", "", filename)


def get_all_names(targets):
    names = {}
    for t in targets:
        for f in targets[t]["fasta"]:
            name = get_name(f)
            names[name] = f
    return names


def merge_groups_mutation_info(groups, input_dir, output_dir):
    nucl_cols = ["A", "C", "G", "T"]

    # load all the data in memory if some groups have intersections; it shouldn't be so much
    with open(join(input_dir, "mutations_by_seq.txt")) as f:
        mutations_by_seq = {
            x.split("\t")[0]: x.split("\t")[1].split(",") for x in f.read().splitlines()
        }
    mutation_info = {
        name: pd.read_csv(join(input_dir, f"{name}.tsv"), sep="\t")
        for name in chain.from_iterable(groups.values())
    }

    mutations_by_seq_g = {}
    for group, names in groups.items():
        # concat all dataframes in group
        df = pd.concat(mutation_info[x] for x in names)
        # get means for all relevant columns
        df = df.groupby("position").mean().reset_index().drop("mutation_count", axis=1)
        # normalize A/C/G/T columns
        df[nucl_cols] = (
            df[nucl_cols].divide(df[nucl_cols].sum(axis=1), axis=0).fillna(0)
        )
        # save df
        df.to_csv(join(output_dir, f"{group}.tsv"), sep="\t", index=False)

        # merge mutations_by_seq
        mutations_by_seq_g[group] = chain.from_iterable(
            mutations_by_seq[x] for x in names
        )

    with open(join(output_dir, "mutations_by_seq.txt"), "w") as f:
        f.write("\n".join(f"{k}\t{','.join(v)}" for k, v in mutations_by_seq_g.items()))


def merge_targets_mutation_info(targets, selected, directory):
    nucl_cols = ["A", "C", "G", "T"]

    with open(join(directory, "mutations_by_seq.txt")) as f:
        mutations_by_seq = {
            x.split("\t")[0]: x.split("\t")[1].split(",") for x in f.read().splitlines()
        }

    for target in selected:
        names = [get_name(x) for x in targets[target]["fasta"]]

        # add entry to mutations_by_seq
        mutations_by_seq[target] = list(chain.from_iterable(
            mutations_by_seq[x] for x in names
        ))

        # concat all dataframes in group
        df = pd.concat(
            pd.read_csv(join(directory, f"{x}.tsv"), sep="\t") for x in names
        )
        # get absolute values
        df[nucl_cols] = df[nucl_cols].multiply(df.mutation_count, axis=0)

        # sum whole dataframe
        df = df.groupby("position").sum().reset_index()

        # re-calculate ratios
        # consider the same number of sequences for all the positions
        df["mutation_probability"] = df.mutation_count / len(mutations_by_seq[target])
        df[nucl_cols] = df[nucl_cols].divide(df.mutation_count, axis=0)
        df.fillna(0, inplace=True)

        # save df
        df.to_csv(join(directory, f"{target}.tsv"), sep="\t", index=False)
    with open(join(directory, "mutations_by_seq.txt"), "w") as f:
        f.write("\n".join(f"{k}\t{','.join(v)}" for k, v in mutations_by_seq.items()))


def merge_summaries(input_dirs, output_dir, summary_names, with_stop=False):
    suffix = ["", "_stop"] if with_stop else [""]

    makedirs(output_dir, exist_ok=True)

    for name in summary_names:
        for output_name in [f"{name}{x}.tsv" for x in suffix]:
            subset = [
                join(d, output_name) for d in input_dirs if isfile(join(d, output_name))
            ]
            df = pd.concat(pd.read_csv(p, sep="\t") for p in subset)
            grouping_cols = [x for x in df.columns if x != "n"]

            df = df.groupby(grouping_cols).agg({"n": "sum"}).reset_index()
            df.to_csv(join(output_dir, output_name), sep="\t", index=False)
