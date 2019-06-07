import re
from itertools import chain
from os import listdir, makedirs
from os.path import join

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


def merge_mutation_info(groups, input_dir, output_dir):
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
        df.groupby("position").mean().reset_index().drop("mutation_count", axis=1)
        # normalize A/C/G/T columns
        df[["A", "C", "G", "T"]] = (
            df[["A", "C", "G", "T"]]
            .divide(df[["A", "C", "G", "T"]].sum(axis=1), axis=0)
            .fillna(0)
        )
        # save df
        df.to_csv(join(output_dir, f"{group}.tsv"), sep="\t", index=False)

        # merge mutations_by_seq
        mutations_by_seq_g[group] = chain.from_iterable(
            mutations_by_seq[x] for x in names
        )

    with open(join(output_dir, "mutations_by_seq.txt"), "w") as f:
        f.write("\n".join(f"{k}\t{' '.join(v)}" for k, v in mutations_by_seq_g.items()))


def merge_summaries(input_dirs, output_dir):
    makedirs(output_dir, exist_ok=True)

    for output_name in listdir(input_dirs[0]):
        df = pd.concat(pd.read_csv(join(d, output_name), sep="\t") for d in input_dirs)
        grouping_cols = [x for x in df.columns if x != "n"]

        df = df.groupby(grouping_cols).agg({"n": "sum"}).reset_index()
        df.to_csv(join(output_dir, output_name), sep="\t", index=False)
