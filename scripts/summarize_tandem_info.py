#!/usr/bin/env python3
# Generates summary tables (counts) to speed up viewing/plotting data later
import sys
from os import makedirs, path
from utils.constants import SUMMARY_DETAILS

import pandas as pd


def main():
    if len(sys.argv) < 3:
        print("./summarize_tandem_info.py <tandem_file.tsv> <output_dir>")
        exit(-1)

    filename, output_dir = sys.argv[1:3]

    makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(filename, sep="\t")

    for x in SUMMARY_DETAILS:
        summarize(df, x[0], output_dir, x[1])


def summarize(df, grouping_vars, output_dir, output_name):
    for v in [False, True]:
        suffix = "_stop" if v is True else ""
        subset = df[df.stop_codon == v]

        if len(subset) > 0:
            df_gr = (
                subset.groupby(grouping_vars)
                .size()
                .reset_index()
                .rename(columns={0: "n"})
            )
            df_gr.to_csv(
                path.join(output_dir, f"{output_name}{suffix}.tsv"),
                sep="\t",
                index=False,
            )


if __name__ == "__main__":
    main()
