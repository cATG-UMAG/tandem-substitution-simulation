#!/usr/bin/env python3
import argparse
import sys
from os import makedirs, path
from typing import List, Optional, Union

import pandas as pd
import yaml


def main():
    args = parse_args()

    if args.output_dir:
        makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.tandem_file, sep="\t")

    with open(args.summary_file) as f:
        summaries = yaml.load(f, Loader=yaml.SafeLoader)

    for x in summaries:
        summarize(df, x["fields"], x["name"], args.output_dir)


def summarize(
    df: pd.DataFrame,
    grouping_vars: List[str],
    output_name: str,
    output_dir: Optional[str] = None,
):
    """
    Makes a summary count table and saves it
    :param df: tandem list dataframe
    :param grouping_vars: list of variables to group by
    :param output_name: output filename base
    :param output_dir: output directory
    """
    if not output_dir:
        output_dir = "."

    for v in [False, True]:
        suffix = "_stop" if v is True else ""
        subset = df[df.stop_codon == v]

        if len(subset) > 0:
            df_gr = (
                subset.groupby(grouping_vars)
                .size()
                .to_frame("n")
                .reset_index()
            )
            df_gr.to_csv(
                path.join(output_dir, f"{output_name}{suffix}.tsv"),
                sep="\t",
                index=False,
            )
            df_gr = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates summary tables (counts) to speed up viewing/plotting data later"
    )

    parser.add_argument(
        "tandem_file",
        metavar="file",
        action="store",
        help="tandem list file (.tsv[.gz])",
    )

    parser.add_argument(
        "summary_file",
        metavar="file",
        action="store",
        help="summaries list file (.yml)",
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        metavar="output_dir",
        action="store",
        help="directory to save output files",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
