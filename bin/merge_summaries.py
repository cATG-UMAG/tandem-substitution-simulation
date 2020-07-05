#!/usr/bin/env python3

import argparse
from os import makedirs
from os.path import isfile, join
from posix import listdir
import shutil

import pandas as pd
import yaml


def main():
    args = parse_arguments()

    # if any file has _stop suffix, then summaries includes stop codons
    if any("_stop" in x for x in listdir(args.summary_directories[0])):
        suffix = ["", "_stop"]
    else:
        suffix = [""]

    if args.output_dir:
        makedirs(args.output_dir, exist_ok=True)
    else:
        args.output_dir = "."

    with open(args.summary_file) as f:
        summary_names = [x["name"] for x in yaml.load(f, Loader=yaml.SafeLoader)]

    for name in summary_names:
        for output_name in [f"{name}{x}.tsv" for x in suffix]:
            if len(args.summary_directories) > 1:
                # concat all dataframes maching output name in input dirs
                df = pd.concat(
                    [
                        pd.read_csv(join(d, output_name), sep="\t")
                        for d in args.summary_directories
                        if isfile(join(d, output_name))
                    ],
                    ignore_index=True,
                )
                grouping_cols = [x for x in df.columns if x != "n"]

                # group by all variables and sum n
                df = df.groupby(grouping_cols).agg({"n": "sum"}).reset_index()
                df.to_csv(join(args.output_dir, output_name), sep="\t", index=False)
            else:
                # if only one summary directory provided, just copy files
                if isfile(join(args.summary_directories[0], output_name)):
                    shutil.copyfile(
                        join(args.summary_directories[0], output_name),
                        join(args.output_dir, output_name),
                    )


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merges summaries")

    parser.add_argument(
        "summary_file",
        metavar="file",
        action="store",
        help="summaries list file (.yml)",
    )

    parser.add_argument(
        "summary_directories",
        metavar="directory",
        action="store",
        nargs="+",
        help="diretories containing summaries",
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
