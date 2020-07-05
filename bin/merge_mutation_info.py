#!/usr/bin/env python3
import argparse
import shutil

import pandas as pd


def main():
    args = parse_arguments()

    if len(args.mutation_info) == 1:
        shutil.copyfile(args.mutation_info[0], args.output)
    else:
        nucl_cols = list("ACTG")

        df = pd.concat(
            pd.read_csv(x, sep="\t").assign(n=i)
            for i, x in enumerate(args.mutation_info)
        )

        # get number of total sequences
        total_seqs = (
            df.assign(total_seqs=lambda x: x.mutation_count / x.mutation_probability)
            .groupby("n")
            .agg({"total_seqs": pd.Series.mode})
            .astype(int)
            .sum()
            .values[0]
        )

        # get absolute values
        df[nucl_cols] = df[nucl_cols].multiply(df.mutation_count, axis=0)

        # sum whole dataframe
        df = df.groupby("position").sum().reset_index()

        # re-calculate ratios
        # consider the same number of sequences for all the positions
        df["mutation_probability"] = df.mutation_count / total_seqs
        df[nucl_cols] = df[nucl_cols].divide(df.mutation_count, axis=0)
        df.fillna(0, inplace=True)
        df.drop(columns=["n"], inplace=True)

        # save df
        df.to_csv(args.output, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merges mutation info files")

    parser.add_argument(
        "mutation_info",
        metavar="file",
        action="store",
        nargs="+",
        help="mutation info files (.tsv)",
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="filename",
        action="store",
        default="mutation_info.tsv",
        help="output filename",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
