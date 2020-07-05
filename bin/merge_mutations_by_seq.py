#!/usr/bin/env python3
import argparse
from itertools import chain
import shutil


def main():
    args = parse_arguments()

    if len(args.mutations_by_seq) == 1:
        shutil.copyfile(args.mutations_by_seq[0], args.output)
    else:
        # load all files
        mutations_by_seq = {}
        for filename in args.mutations_by_seq:
            with open(filename) as f:
                for line in f.read().splitlines():
                    line_splitted = line.split("\t")
                    mutations_by_seq[line_splitted[0]] = line_splitted[1].split(",")

        name = args.name if args.name else "|".join(mutations_by_seq.keys())

        with open(args.output, "w") as f:
            f.write(
                f"{name}\t{','.join(chain.from_iterable(mutations_by_seq.values()))}"
            )


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merges mutations by sequence files")

    parser.add_argument(
        "mutations_by_seq",
        metavar="file",
        action="store",
        nargs="+",
        help="mutations by sequence files (.txt)",
    )

    parser.add_argument(
        "--name", "-n", action="store", help="target name in the output file",
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="filename",
        action="store",
        default="mutations_by_seq.txt",
        help="output filename",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
