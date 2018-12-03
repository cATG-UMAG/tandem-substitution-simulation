#!/usr/bin/env python3
# Generates summary tables to speed up viewing/plotting data later
import re
from multiprocessing import Pool
from glob import glob
from os import path, makedirs
import pandas as pd

INPUT_DIR = "tandem_info/simulated"
OUTPUT_DIR = "tandem_info_summarized"


def main():
    with Pool() as p:
        p.map(summarize, sorted(glob(path.join(INPUT_DIR, "*/*.tsv.gz"))))


def summarize(filename):
    print(filename)
    name = re.sub(r'^.*/|\..*$', '', filename)
    method = filename.split('/')[-2]

    makedirs(path.join(OUTPUT_DIR, method, name), exist_ok=True)

    df = pd.read_table(filename)

    for v in [False, True]:
        suffix = '_stop' if v is True else ''
        subset = df[df.stop_codon == v] 

        if len(subset) > 0:
            # count by simulation/size
            df_gr = pd.DataFrame(subset.groupby(['simulation', 'size']).size()).rename(columns={0: 'n'})
            df_gr.to_csv(path.join(OUTPUT_DIR, method, name, "count_sim_size{}.tsv".format(suffix)), sep='\t')

            # count by simulation/mutation
            df_gr = pd.DataFrame(subset.groupby(['simulation', 'ref', 'alt', 'size']).size()).rename(columns={0: 'n'})
            df_gr.to_csv(path.join(OUTPUT_DIR, method, name, "count_sim_mut{}.tsv".format(suffix)), sep='\t')

            # count by simulation/sequence
            df_gr = pd.DataFrame(subset.groupby(['simulation', 'seq_id']).size()).rename(columns={0: 'n'})
            df_gr.to_csv(path.join(OUTPUT_DIR, method, name, "count_sim_seq{}.tsv".format(suffix)), sep='\t')

            # count by simulation/position
            df_gr = pd.DataFrame(subset.groupby(['simulation', 'pos', 'ref', 'size']).size()).rename(columns={0: 'n'})
            df_gr.to_csv(path.join(OUTPUT_DIR, method, name, "count_sim_pos{}.tsv".format(suffix)), sep='\t')


if __name__ == '__main__':
    main()
