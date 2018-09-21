#!/usr/bin/env python3
# Generates summary tables to speed up viewing/plotting data later
import re
from multiprocessing import Pool
from glob import glob
from os import mkdir, path
import pandas as pd

INPUT_DIR = "tandem_info/simulated"
OUTPUT_DIR = "tandem_info_summarized"


def main():
    if not path.isdir(OUTPUT_DIR):
        mkdir(OUTPUT_DIR)

    with Pool(4) as p:
        p.map(summarize, sorted(glob(path.join(INPUT_DIR, "*2/*.tsv.gz"))))


def summarize(filename):
    print(filename)
    name = re.sub(r'^.*/|\..*$', '', filename)
    method = filename.split('/')[-2]

    if not path.isdir(path.join(OUTPUT_DIR, name)):
        mkdir(path.join(OUTPUT_DIR, name))

    df = pd.read_table(filename)

    # count by simulation/size
    df_gr = pd.DataFrame(df.groupby(['simulation', 'size']).size()).rename(columns={0: 'n'})
    df_gr.to_csv(path.join(OUTPUT_DIR, name, "count_sim_size_{}.tsv".format(method)), sep='\t')

    # count by simulation/mutation
    df_gr = pd.DataFrame(df.groupby(['simulation', 'ref', 'alt', 'size']).size()).rename(columns={0: 'n'})
    df_gr.to_csv(path.join(OUTPUT_DIR, name, "count_sim_mut_{}.tsv".format(method)), sep='\t')

    # count by simulation/sequence
    df_gr = pd.DataFrame(df.groupby(['simulation', 'seq_id']).size()).rename(columns={0: 'n'})
    df_gr.to_csv(path.join(OUTPUT_DIR, name, "count_sim_seq_{}.tsv".format(method)), sep='\t')


if __name__ == '__main__':
    main()
