#!/usr/bin/env python3
import sys
from os import listdir, makedirs
from os.path import join

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import rc

sys.path.append("../")
from utils.plots import (
    cluster_size_distribution,
    mutation_probability,
    mutations_byseq_distplot,
    tandem_distplot,
    tandem_heatmap,
    tandem_heatmap_percentage,
    stop_codon_distribution,
)


def main():
    if len(sys.argv) < 2:
        print("./make_plots.py <base_dir>")
        exit(-1)

    base_dir = sys.argv[1]

    # settings
    sns.set_context("notebook", font_scale=1.3, rc={"lines.linewidth": 2})
    rc("font", **{"family": "sans-serif", "sans-serif": ["Roboto"]})
    formats = ["pdf"]
    n_simulations = 100_000

    makedirs(join(base_dir, "output/figures/"), exist_ok=True)

    for kind in listdir(join(base_dir, "tandem_info_summarized/")):
        for name in listdir(join(base_dir, "tandem_info_summarized/", kind)):
            name_display = name.replace("_", "*")
            print(f"{name_display}...")
            makedirs(join(base_dir, "output/figures", name), exist_ok=True)

            savefig_multiformat(
                tandem_distplot(
                    join(base_dir, f"tandem_info/{kind}/real/{name}.tsv"),
                    join(
                        base_dir,
                        f"tandem_info_summarized/{kind}/{name}/count_sim_size.tsv",
                    ),
                    name_display,
                ),
                join(base_dir, f"output/figures/{name}/tandem_dist"),
                formats,
            )

            savefig_multiformat(
                mutations_byseq_distplot(
                    join(base_dir, f"mutation_info/{kind}/mutations_by_seq.txt"),
                    name_display,
                ),
                join(base_dir, f"output/figures/{name}/mutations_byseq"),
                formats,
            )

            savefig_multiformat(
                mutation_probability(
                    join(base_dir, f"mutation_info/{kind}/{name}.tsv"), name_display
                ),
                join(base_dir, f"output/figures/{name}/mutation_probability"),
                formats,
            )

            savefig_multiformat(
                cluster_size_distribution(
                    join(base_dir, f"tandem_info/{kind}/real/{name}.tsv"),
                    join(
                        base_dir, f"tandem_info_summarized/{kind}/{name}/count_size.tsv"
                    ),
                    name_display,
                    n_simulations,
                ),
                join(base_dir, f"output/figures/{name}/cluster_size"),
                formats,
            )

            savefig_multiformat(
                stop_codon_distribution(
                    join(
                        base_dir, f"tandem_info_summarized/{kind}/{name}/count_size.tsv"
                    ),
                    join(
                        base_dir,
                        f"tandem_info_summarized/{kind}/{name}/count_size_stop.tsv",
                    ),
                    name_display,
                    n_simulations,
                ),
                join(base_dir, f"output/figures/{name}/stop_codons"),
                formats,
            )

            savefig_multiformat(
                tandem_heatmap(
                    join(base_dir, f"tandem_info/{kind}/real/{name}.tsv"),
                    name_display,
                    "real",
                ),
                join(base_dir, f"output/figures/{name}/tandem_heatmap_real"),
                formats,
            )

            savefig_multiformat(
                tandem_heatmap(
                    join(
                        base_dir,
                        f"tandem_info_summarized/{kind}/{name}/count_mutation.tsv",
                    ),
                    name_display,
                    "simulated",
                    n_simulations,
                ),
                join(base_dir, f"output/figures/{name}/tandem_heatmap_simulated"),
                formats,
            )

            savefig_multiformat(
                tandem_heatmap(
                    join(
                        base_dir,
                        f"tandem_info_summarized/{kind}/{name}/count_mutation_stop.tsv",
                    ),
                    name_display,
                    "simulated",
                    n_simulations,
                ),
                join(base_dir, f"output/figures/{name}/tandem_heatmap_simulated_stop"),
                formats,
            )

            savefig_multiformat(
                tandem_heatmap_percentage(
                    join(base_dir, f"tandem_info/{kind}/real/{name}.tsv"),
                    join(
                        base_dir,
                        f"tandem_info_summarized/{kind}/{name}/count_mutation.tsv",
                    ),
                    name_display,
                    n_simulations,
                    True,
                ),
                join(base_dir, f"output/figures/{name}/tandem_heatmap_percentage"),
                formats,
            )

            plt.close("all")


def savefig_multiformat(fig, path, formats):
    for f in formats:
        fig.savefig(f"{path}.{f}", bbox_inches="tight")


if __name__ == "__main__":
    main()
