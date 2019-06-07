configfile: "config.yaml"

from utils.helpers import get_name, get_all_names, merge_mutation_info, merge_summaries
from utils.constants import SUMMARY_DETAILS

NAMES = get_all_names(config["targets"])
SUMMARIES = [x[1] for x in SUMMARY_DETAILS]


rule all:
    input:
        expand("tandem_info/single/simulated/{target}.tsv.gz",
            target=config["targets"].keys()),
        expand("tandem_info_summarized/single/{target}/{summary}.tsv",
            target=config["targets"].keys(), summary=SUMMARIES),
        expand("tandem_info_summarized/grouped/{group}/{summary}.tsv",
            group=config["groups"].keys(), summary=SUMMARIES),
        expand("tandem_info/grouped/real/{group}.tsv", group=config["groups"].keys()),
        expand("mutation_info/grouped/{name}.tsv", name=config["groups"].keys()),
        "mutation_info/grouped/mutations_by_seq.txt"


rule get_real_tandems:
    input:
        fasta_file = lambda w: NAMES[w.name],
        reference_file = config["references"]
    output:
        "tandem_info/single/real/{name}.tsv"
    shell:
        "python3 scripts/get_tandems.py {input.fasta_file} {input.reference_file} {output}"


rule get_mutation_probability:
    input:
        fasta_files = NAMES.values(),
        reference_file = config["references"]
    output:
        expand("mutation_info/single/{name}.tsv", name=NAMES.keys()),
        "mutation_info/single/mutations_by_seq.txt"
    params:
        outdir = "mutation_info/single/"
    shell:
        "python3 scripts/get_probabilities.py {input.reference_file} {params.outdir} {input.fasta_files}"


rule make_simulation:
    input:
        mutation_info = "mutation_info/single/{target}.tsv",
        references = config["references"],
        mutations_by_seq = "mutation_info/single/mutations_by_seq.txt"
    output:
        "tmp/tandem_info_simulated/{target}.tsv"
    threads: 4
    params:
        reference_names = lambda w: [get_name(x) for x in config["targets"][w.target]["fasta"]],
        simulations = config["simulations"]
    shell:
        "python3 scripts/tandem_generation.py {input.mutation_info} {input.references} {input.mutations_by_seq} \
         {params.simulations} {output} -r {params.reference_names} -t {threads}"


rule sort_and_compress:
    input:
        "tmp/tandem_info_simulated/{target}.tsv"
    output:
        "tandem_info/single/simulated/{target}.tsv.gz"
    shell:
        "sort -nk1 -nk2 {input}|gzip > {output}"


rule group_mutation_info:
    input:
        expand("mutation_info/single/{name}.tsv", name=NAMES.keys()),
        "mutation_info/single/mutations_by_seq.txt"
    output:
        expand("mutation_info/grouped/{name}.tsv", name=config["groups"].keys()),
        "mutation_info/grouped/mutations_by_seq.txt"
    params:
        input_dir = "mutation_info/single/",
        output_dir = "mutation_info/grouped/"
    run:
        merge_mutation_info(config["groups"], params.input_dir, params.output_dir)


rule group_real_tandem_info:
    input:
        lambda w: [f"tandem_info/single/real/{x}.tsv" for x in config["groups"][w.group]]
    output:
        "tandem_info/grouped/real/{group}.tsv"
    shell:
        "head -n 1 {input[0]} > {output} && \
        for f in {input}; do tail -q -n +2 $f >> {output}; done"


rule summarize:
    input:
        "tmp/tandem_info_simulated/{target}.tsv"
    output:
        expand("tandem_info_summarized/single/{{target}}/{summary}.tsv", summary=SUMMARIES)
    threads: 2
    params:
        output_dir = "tandem_info_summarized/single/{target}"
    shell:
        "PYTHONPATH=. python3 scripts/summarize_tandem_info.py {input} {params.output_dir} {threads}"


rule group_summaries:
    input:
        lambda w: expand(
            "tandem_info_summarized/single/{target}/{summary}.tsv",
            target=config["groups"][w.group], summary=SUMMARIES)
    output:
        expand("tandem_info_summarized/grouped/{{group}}/{summary}.tsv", summary=SUMMARIES)
    params:
        input_dirs = lambda w: [f"tandem_info_summarized/single/{x}" for x in config["groups"][w.group]],
        output_dir = "tandem_info_summarized/grouped/{group}"
    run:
        merge_summaries(params.input_dirs, params.output_dir)
