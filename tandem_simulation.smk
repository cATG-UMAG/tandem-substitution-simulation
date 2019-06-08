configfile: "config.yaml"

from os.path import join
from utils.helpers import get_name, get_all_names, merge_mutation_info, merge_summaries
from utils.constants import SUMMARY_DETAILS

BASEDIR = config["base_dir"]
NAMES = get_all_names(config["targets"])
SUMMARIES = [x[1] for x in SUMMARY_DETAILS]
GROUPS = config["groups"] if "groups" in config and config["groups"] else {}


rule all:
    input:
        expand(join(BASEDIR, "tandem_info/single/simulated/{target}.tsv.gz"),
            target=config["targets"].keys()),
        expand(join(BASEDIR, "tandem_info_summarized/single/{target}/{summary}.tsv"),
            target=config["targets"].keys(), summary=SUMMARIES),
        expand(join(BASEDIR, "tandem_info_summarized/grouped/{group}/{summary}.tsv"),
            group=GROUPS.keys(), summary=SUMMARIES),
        expand(join(BASEDIR, "tandem_info/grouped/real/{group}.tsv"), group=GROUPS.keys()),
        expand(join(BASEDIR, "mutation_info/grouped/{name}.tsv"), name=GROUPS.keys()),
        join(BASEDIR, "mutation_info/grouped/mutations_by_seq.txt") if len(GROUPS) > 0 else []


rule get_real_tandems:
    input:
        fasta_file = lambda w: join(BASEDIR, NAMES[w.name]),
        reference_file = join(BASEDIR, config["references"])
    output:
        join(BASEDIR, "tandem_info/single/real/{name}.tsv")
    shell:
        "python3 scripts/get_tandems.py {input.fasta_file} {input.reference_file} {output}"


rule get_mutation_probability:
    input:
        fasta_files = [join(BASEDIR, x) for x in NAMES.values()],
        reference_file = join(BASEDIR, config["references"])
    output:
        expand(join(BASEDIR, "mutation_info/single/{name}.tsv"), name=NAMES.keys()),
        join(BASEDIR, "mutation_info/single/mutations_by_seq.txt")
    params:
        outdir = join(BASEDIR, "mutation_info/single/")
    shell:
        "python3 scripts/get_probabilities.py {input.reference_file} {params.outdir} {input.fasta_files}"


rule make_simulation:
    input:
        mutation_info = join(BASEDIR, "mutation_info/single/{target}.tsv"),
        references = join(BASEDIR, config["references"]),
        mutations_by_seq = join(BASEDIR, "mutation_info/single/mutations_by_seq.txt")
    output:
        temp(join(BASEDIR, "tmp/tandem_info_simulated/{target}.tsv"))
    threads: 4
    params:
        reference_names = lambda w: [get_name(x) for x in config["targets"][w.target]["fasta"]],
        simulations = config["simulations"]
    shell:
        "python3 scripts/tandem_generation.py {input.mutation_info} {input.references} {input.mutations_by_seq} \
         {params.simulations} {output} -r {params.reference_names} -t {threads}"


rule sort_and_compress:
    input:
        join(BASEDIR, "tmp/tandem_info_simulated/{target}.tsv")
    output:
        join(BASEDIR, "tandem_info/single/simulated/{target}.tsv.gz")
    shell:
        "sort -nk1 -nk2 {input}|gzip > {output}"


rule group_mutation_info:
    input:
        expand(join(BASEDIR, "mutation_info/single/{name}.tsv"), name=NAMES.keys()),
        join(BASEDIR, "mutation_info/single/mutations_by_seq.txt")
    output:
        expand(join(BASEDIR, "mutation_info/grouped/{name}.tsv"), name=GROUPS.keys()),
        join(BASEDIR, "mutation_info/grouped/mutations_by_seq.txt") if len(GROUPS) > 0 else []
    params:
        input_dir = join(BASEDIR, "mutation_info/single/"),
        output_dir = join(BASEDIR, "mutation_info/grouped/")
    run:
        if len(GROUPS) > 0:
            merge_mutation_info(GROUPS, params.input_dir, params.output_dir)


rule group_real_tandem_info:
    input:
        lambda w: [join(BASEDIR, f"tandem_info/single/real/{x}.tsv") for x in GROUPS[w.group]]
    output:
        BASEDIR + "tandem_info/grouped/real/{group}.tsv"
    shell:
        "head -n 1 {input[0]} > {output} && \
        for f in {input}; do tail -q -n +2 $f >> {output}; done"


rule summarize:
    input:
        join(BASEDIR, "tmp/tandem_info_simulated/{target}.tsv")
    output:
        expand(join(BASEDIR, "tandem_info_summarized/single/{{target}}/{summary}.tsv"), summary=SUMMARIES)
    threads: 2
    params:
        output_dir = join(BASEDIR, "tandem_info_summarized/single/{target}")
    shell:
        "PYTHONPATH=. python3 scripts/summarize_tandem_info.py {input} {params.output_dir} {threads}"


rule group_summaries:
    input:
        lambda w: expand(
            join(BASEDIR, "tandem_info_summarized/single/{target}/{summary}.tsv"),
            target=GROUPS[w.group], summary=SUMMARIES)
    output:
        expand(join(BASEDIR, "tandem_info_summarized/grouped/{{group}}/{summary}.tsv"), summary=SUMMARIES)
    params:
        input_dirs = lambda w: [join(BASEDIR, f"tandem_info_summarized/single/{x}") for x in GROUPS[w.group]],
        output_dir = join(BASEDIR, "tandem_info_summarized/grouped/{group}")
    run:
        merge_summaries(params.input_dirs, params.output_dir)
