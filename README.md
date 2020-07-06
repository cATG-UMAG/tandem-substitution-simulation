# Tandem Substitutions Simulation

The purpose of the project is to provide a easy way to simulate the presence of "false" tandem substitutions, that is, when multiple adjacent nucleotides mutate on independent events. This can be used improve accuracy in the identification of real tandems. It was designed to work with immune repertoire data, but it should work with other data too.

It consists on a set of python scripts organized in a nextflow pipeline. This allows to run the whole process easily with any set of data and also executing it in a distributed computing environment.

As a brief summary, this pipeline does the following:

1. Compare a list of sequences against its reference by position, collecting all the substitutions found.
2. The list of substitutions are summarized in a mutational frequency value and a nucleotide frequency for each position.
3. The mutational frequencies are then used as probability distribution to simulate the exact input dataset (with the same number of sequences and the same same number of substitutions by sequences) generating each substitution as an independent event. This is done a large number of times to provide robustness in the results.

## Quick Start

1. [Prepare you data](#input-data)
2. Install [nextflow](https://www.nextflow.io)
3. Clone this repository (it's possible to run the pipeline without cloning it too, more info in [this link](https://www.nextflow.io/docs/latest/sharing.html#running-a-pipeline))
4. Modify `params.yml` [to your needs](#pipeline-parameters)
5. Choose a profile according to your execution environment
6. Run the pipeline (for example with singularity profile):
   ```
   nextflow run -params-file params.yml -profile singularity main.nf
   ```
7. Check `output/` directory after the pipeline execution ends

## Input Data

The data needed to run the pipeline is basically a one or more FASTA files each one with a set of sequences aligned against a reference sequence. It's important for the sequences to be previously aligned and even have the same length because the pipeline doesn't do alignment itself.

Another FASTA file must have all the references. The filenames of each of the sequences must be the same one that has the reference sequence, for example, if a sequence in the references FASTA has for name `IGHV3-7_01` the file that has the sequences aligned against it must have for name `IGHV3-7_01.fasta`.

The input sequences also must start in the correct reading frame to being able to traduce them properly (and detect amino acid changes).

## Pipeline Parameters

In the `params.yml` file you will need to specify a set of configurations to adapt the pipeline to suit your use case. These are:

- `references`: the path of the FASTA file containing the references.
- `n_simulations`: the number of times that the whole input dataset is simulated.
- `with_stopcodons`: when `False`, only productive rearrangements will be produced in the simulation.
- `with_summaries`: given that the simulation output is a big list of tandems, you might want to summarize this data. When this param is set to `True`, this will be done as a step of the pipeline, that does basically `group by`+`count` by a set of fields. Check next param to know how to specify them.
- `summary_list`: the path of the summaries specification files. This must be a YAML file, having the following structure:
  ```yaml
  - name: count_sim_size
    fields:
      - simulation
      - size
  - name: count_size
    fields:
      - size
  ```
  Here, each entry must have a `name` property that will be used for naming output files, and a `fields` property having a list of fields to make the groups (these obviously must be a subset of the columns of the tandem output list).
- `targets`: in this param you must specify the simulation targets and their corresponding files. It should look like this:

  ```yaml
  targets:
    IGHV1-2_02:
      - input/fasta/IGHV1-2_02.fasta
    IGHV1-69_01:
      - input/fasta/IGHV1-69_01.fasta
    IGHV4-30:
      - input/fasta/IGHV4-30-2_01.fasta
      - input/fasta/IGHV4-30-4_01.fasta
    IGKV3-15_01:
      - input/fasta/IGKV3-15_01.fasta
  ```

  As you can see, it's a mapping with a name as key and a list of files as values. The "key" will be used to name output files, and the files associated with it are obviously the files that will be used to simulate this target.

  It's possible to specify multiple files for a target, this could help you don't have enough sequences and want to improve the robustness of the simulation. Each one of this files will use a different reference to get the mutation metrics, and then these mutation metrics will be joined in a single one that will be used as input for the simulation. But it's important that all of the input sequences groups in the target must be aligned against each other and also have the same length for the mutational metrics to be joined consistently.

- `groups`: if you want to summarize the output information by a higher grouping level, you can use this parameter to specify the groups using this structure:
  ```yaml
  groups:
    IGH:
      - IGHV1-2_02
      - IGHV1-69_01
      - IGHV4-30
    IGV:
      - IGHV1-2_02
      - IGHV1-69_01
      - IGHV4-30
      - IGKV3-15_01
  ```
  Each group must have by value a list of targets. For each group, observed tandems, mutational metrics and summaries of simulated tandems will be joined. If you don't want groups, then delete the `groups` entry from params file or leave it empty.

It is recommended to use absolute paths for all the params involving files, but it works with relative paths if you specify them correctly.

## Pipeline Output

If the pipeline ended successfully, you will find the following:

- Mutational metrics used for the simulation (`mutation_info/` directory). Here you will find a `.tsv` file for each target and a `mutations_by_seq.txt` file containing the number of substitutions of each sequence on the dataset.
- Observed tandems in the input data (`tandem_info/observed` directory). For each target, you will find a list (`.tsv`) of tandems with:
  - sequence id
  - position
  - reference allele
  - alternative allele
  - tandem size
  - context (two reference nucleotides at each side of the substitution)
- Simulated tandems (`tandem_info/simulated` directory). For each target, you will find a compressed `.tsv` file containing the list of the tandem substitutions produced by the simulations. The fields of this table are:
  - simulation id
  - number of substitution in the original sequence that correspond this tandem
  - sequence id (a correlative number for each of the sequences in the dataset)
  - position
  - reference allele
  - alternative allele
  - tandem size
  - context (two reference nucleotides at each side of the substitution)
  - frame position of the tandem (1 is at the start of the codon, 3 at the end)
  - amino acid change
  - if the substitution produces a stop codon
- A list of summary tables for each target if you enabled summaries (`tandem_info_summarized` directory).

For each of these items you will find a "single" directory and also you could find a "grouped" one if you provided groups in the configuration.

## Reproducibility

In the `env/` directory you will find several ways to reproduce a working environment for executing this pipeline. These include:

- A specification for a Conda environment (`environment.yml`).
- Specifications for building a Docker or a Singularity container. These also have automated builds hosted in DockerHub/SingularityHub respectively, you can find the corresponding URLs in the provided `nextflow.config` file.
- And finally a `requirements.txt` file if you want to use simply a Python virtualenv.

The config provided has defined some profiles to use Conda, Docker or Singularity directly from nextflow, you only need to specify it on execution (for example `-profile conda`).
