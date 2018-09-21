# Scripts to simulate tandem generation in IG sequences

This project aims to simulate and identify "false tandems" produced by the presence of multiple contiguous mutations. The main script is `tandem_generation.py`, it simulates mutations and collect information about clusters multiple mutations. Outputs a table with the information for each cluster found.

Also, in utils directory are some aditional scripts:

- `get_references.py`  finds refereces by name in germline and output them to a single fasta file.

- `get_probabilities.py` counts mutations comparing the sequences with their reference, generating a vector with the probabilites (for each file used as input). Also output another file with the number mutations by sequence (`mutations_by_seq.txt`).

- `get_tandems.py` gets tandems from fasta files comparing the sequences against their reference.

- `summarize_tandem_info.py` creates some summary tables to extract relevant information from the whole tandem simulation files (these files can be big with a high number of simulations).

- `plots.py` some plot functions to show data, designed as a module to be imported by another script or jupyter notebook.

## Usage example

- Install dependencies: `pip install -r requirements.txt`
- Get references: `./get_references.py`
- Get mutation probabilities: `./get_probabilities.py`
- Get tandems from real data: `./get_tandems.py`
- Generate tandem information (1000 simulations): `./tandem_generation.py mutation_info/IGHV1-2-02.tsv fasta/references.fasta mutation_info/mutations_by_seq.txt 1000 tandem_info/simulated/IGHV1-2-02.tsv`
- Summarize simulated info: `./summarize_tandem_info.py`
- Plot using the functions from `plots.py`