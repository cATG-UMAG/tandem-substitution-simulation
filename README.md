# Scripts to simulate tandem generation in IGH sequences

- `get_references.py`  finds refereces by name in germline and output them to a single fasta file.

- `get_probabilities.py` counts mutations comparing the sequences with their reference, generating a vector with the probabilites (for each file used as input). Also output another file with the number mutations by sequence (`mutations_by_seq.txt`).

- `tandem_generation.py` simulates mutations and collect information about clusters of mutations with size 2 or higher. Outputs a table with the information for each cluster found


## Usage example
- Install dependencies: `pip install -r requirements.txt`
- Get references: `./get_references.py`
- Get probabilities: `./get_probabilities.py`
- Generate tandem information (1000 simulations): `./tandem_generation.py mutation_info/IGHV1-2-02.tsv fasta/references.fasta mutation_info/mutations_by_seq.txt 1000 mutation_info/IGHV1-2-02.tsv`

- To view the notebooks `jupyter` must be installed. If `jupyterlab` is used, [jupyterlabs' plotly extension](https://github.com/jupyterlab/jupyter-renderers/tree/master/packages/plotly-extension) is needed too.