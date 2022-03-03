# AFPAP
AlphaFold-based Protein Analysis Pipeline

## Installation
1. AlphaFold
    https://github.com/deepmind/alphafold
    https://github.com/kalininalab/alphafold_non_docker
2. P2RANK
    https://github.com/rdk/p2rank
3. ADFR suite - MSMS not required
    https://ccsb.scripps.edu/adfr/downloads/
4. AFPAP
```
conda install -c conda-forge -c anaconda -c bioconda numpy pandas biopython multiqc pymol-open-source pillow
pip install vina
sudo apt-get install dssp
```