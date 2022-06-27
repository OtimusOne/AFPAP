<p align="center">
  <a href="" rel="noopener">
 <img width=400px src="./docs/logo2.png" alt="Project logo"></a>
</p>

<h3 align="center">AlphaFold-based Protein Analysis Pipeline</h3>

<div align="center">

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

## 📝 Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## Installation <a name = "installation"></a>

### Prerequisites

<ul>
<li>Python  ≥ 3.8</li>
<li>Java  ≥ 11</li>
<li>git, pip, conda/mamba</li>
</ul>

### Mandatory installation:
```
conda config --env --add channels conda-forge
conda config --env --add channels anaconda
conda config --env --add channels bioconda

conda install 'numpy>=1.18.5' 'pandas>=1.4.1' 'biopython>=1.76' 'multiqc>=1.12' pymol-open-source=2.5.0
conda install -c salilab dssp=3.0.0
```
- Install Nextflow and add the executable to PATH: https://github.com/nextflow-io/nextflow
```
curl -fsSL https://get.nextflow.io | bash
```
- Clone this repository:
```
git clone https://github.com/OtimusOne/AFPAP.git
```
### Optional component - AlphaFold:
- AlphaFold is used to predicted the protein structure if a PDB file is not provided. We use the LocalColabFold implementation in order to avoid the large databases used by native AlphaFold.
- Follow the install intructions:
    - https://github.com/YoshitakaMo/localcolabfold
- If ColabFold is not installed the user must provide a PDB structure using the *--pdb* argument if structural analysis is desired. Set *skipAlphaFold = true* inside **nextflow.config**.

### Optional component - Pfam:
- The Pfam database is used to match the protein sequence agains protein families.
- Download the Pfam database(~1.5GB) 
http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
- Install the following packages:
```
conda install pfam_scan perl-json
```
- Set path variable inside **nextflow.config**
```
params {
    ...
    pfam_path="/path/to/Pfam/directory"
    ...
}
```
- If Pfam is not used set *skipPfamSearch = true* inside **nextflow.config**.

### Optional component - Conservation MSA:
- Install Blast, Muscle and CD-HIT
```
conda install blast muscle cd-hit
```
- Create a BLAST database from which to generate the MSA file. You can use any FASTA database such as SwissProt/TrEMBL/UniRef50/UniRef90 (*warning: SwissProt might not return enough hits to calculate the sequence conservation*). Example setting UniRef90 up as a Blast database: https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/
```
curl {https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz} | gunzip | makeblastdb -out uniref90 -dbtype prot -title UniRef90 -parse_seqids
```
- Set path variable inside **nextflow.config**
```
params {
    ...
    blast_path="/path/to/BLASTDB/uniref90"
    ...
}
```
- If this component is not used set *skipConservationMSA = true* inside **nextflow.config**.

### Optional component - SimBa2:
- SimBa2 is used to predict the effect of point mutation on protein stability.
```
python -m pip install https://github.com/kasperplaneta/SimBa2/archive/main.tar.gz
```
- If SimBa2 is not installed set *skipPointMutations = true* inside **nextflow.config**.

### Optional component - P2Rank:
- P2Rank is used to predict ligand-binding pockets from the protein structure.
- Follow the install intructions:
https://github.com/rdk/p2rank
- Pillow is required for generating the pocket visualizations:
```
conda install pillow
```
- If P2Rank is not installed set *skipPocketPrediction = true* inside **nextflow.config**.

### Optional component - AutoDock Vina:
- AutoDock Vina is used to dock the ligands provided with the *--ligands* argument.
- Install ADFR suite - MSMS not required
https://ccsb.scripps.edu/adfr/downloads/
- Install Vina python bindings:
```
pip install vina
```
- If P2Rank is installed Vina will dock the ligands to each of the predicted pockets, otherwise it will execute only blind docking.
- If Vina is not installed set *skipMolecularDocking = true* inside **nextflow.config**.

## Usage <a name="usage"></a>

Usage example:
```
nextflow run main.nf --fasta input.fasta
nextflow run main.nf --pdb input.pdb --ligands ligand.mol2
```
For a full list of parameters run:
```
nextflow run main.nf --help
```


## Authors <a name = "authors"></a>

- [@Maghiar Octavian](https://github.com/OtimusOne) 

## Acknowledgements <a name = "acknowledgement"></a>
If you find this work useful please properly cite each of the used tools.
