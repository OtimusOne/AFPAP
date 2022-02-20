#! /bin/bash

python scripts/AFPAP_sequence_analysis.py -i input/LGB.fasta
python scripts/AFPAP_secondary_structure.py -i input/ranked_0.pdb
python scripts/AFPAP_structure_analysis.py -i ./output/work/ranked_0_secondary_structure.pdb
multiqc -f -c ./config/multiqc_config.yaml --custom-css-file ./config/multiqc_custom_css.css -o output/ .