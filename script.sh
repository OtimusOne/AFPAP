#! /bin/bash


while getopts i:p: flag
do
    case "${flag}" in
        i) inputFile=${OPTARG};;
        p) pdbFile=${OPTARG};;
    esac
done
if [ ! -f "$inputFile" ] || [ ! -f "$pdbFile" ] ; then
    echo "No arguments supplied"
    exit 1
fi


echo "Input: $inputFile";
echo "PDB: $pdbFile";

baseDir="$PWD"
mkdir -p output/work

python scripts/AFPAP_sequence_analysis.py -i $inputFile
python scripts/AFPAP_secondary_structure.py -i $pdbFile
python scripts/AFPAP_structure_analysis.py -i ./output/work/proteinStructure_ss.pdb

prank predict -f output/work/proteinStructure_ss.pdb -o output/work/ -c alphafold

python scripts/AFPAP_p2rank_visualization.py -p output/work/visualizations/proteinStructure_ss.pdb.pml -c output/work/proteinStructure_ss.pdb_predictions.csv

cd output/work/visualizations
pymol -cq proteinStructure_ss.pdb.pml

cd $baseDir

multiqc -f -c ./config/multiqc_config.yaml --custom-css-file ./config/multiqc_custom_css.css -o output/ .