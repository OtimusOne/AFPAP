#! /bin/bash

while getopts i:p:l: flag
do
    case "${flag}" in
        i) inputFile=${OPTARG};;
        p) pdbFile=${OPTARG};;
        l) ligandFile=${OPTARG};;
    esac
done
if [ ! -f "$inputFile" ] || [ ! -f "$pdbFile" ] ; then
    echo "No arguments supplied"
    exit 1
fi

echo "Input: $inputFile";
echo "PDB: $pdbFile";

baseDir="$PWD"
mkdir -p output/work/docking

python scripts/AFPAP_sequence_analysis.py -i $inputFile
python scripts/AFPAP_secondary_structure.py -i $pdbFile
python scripts/AFPAP_structure_analysis.py -i ./output/work/proteinStructure_ss.pdb

prank predict -f output/work/proteinStructure_ss.pdb -o output/work/ -c alphafold

python scripts/AFPAP_p2rank_visualization.py -p output/work/visualizations/proteinStructure_ss.pdb.pml -c output/work/proteinStructure_ss.pdb_predictions.csv

cd output/work/visualizations
pymol -cq proteinStructure_ss.pdb.pml
cd $baseDir
python scripts/AFPAP_p2rank_gallery.py

if [ -f "$ligandFile" ] ; then
    cp ./output/work/proteinStructure_ss.pdb ./output/work/docking/receptor.pdb
    ligand_base_name=$(basename ${ligandFile})    
    cp $ligandFile ./output/work/docking/
    cd ./output/work/docking
    prepare_receptor -r receptor.pdb -o receptor.pdbqt -A "hydrogens"
    prepare_ligand -l ./output/work/docking/$ligand_base_name -o ligand.pdbqt
    cd $baseDir
    python scripts/AFPAP_molecular_docking.py -r ./output/work/docking/receptor.pdbqt -l ./output/work/docking/ligand.pdbqt -n $ligand_base_name
fi

multiqc -f -c ./config/multiqc_config.yaml --custom-css-file ./config/multiqc_custom_css.css -o output/ output/