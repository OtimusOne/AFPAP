#! /bin/bash

AFPAP_PATH="/mnt/x/Disertatie/AFPAP"
Pfam_PATH="/mnt/x/Disertatie/Pfam"
md_exhaustiveness=8

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

python "$AFPAP_PATH/bin/AFPAP_sequence_analysis.py" -i $inputFile --template "$AFPAP_PATH/config/sequenceViewer_template.html"

pfam_scan.pl -clan_overlap -align -json pretty -fasta $inputFile -dir $Pfam_PATH > ./output/work/pfam.json
python "$AFPAP_PATH/bin/AFPAP_pfam.py" -j output/work/pfam.json --template "$AFPAP_PATH/config/pfam_template.html"

python "$AFPAP_PATH/bin/AFPAP_secondary_structure.py" -i $pdbFile
python "$AFPAP_PATH/bin/AFPAP_structure_analysis.py" -i ./output/work/proteinStructure_ss.pdb --template "$AFPAP_PATH/config/pdbViewer_template.html"

prank predict -f output/work/proteinStructure_ss.pdb -o output/work/ -c alphafold

python "$AFPAP_PATH/bin/AFPAP_p2rank_visualization.py" -p output/work/visualizations/proteinStructure_ss.pdb.pml -c output/work/proteinStructure_ss.pdb_predictions.csv --template "$AFPAP_PATH/config/pymol.pml"

cd output/work/visualizations
pymol -cq proteinStructure_ss.pdb.pml
cd $baseDir
python "$AFPAP_PATH/bin/AFPAP_p2rank_gallery.py" --template "$AFPAP_PATH/config/pocketViewer_template.html"

if [ -f "$ligandFile" ] ; then
    cp ./output/work/proteinStructure_ss.pdb ./output/work/docking/receptor.pdb
    ligand_base_name=$(basename ${ligandFile})    
    cp $ligandFile ./output/work/docking/
    cd ./output/work/docking
    prepare_receptor -r receptor.pdb -o receptor.pdbqt -A "hydrogens"
    prepare_ligand -l ./output/work/docking/$ligand_base_name -o ligand.pdbqt
    cd $baseDir
    python "$AFPAP_PATH/bin/AFPAP_molecular_docking.py" -r ./output/work/docking/receptor.pdbqt -l ./output/work/docking/ligand.pdbqt -n $ligand_base_name -e $md_exhaustiveness --template "$AFPAP_PATH/config/molecularDocking_template.txt"
fi

multiqc -f -c "$AFPAP_PATH/config/multiqc_config.yaml" --custom-css-file "$AFPAP_PATH/config/multiqc_custom_css.css" -o output/ output/