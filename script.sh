#! /bin/bash

AFPAP_PATH="/mnt/x/Disertatie/AFPAP"
Pfam_PATH="/mnt/x/Disertatie/Pfam"
md_exhaustiveness=4

while getopts i:p:l:o: flag
do
    case "${flag}" in
        i) inputFile=${OPTARG};;
        p) pdbFile=${OPTARG};;
        l) ligandFile=${OPTARG};;
        o) outputDir=${OPTARG};;
    esac
done
if [ ! -f "$inputFile" ] || [ ! -f "$pdbFile" ] ; then
    echo "No arguments supplied"
    exit 1
fi
outputDir="${outputDir:="./output"}"
echo "Input: $inputFile";
echo "PDB: $pdbFile";
echo "Output dir: $outputDir";

baseDir="$PWD"
mkdir -p "$outputDir/work/docking"
mkdir -p "$outputDir/work/multiqc_files"

python "$AFPAP_PATH/bin/AFPAP_sequence_analysis.py" -i $inputFile -o $outputDir --AFPAPpath $AFPAP_PATH

pfam_scan.pl -clan_overlap -align -json pretty -fasta $inputFile -dir $Pfam_PATH > "$outputDir/work/pfam.json"
python "$AFPAP_PATH/bin/AFPAP_pfam.py" -j "$outputDir/work/pfam.json" -o $outputDir --AFPAPpath $AFPAP_PATH

python "$AFPAP_PATH/bin/AFPAP_secondary_structure.py" -i $pdbFile -o $outputDir --AFPAPpath $AFPAP_PATH
python "$AFPAP_PATH/bin/AFPAP_structure_analysis.py" -i "$outputDir/work/proteinStructure_ss.pdb" -o $outputDir --AFPAPpath $AFPAP_PATH

prank predict -f "$outputDir/work/proteinStructure_ss.pdb" -o "$outputDir/work" -c alphafold

python "$AFPAP_PATH/bin/AFPAP_p2rank_visualization.py" -p "$outputDir/work/visualizations/proteinStructure_ss.pdb.pml" -c "$outputDir/work/proteinStructure_ss.pdb_predictions.csv" -o $outputDir --AFPAPpath $AFPAP_PATH

cd "$outputDir/work/visualizations"
pymol -cq proteinStructure_ss.pdb.pml
cd $baseDir
python "$AFPAP_PATH/bin/AFPAP_p2rank_gallery.py"  -o $outputDir --AFPAPpath $AFPAP_PATH 

if [ -f "$ligandFile" ] ; then
    cp "$outputDir/work/proteinStructure_ss.pdb" "$outputDir/work/docking/receptor.pdb"
    ligand_base_name=$(basename ${ligandFile})    
    cp $ligandFile "$outputDir/work/docking/"
    cd "$outputDir/work/docking"
    prepare_receptor -r receptor.pdb -o receptor.pdbqt -A "hydrogens"
    prepare_ligand -l "$outputDir/work/docking/$ligand_base_name" -o ligand.pdbqt
    cd $baseDir
    python "$AFPAP_PATH/bin/AFPAP_molecular_docking.py" -r "$outputDir/work/docking/receptor.pdbqt" -l "$outputDir/work/docking/ligand.pdbqt" -n $ligand_base_name -e $md_exhaustiveness -o $outputDir --AFPAPpath $AFPAP_PATH
fi

multiqc -f -c "$AFPAP_PATH/config/multiqc_config.yaml" --custom-css-file "$AFPAP_PATH/config/multiqc_custom_css.css" -o $outputDir $outputDir