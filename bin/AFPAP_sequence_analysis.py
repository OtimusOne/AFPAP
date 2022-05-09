'''
# File name: AFPAP_sequence_analysis.py
# Description: Python script for generating MultiQC report sequence analysis section.
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import logging
import os
import pathlib
import numpy as np
from Bio import SeqIO
from Bio.PDB.Polypeptide import one_to_three
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData


def section_protein_sequence_viewer(protein_sequence, args):
    '''
    Generates MultiQC report sequence viewer.

    Parameters:
        protein_sequence: amino acids sequence string.
    '''
    residue_color_classes = {"GAST": "LeskSmallNonpolar", "CVILPFYMW": "LeskHydrophobic",
                             "NQH": "LeskPolar", "DE": "LeskNegativeCharge", "KR": "LeskPositiveCharge"}
    residue_properties = {
        "A": ["Alanine", "Small", "Hydrophobic", "Non-polar"],
        "R": ["Arginine", "Hydrophilic", "Polar positive"],
        "N": ["Asparagine", "Small", "Hydrophilic", "Polar"],
        "D": ["Aspartic acid ", "Small", "Hydrophilic", "Polar negative"],
        "C": ["Cysteine", "Small", "Hydrophobic", "Special"],
        "E": ["Glutamic acid", "Hydrophilic", "Polar negative"],
        "Q": ["Glutamine", "Hydrophilic", "Polar"],
        "G": ["Glycine", "Small", "Non-polar"],
        "H": ["Histidine", "Polar positive"],
        "I": ["Isoleucine", "Hydrophobic", "Non-polar"],
        "L": ["Leucine", "Hydrophobic", "Non-polar"],
        "K": ["Lysine", "Hydrophilic", "Polar positive"],
        "M": ["Methionine", "Hydrophobic", "Non-polar"],
        "F": ["Phenylalanine", "Hydrophobic", "Non-polar"],
        "P": ["Proline", "Small", "Non-polar"],
        "S": ["Serine", "Small", "Polar"],
        "T": ["Thereonine", "Small", "Polar"],
        "W": ["Tryptophan", "Hydrophobic", "Non-polar"],
        "Y": ["Tyrosine", "Polar"],
        "V": ["Valine", "Hydrophobic", "Non-polar"],
        "X": [],
    }
    with open(f'{args.outputDir}/work/multiqc_files/seq_prot_mqc.html', 'w', encoding="utf8") as mqc_file:
        seq_split = 50
        residue_number = 1
        sequence_rows = "\n"
        for i in range(0, len(protein_sequence), seq_split):
            sequence_rows += "<code>"
            for residue in protein_sequence[i:i+seq_split]:
                residue_class = ""
                for key, value in residue_color_classes.items():
                    if residue in key:
                        residue_class = f"class='sequenceResidue {value}'"
                        break
                residue_description = '\n'.join(residue_properties[residue])
                sequence_rows += f"<span {residue_class}>{residue}<span class='sequenceResidueTooltip'>{one_to_three(residue)}{residue_number}<div class='sequenceResidueTooltipDescription'>{residue_description}</div></span></span>"
                residue_number += 1
            sequence_rows += "</code>\n"

        with open(f'{args.AFPAPpath}/config/sequenceViewer_template.html', 'r', encoding="utf8") as template:
            template_data = template.read()
            template_data = template_data.replace("--sequence--", sequence_rows)
            print(template_data, file=mqc_file)


def section_sequence_properties(protein_sequence, args):
    '''
    Generates MultiQC report sequence properties.

    Parameters:
        protein_sequence: amino acids sequence string.
    '''

    residue_atom_count = {
        "A": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 0},
        "C": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 1},
        "D": {"C": 4, "H": 7, "N": 1, "O": 4, "S": 0},
        "E": {"C": 5, "H": 9, "N": 1, "O": 4, "S": 0},
        "F": {"C": 9, "H": 11, "N": 1, "O": 2, "S": 0},
        "G": {"C": 2, "H": 5, "N": 1, "O": 2, "S": 0},
        "H": {"C": 6, "H": 9, "N": 3, "O": 2, "S": 0},
        "I": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
        "K": {"C": 6, "H": 14, "N": 2, "O": 2, "S": 0},
        "L": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
        "M": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 1},
        "N": {"C": 4, "H": 8, "N": 2, "O": 3, "S": 0},
        "P": {"C": 5, "H": 9, "N": 1, "O": 2, "S": 0},
        "Q": {"C": 5, "H": 10, "N": 2, "O": 3, "S": 0},
        "R": {"C": 6, "H": 14, "N": 4, "O": 2, "S": 0},
        "S": {"C": 3, "H": 7, "N": 1, "O": 3, "S": 0},
        "T": {"C": 4, "H": 9, "N": 1, "O": 3, "S": 0},
        "V": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 0},
        "W": {"C": 11, "H": 12, "N": 2, "O": 2, "S": 0},
        "Y": {"C": 9, "H": 11, "N": 1, "O": 3, "S": 0},
    }
    atom_count = {"C": 0, "H": 2, "N": 0, "O": 1, "S": 0}
    window_scale = 9

    prot_param = ProteinAnalysis(protein_sequence)
    aa_dict = prot_param.count_amino_acids()
    for amino_acid in aa_dict:
        atom_count["C"] += residue_atom_count[amino_acid]["C"]*aa_dict[amino_acid]
        atom_count["H"] += (residue_atom_count[amino_acid]["H"]-2)*aa_dict[amino_acid]
        atom_count["N"] += residue_atom_count[amino_acid]["N"]*aa_dict[amino_acid]
        atom_count["O"] += (residue_atom_count[amino_acid]["O"]-1)*aa_dict[amino_acid]
        atom_count["S"] += residue_atom_count[amino_acid]["S"]*aa_dict[amino_acid]

    with open(f'{args.outputDir}/work/multiqc_files/seq_stats.txt', 'w', encoding="utf8") as mqc_file:

        negative_charge = protein_sequence.count('D')+protein_sequence.count('E')
        positive_charge = protein_sequence.count('R')+protein_sequence.count('K')
        print("Input", "Molecular weight", "Atoms", "Aromaticity", "Instability index", "Isoelctric point",
              "Extinction coefficients", "GRAVY", "Negatively charged residues", "Positively charged residues",  sep='\t', file=mqc_file)
        print(os.path.basename(args.input), prot_param.molecular_weight(), atom_count["C"]+atom_count["H"]+atom_count["N"]+atom_count["O"]+atom_count["S"], prot_param.aromaticity(), prot_param.instability_index(), prot_param.isoelectric_point(
        ), f"{prot_param.molar_extinction_coefficient()[0]} : {prot_param.molar_extinction_coefficient()[1]}", prot_param.gravy(), negative_charge, positive_charge, sep='\t', file=mqc_file)

    with open(f'{args.outputDir}/work/multiqc_files/seq_flexibility.txt', 'w', encoding="utf8") as mqc_file:
        for i, flexibility in zip([x for x in range(1, len(protein_sequence)+1)], prot_param.protein_scale(window=window_scale, param_dict=ProtParamData.Flex)):
            print(f"{i+window_scale//2}\t{flexibility}", file=mqc_file)

    with open(f'{args.outputDir}/work/multiqc_files/seq_hydrophobicity.txt', 'w', encoding="utf8") as mqc_file:
        for i, hydropathy in zip([x for x in range(1, len(protein_sequence)+1)], prot_param.protein_scale(window=window_scale, param_dict=ProtParamData.kd)):
            print(f"{i+window_scale//2}\t{hydropathy}", file=mqc_file)

    with open(f'{args.outputDir}/work/multiqc_files/seq_ph.txt', 'w', encoding="utf8") as mqc_file:
        for i in np.arange(1, 14.01, 0.05):
            print(i, prot_param.charge_at_pH(i), sep="\t", file=mqc_file)

    with open(f'{args.outputDir}/work/multiqc_files/seq_atom_count.txt', 'w', encoding="utf8") as mqc_file:
        atom_sting = "Atom\t"
        count_string = "Count\t"
        for atom, count in atom_count.items():
            atom_sting += f"{atom}\t"
            count_string += f"{count}\t"
        print(atom_sting, file=mqc_file)
        print(count_string, file=mqc_file)

    with open(f'{args.outputDir}/work/multiqc_files/seq_residue_count.txt', 'w', encoding="utf8") as mqc_file:
        print("Residue", end="\t", file=mqc_file)
        for amino_acid in aa_dict:
            print(amino_acid, sep='\t', end="\t", file=mqc_file)
        print("", file=mqc_file)
        print("Residue", end="\t", file=mqc_file)
        for amino_acid in aa_dict:
            print(aa_dict[amino_acid], sep='\t', end="\t", file=mqc_file)


def main():
    '''
    Python script for generating MultiQC report sequence analysis section.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-i", "--input", help="Input file", required=True)

    args = parser.parse_args()

    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    file_logger = logging.FileHandler(f"{args.outputDir}/workflow.log", mode='a')
    file_logger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)-8.8s - %(message)s",
        handlers=[
            file_logger,
            console_logger
        ]
    )

    logging.info("Sequence analysis...")

    records = list(SeqIO.parse(args.input, "fasta"))
    if len(records) == 0:
        logging.error("No records in input fasta file")
        raise AssertionError()
    protein_sequence = str(records[0].seq)

    section_protein_sequence_viewer(protein_sequence, args)
    section_sequence_properties(protein_sequence, args)


if __name__ == '__main__':
    main()
