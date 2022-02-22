import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from Bio.PDB.Polypeptide import one_to_three

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
args = parser.parse_args()

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
atomCount = {"C": 0, "H": 0, "N": 0, "O": 0, "S": 0}

records = list(SeqIO.parse(args.input, "fasta"))
if not len(records):
    raise AssertionError()
record = records[0]
protein_sequence = str(record.seq)
protParam = ProteinAnalysis(protein_sequence)
AAdict = protParam.count_amino_acids()
for aa in AAdict:
    atomCount["C"] += residue_atom_count[aa]["C"]*AAdict[aa]
    atomCount["H"] += residue_atom_count[aa]["H"]*AAdict[aa]
    atomCount["N"] += residue_atom_count[aa]["N"]*AAdict[aa]
    atomCount["O"] += residue_atom_count[aa]["O"]*AAdict[aa]
    atomCount["S"] += residue_atom_count[aa]["S"]*AAdict[aa]


def section_proteinSequence():
    with open('./output/work/seq_prot_mqc.html', 'w') as f:
        seq_split = 50
        residueColorClass = {"GAST": "smallNonpolar", "CVILPFYMW": "hydrophobic",
                             "NQH": "polar", "DE": "negativeCharge", "KR": "positiveCharge"}
        resNumber = 1
        residueLines = "\n"
        for i in range(0, len(protein_sequence), seq_split):
            residueLines += "<code>"
            for residue in protein_sequence[i:i+seq_split]:
                residueClass = ""
                for key, value in residueColorClass.items():
                    if residue in key:
                        residueClass = f"class='residueTooltip {value}'"
                        break
                residueLines += f"<span {residueClass}>{residue}<span class='residueTooltipText'>{one_to_three(residue)}{resNumber}</span></span>"
                resNumber += 1
            residueLines += "</code>\n"

        with open('config/sequenceViewer_template.html', 'r') as template:
            templateData = template.read()
            templateData = templateData.replace("--sequence--", residueLines)
            print(templateData, file=f)


def section_sequenceProperties():
    with open('./output/work/seq_stats.txt', 'w') as f:

        negativeCharge = protein_sequence.count(
            'D')+protein_sequence.count('E')
        positiveCharge = protein_sequence.count(
            'R')+protein_sequence.count('K')
        print("Input", "Molecular weight", "Atoms", "Aromaticity", "Instability index", "Isoelctric point",
              "Extinction coefficients", "GRAVY", "Negatively charged residues", "Positively charged residues",  sep='\t', file=f)
        print(os.path.basename(args.input), protParam.molecular_weight(), atomCount["C"]+atomCount["H"]+atomCount["N"]+atomCount["O"]+atomCount["S"], protParam.aromaticity(), protParam.instability_index(), protParam.isoelectric_point(
        ), f"{protParam.molar_extinction_coefficient()[0]} : {protParam.molar_extinction_coefficient()[1]}", protParam.gravy(), negativeCharge, positiveCharge, sep='\t', file=f)


def section_sequenceFlexibility(scaleWindow=9):
    with open('./output/work/seq_stats_flex.txt', 'w') as f:

        for x, y in zip([x for x in range(1, len(protein_sequence)+1)], protParam.protein_scale(window=scaleWindow, param_dict=ProtParamData.Flex)):
            print(f"{x+scaleWindow//2}\t{y}", file=f)


def section_sequenceHydrophobicity(scaleWindow=9):
    with open('./output/work/seq_stats_hydro.txt', 'w') as f:
        for x, y in zip([x for x in range(1, len(protein_sequence)+1)], protParam.protein_scale(window=scaleWindow, param_dict=ProtParamData.kd)):
            print(f"{x+scaleWindow//2}\t{y}", file=f)


def section_sequencePH():
    with open('./output/work/seq_stats_ph.txt', 'w') as f:
        for x in np.arange(1, 14.01, 0.05):
            print(x, protParam.charge_at_pH(x), sep="\t", file=f)


def section_atomCount():
    with open('./output/work/seq_stats_atomcount.txt', 'w') as f:

        print("Atom", end="\t", file=f)
        for x in atomCount:
            print(x, end="\t", file=f)

        print("\nCount", end="\t", file=f)
        for x in atomCount:
            print(atomCount[x], end="\t", file=f)


def section_residueCount():
    with open('./output/work/seq_stats_aa.txt', 'w') as f:

        print("Residue", end="\t", file=f)
        for aa in AAdict:
            print(aa, sep='\t', end="\t", file=f)
        print("\n", file=f)
        print("Residue", end="\t", file=f)
        for aa in AAdict:
            print(AAdict[aa], sep='\t', end="\t", file=f)


print("Sequence analysis...")
section_proteinSequence()
section_sequenceProperties()
section_sequenceFlexibility()
section_sequenceHydrophobicity()
section_sequencePH()
section_atomCount()
section_residueCount()
