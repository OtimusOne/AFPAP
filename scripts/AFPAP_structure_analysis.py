import sys
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
args = parser.parse_args()

with open(args.input, 'r') as pdbFile:
    pdbData = pdbFile.read()
    with open('config/pdbViewer_template.html', 'r') as template:
        templateData = template.read()
        templateData = templateData.replace("--pdb--", pdbData)
        with open('./output/work/struct_viewer_mqc.html', 'w') as f:
            f.write(templateData)

