'''
# File name: AFPAP_FASTA_from_PDB.py
# Description: Python script for extracting sequence from PDB
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import logging
import pathlib
import warnings
from Bio import SeqIO
from Bio import BiopythonWarning


def main():
    '''
    Python script for extracting sequence from PDB
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-n", "--name", default='protein', help="PDB name")

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

    logging.info("Sequence extraction...")
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        with open(f'{args.outputDir}/{args.name}', 'w', encoding="utf8") as fasta_file:
            for record in SeqIO.parse(args.input, 'pdb-atom'):
                print('>' + record.id, file=fasta_file)
                print(record.seq, file=fasta_file)


if __name__ == '__main__':
    main()
