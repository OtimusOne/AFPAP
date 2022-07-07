'''
# File name: AFPAP_validate_FASTA.py
# Description: Python script for validating the protein sequence
# Author: Maghiar Octavian
# Date: 30-06-2022
'''
import argparse
import logging
import pathlib
import warnings
from Bio import SeqIO
from Bio import BiopythonWarning


def main():
    '''
    Python script for validating the protein sequence
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-n", "--name", default='protein', help="Protein name")
    parser.add_argument("-t", "--type", default=0, help="Input file type 0-FASTA, 1-PDB")

    args = parser.parse_args()

    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    file_logger = logging.FileHandler(f"{args.outputDir}/workflow.log", mode='a')
    file_logger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s - %(message)s",
        handlers=[
            file_logger,
            console_logger
        ]
    )

    logging.info("Extracting protein sequences from input file...")
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        with open(f'{args.outputDir}/{args.name}.fasta', 'w', encoding="utf8") as fasta_file:
            records = list(SeqIO.parse(args.input, "pdb-atom" if int(args.type) == 1 else "fasta"))
            if len(records) > 1:
                logging.warning("Input file contains multiple sequences, extracting one...")
            if len(records):
                if 'X' in str(records[0].seq).upper():
                    logging.error("Unknown residue in protein sequence!")
                    raise Exception("Unknown residue in protein sequence!")
                AA_ALPHABET = frozenset("ARNDCEQGHILKMFPSTWYV")
                if (not set(str(records[0].seq).upper()) <= AA_ALPHABET):
                    logging.error("Protein sequence contains characters different than the 20 standard amino-acids!")
                    raise Exception("Protein sequence contains characters different than the 20 standard amino-acids!")
                print('>' + str(records[0].id).replace("????", args.name), file=fasta_file)
                print(str(records[0].seq).upper(), file=fasta_file)
                logging.info("Protein sequence extraction completed.")
            else:
                logging.error("No sequences found in input file!")
                raise Exception("No sequences found in input file!")
    return 0


if __name__ == '__main__':
    main()
