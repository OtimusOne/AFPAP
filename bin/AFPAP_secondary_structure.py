'''
# File name: AFPAP_secondary_structure.py
# Description: Python script for generating secondary structure of pdb file
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import logging
import pathlib
import re
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import one_to_three


def prepend_whitespace(string="", length=1):
    '''
    Prepend whitespace to string.

    Parameters:
        string(str)
        length(int)
    Returns:
        string with added whitespace
    '''
    return ' '*(length-len(string))+string if len(string) < length else string


def get_dssp_df(pdb_file=None):
    '''
    Create DSSP dataframe from pdb file.

    Parameters:
        pdb_file: PDB file.
    Returns:
        Dataframe containing DSSP prediction.
    '''
    pdb_parser = PDBParser()
    pdb_structure = pdb_parser.get_structure("APLHA", pdb_file)
    model = pdb_structure[0]
    dssp = DSSP(model=model, in_file=pdb_file)
    appender = []
    for k in dssp.property_keys:
        to_append = []
        properties = dssp.property_dict[k]
        chain = k[0]
        residue = k[1]
        resnum = residue[1]
        to_append.extend([chain, resnum])
        to_append.extend(properties[1:3])
        appender.append(to_append)

    cols = ['chain', 'resnum', 'aa', 'ss']

    dataframe = pd.DataFrame.from_records(appender, columns=cols)
    dataframe['aa_three'] = dataframe['aa'].apply(one_to_three)
    return dataframe


def main():
    '''
    Python script for generating secondary structure of pdb file
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True,  help="Path to AFPAP home")
    parser.add_argument("-i", "--input", required=True, help="Input file")

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
    logging.info("Secondary structure calculation...")

    with open(args.input, 'r', encoding="utf8") as pdb_file:
        pdb_data = pdb_file.read()
        with open(f"{args.outputDir}/work/proteinStructure.pdb", 'w', encoding="utf8") as pdb_secondary_structure_file:

            if pdb_data.find("HELIX") == -1 and pdb_data.find("SHEET") == -1:

                dssp_df = get_dssp_df(args.input)
                secondary_structure = "".join(dssp_df["ss"].tolist())
                for symbol_pair in (("G", "H"), ("I", "H"), ("S", "-"), ("T", "-"), ("B", "E")):
                    secondary_structure = secondary_structure.replace(*symbol_pair)

                k_helix = 1
                ss_rows = []
                helixes = [match.span() for match in re.finditer('[H]{3,}', secondary_structure)]
                for helix_start, helix_end in helixes:
                    helix_start, helix_end = dssp_df.iloc[helix_start], dssp_df.iloc[helix_end-1]
                    helix_row = f"HELIX  {prepend_whitespace(str(k_helix),3)} {prepend_whitespace(str(k_helix),3)} {prepend_whitespace(str(helix_start['aa_three']),3)} {helix_start['chain']} {prepend_whitespace(str(helix_start['resnum']),4)}  {prepend_whitespace(str(helix_end['aa_three']),3)} {helix_end['chain']} {prepend_whitespace(str(helix_end['resnum']),4)}  1                               {prepend_whitespace(str(helix_end['resnum']-helix_start['resnum']+1),5)}"
                    ss_rows.append(helix_row)
                    k_helix += 1

                k_sheets = 1
                sheets = [match.span() for match in re.finditer('[E]{3,}', secondary_structure)]
                for sheet_start, sheet_end in sheets:
                    sheet_start, sheet_end = dssp_df.iloc[sheet_start], dssp_df.iloc[sheet_end-1]
                    sheet_row = f"SHEET  {prepend_whitespace(str(k_sheets),3)} {prepend_whitespace(str(k_sheets),3)} 1 {prepend_whitespace(str(sheet_start['aa_three']),3)} {sheet_start['chain']}{prepend_whitespace(str(sheet_start['resnum']),4)}  {prepend_whitespace(str(sheet_end['aa_three']),3)} {sheet_end['chain']}{prepend_whitespace(str(sheet_end['resnum']),4)}  0"
                    ss_rows.append(sheet_row)
                    k_sheets += 1

                atom_records = re.search("ATOM      1", pdb_data)
                if atom_records:
                    if atom_records.span()[0] != 0:
                        print(pdb_data[:atom_records.span()[0]], file=pdb_secondary_structure_file, end="")
                    for secondary_structure in ss_rows:
                        print(secondary_structure, file=pdb_secondary_structure_file)
                    print(pdb_data[atom_records.span()[0]:], file=pdb_secondary_structure_file)
                else:
                    print(pdb_data, file=pdb_secondary_structure_file)

            else:
                print(pdb_data, file=pdb_secondary_structure_file)


if __name__ == '__main__':
    main()
