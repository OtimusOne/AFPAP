'''
# File name: AFPAP_mutations.py
# Description: Python script for predicting point mutations effect on protein stability
# Author: Maghiar Octavian
# Date: 12-04-2022
'''
import argparse
import logging
import pathlib
import warnings
from Bio import BiopythonWarning
from Bio.PDB.Polypeptide import one_to_three
from simba2 import methods as sb


def generate_mutation_effect_color(score=0):
    '''
    Calculate red-green scpectrum color based on mutation effect score.

    Parameters:
        score: Mutation effect on protein stability.
    Returns:
        rgb(r,g,b) string
    '''
    limit = 8
    if score < 0:
        offset = (limit+score)/limit if score > -limit else 0
        return f'rgb(255,{offset*255},{offset*255})'
    offset = score/limit if score < limit else 1
    return f'rgb({255*(1-offset)},255,{255*(1-offset)})'


def main():
    '''
    Python script for predicting point mutations effect on protein stability.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-i", "--input", help="Input file")

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

    logging.info("Protein stability...")

    residues = ["A", "R", "N", "D", "C", "E", "Q",  "G", "H",  "I",  "L",  "K",  "M",  "F",  "P",   "S", "T",  "W",  "Y",  "V", ]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        mutation_result = sb.simba2_predict('PDB', args.input)
        with open(f"{args.AFPAPpath}/config/mutation_template.html", 'r', encoding="utf8") as template:
            template_data = template.read()
            header = ""
            for i, residue in enumerate(mutation_result[mutation_result["Mutated"] == 'A']["Wild"].to_list()):
                header += f"<th><span data-toggle='tooltip' title='{one_to_three(residue)}{i+1}'>{residue}</span></th>\n"
            template_data = template_data.replace("---tableHeader---", header)

            ib_table_body = ""
            sym_table_body = ""
            for residue in residues:
                ib_table_body += f"<tr><td>{residue}</td>"
                sym_table_body += f"<tr><td>{residue}</td>"
                mutation_result_row = mutation_result[mutation_result["Mutated"] == residue]
                for i, row in mutation_result_row.iterrows():
                    ib_table_body += f"<td style='background-color:{generate_mutation_effect_color(row['ddG_SimBa_IB'])}'><span data-toggle='tooltip' title='{one_to_three(row['Wild'])}{i+1} to {one_to_three(residue)}'>{round(row['ddG_SimBa_IB'],3)}</span></td>"
                    sym_table_body += f"<td style='background-color:{generate_mutation_effect_color(row['ddG_SimBa_SYM'])}'><span data-toggle='tooltip' title='{one_to_three(row['Wild'])}{i+1} to {one_to_three(residue)}'>{round(row['ddG_SimBa_SYM'],3)}</span></td>"
                ib_table_body += "</tr>\n"
                sym_table_body += "</tr>\n"

            template_data = template_data.replace("---tableSimBaIB---", ib_table_body)
            template_data = template_data.replace("---tableSimBaSYM---", sym_table_body)
            with open(f'{args.outputDir}/work/multiqc_files/mutation_mqc.html', 'w', encoding="utf8") as mqc_file:
                print(template_data, file=mqc_file)
        mutation_result.to_csv(f"{args.outputDir}/work/SimBa_predictions.csv", index=False)


if __name__ == '__main__':
    main()
