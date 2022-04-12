# File name: AFPAP_mutations.py
# Description: Python script for point mutations effect on protein stability
# Author: Maghiar Octavian
# Date: 12-04-2022
import argparse
import logging
import pathlib
from simba2 import methods as sb
from Bio.PDB.Polypeptide import one_to_three
from Bio import BiopythonWarning
import warnings


def generateMutationColor(score):
    limit = 8
    if score < 0:
        offset = (limit+score)/limit if score > -limit else 0
        return(f'rgb(255,{offset*255},{offset*255})')
    else:
        offset = score/limit if score < limit else 1
        return(f'rgb({255*(1-offset)},255,{255*(1-offset)})')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-i", "--input", help="Input file")

    args = parser.parse_args()

    consoleLogger = logging.StreamHandler()
    consoleLogger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    fileLogger = logging.FileHandler(f"{args.outputDir}/workflow.log", mode='a')
    fileLogger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)-8.8s - %(message)s",
        handlers=[
            fileLogger,
            consoleLogger
        ]
    )

    logging.info("Protein stability...")

    aa = ["A", "R", "N", "D", "C", "E", "Q",  "G", "H",  "I",  "L",  "K",  "M",  "F",  "P",   "S", "T",  "W",  "Y",  "V", ]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        result = sb.simba2_predict('PDB', args.input)
        with open(f"{args.AFPAPpath}/config/mutation_template.html", 'r') as template:
            templateData = template.read()
            header = ""
            for i, a in enumerate(result[result["Mutated"] == 'A']["Wild"].to_list()):
                header += f"<th><span data-toggle='tooltip' title='{one_to_three(a)}{i+1}'>{a}</span></th>\n"
            templateData = templateData.replace("---tableHeader---", header)

            IBbody = ""
            SYMbody = ""
            for a in aa:
                IBbody += f"<tr><td>{a}</td>"
                SYMbody += f"<tr><td>{a}</td>"
                df = result[result["Mutated"] == a]
                for i, row in df.iterrows():
                    IBbody += f"<td style='background-color:{generateMutationColor(row['ddG_SimBa_IB'])}'><span data-toggle='tooltip' title='{one_to_three(row['Wild'])}{i+1} to {one_to_three(a)}'>{round(row['ddG_SimBa_IB'],3)}</span></td>"
                    SYMbody += f"<td style='background-color:{generateMutationColor(row['ddG_SimBa_SYM'])}'><span data-toggle='tooltip' title='{one_to_three(row['Wild'])}{i+1} to {one_to_three(a)}'>{round(row['ddG_SimBa_SYM'],3)}</span></td>"
                IBbody += "</tr>"
                SYMbody += "</tr>"

            templateData = templateData.replace("---tableSimBaIB---", IBbody)
            templateData = templateData.replace("---tableSimBaSYM---", SYMbody)
            with open(f'{args.outputDir}/work/multiqc_files/mutation_mqc.html', 'w') as f:
                print(templateData, file=f)
        result.to_csv(f"{args.outputDir}/work/SimBa_predictions.csv", index=False)


if __name__ == '__main__':
    main()
