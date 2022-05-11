'''
# File name: AFPAP_p2rank_visualization.py
# Description: Python script for parsing P2Rank output
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import logging
import pathlib
import pandas as pd


def main():
    '''
    Python script for parsing P2Rank output.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-p", "--pymol_script", help="PyMol script")
    parser.add_argument("--csv", help="P2Rank prediction CSV file")

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
    logging.info("P2Rank visualization...")
    with open(args.pymol_script, 'a', encoding="utf8") as pymol_script:
        with open(f"{args.AFPAPpath}/config/pymol.pml", 'r', encoding="utf8") as pymol_template:
            pymol_data = pymol_template.read()
            pymol_script.write(pymol_data)

    p2rank_predictions = pd.read_csv(args.csv)
    p2rank_predictions = p2rank_predictions.rename(columns=lambda x: x.strip())
    p2rank_predictions = p2rank_predictions.drop(["residue_ids", "surf_atom_ids"], axis=1)
    p2rank_predictions = p2rank_predictions.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    p2rank_predictions.columns = ["Name", "Rank", "Score", "Probability", "SAS points", "Surface atoms", "Center X", "Center Y", "Center Z"]
    p2rank_predictions.to_csv(f"{args.outputDir}/work/p2rank_predictions.csv", index=False)


if __name__ == '__main__':
    main()
