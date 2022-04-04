# File name: AFPAP_p2rank_visualization.py
# Description: Python script for parsing P2Rank output
# Author: Maghiar Octavian
# Date: 04-04-2022
import argparse
import logging
import pathlib
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count",
                        help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output",
                        help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True,
                        help="Path to AFPAP home")
    parser.add_argument("-p", "--pymol_script", help="PyMol script")
    parser.add_argument("-c", "--csv", help="CSV")

    args = parser.parse_args()
    consoleLogger = logging.StreamHandler()
    consoleLogger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    fileLogger = logging.FileHandler(
        f"{args.outputDir}/workflow.log", mode='a')
    fileLogger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)-8.8s - %(message)s",
        handlers=[
            fileLogger,
            consoleLogger
        ]
    )
    logging.info("P2Rank visualization...")
    with open(args.pymol_script, 'a') as file:
        with open(f"{args.AFPAPpath}/config/pymol.pml", 'r') as pymol:
            pymolData = pymol.read()
            file.write(pymolData)

    df = pd.read_csv(args.csv)
    df = df.rename(columns=lambda x: x.strip())
    df = df.drop(["residue_ids", "surf_atom_ids"], axis=1)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    df.columns = ["Name", "Rank", "Score", "Probability", "SAS points",
                  "Surface atoms", "Center X", "Center Y", "Center Z"]
    df.to_csv(f"{args.outputDir}/work/p2rank_predictions.csv", index=False)


if __name__ == '__main__':
    main()
