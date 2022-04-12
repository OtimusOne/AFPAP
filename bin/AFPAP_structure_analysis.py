# File name: AFPAP_structure_analysis.py
# Description: Python script for generating iCn3D viewer section
# Author: Maghiar Octavian
# Date: 04-04-2022
import argparse
import logging
import pathlib


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

    logging.info("Structure Viewer...")
    with open(args.input, 'r') as pdbFile:
        pdbData = pdbFile.read()
        with open(f"{args.AFPAPpath}/config/pdbViewer_template.html", 'r') as template:
            templateData = template.read()
            templateData = templateData.replace("--pdb--", pdbData)
            with open(f'{args.outputDir}/work/multiqc_files/struct_viewer_mqc.html', 'w') as f:
                f.write(templateData)


if __name__ == '__main__':
    main()
