'''
# File name: AFPAP_structure_analysis.py
# Description: Python script for generating iCn3D viewer section
# Author: Maghiar Octavian
# Date: 30-06-2022
'''
import argparse
import logging
import pathlib


def main():
    '''
    Python script for generating MultiQC report iCn3D viewer section.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("--pdb_type", type=int, help="AlpahFold prediction PDB")

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

    logging.info("Generating Structure Viewer...")
    with open(args.input, 'r', encoding="utf8") as pdb_file:
        pdb_data = pdb_file.read()
        with open(f"{args.AFPAPpath}/config/pdbViewer_template.html", 'r', encoding="utf8") as template:
            template_data = template.read()
            template_data = template_data.replace("--pdb--", pdb_data)
            template_data = template_data.replace("--pdb_type--", 'alphafold' if int(args.pdb_type) == 0 else 'experimental')

            with open(f'{args.outputDir}/work/multiqc_files/structure_viewer_mqc.html', 'w', encoding="utf8") as mqc_file:
                mqc_file.write(template_data)
    logging.info("Structure Viewer completed.")
    return 0


if __name__ == '__main__':
    main()
