# File name: AFPAP_molecular_docking.py
# Description: Python script for AutoDock Vina
# Author: Maghiar Octavian
# Date: 04-04-2022
import argparse
import logging
import pathlib
import pandas as pd
import numpy as np
from vina import Vina
from Bio.PDB import PDBParser


def dock(receptor=None, flexible=None, ligand=None, center=[0, 0, 0], box_size=[20, 20, 20], spacing=0.375, exhaustiveness=32, n_poses=20):
    v = Vina(sf_name='vina')
    v.set_receptor(rigid_pdbqt_filename=receptor, flex_pdbqt_filename=flexible)

    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=center, box_size=box_size, spacing=spacing)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    return v


def wideBox(structureFile):
    p = PDBParser(QUIET=True)
    structure = p.get_structure("APLHA", structureFile)
    model = structure[0]
    coords = np.array([atom.get_coord() for atom in model.get_atoms()])
    x_min = np.min(coords[:, 0])
    x_max = np.max(coords[:, 0])
    y_min = np.min(coords[:, 1])
    y_max = np.max(coords[:, 1])
    z_min = np.min(coords[:, 2])
    z_max = np.max(coords[:, 2])
    x_box = np.ceil(np.sqrt((x_min-x_max)**2))
    y_box = np.ceil(np.sqrt((y_min-y_max)**2))
    z_box = np.ceil(np.sqrt((z_min-z_max)**2))
    # print(x_min,x_max,(x_min+x_max)/2,x_box)
    # print(y_min,y_max,(y_min+y_max)/2,y_box)
    # print(z_min,z_max,(z_min+z_max)/2,z_box)
    return((x_min+x_max)/2, x_box, (y_min+y_max)/2, y_box, (z_min+z_max)/2, z_box)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count",
                        help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output",
                        help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True,
                        help="Path to AFPAP home")
    parser.add_argument("-r", "--receptor", help="Receptor")
    parser.add_argument("-l", "--ligand", help="Ligand")
    parser.add_argument("-n", "--name", help="Ligand name")
    parser.add_argument("-e", "--exhaustiveness", help="Exhaustiveness")

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
    logging.info(f"Molecular docking {args.name}...")
    with open(f'{args.outputDir}/work/multiqc_files/molecular_docking_{args.name}_mqc.txt', 'w') as f:
        with open(f'{args.AFPAPpath}/config/molecularDocking_template.txt', 'r') as template:
            exh = int(args.exhaustiveness)
            templateData = template.read()
            templateData = templateData.replace('--ligand--', args.name)
            print(templateData, file=f)
            df = pd.read_csv(f"{args.outputDir}/work/p2rank_predictions.csv")
            for i, row in df.iterrows():
                logging.info(f"Docking pocket {i+1}/{df.shape[0]}...")

                vinaObject = dock(receptor=args.receptor, ligand=args.ligand, center=[
                    row['Center X'], row['Center Y'], row['Center Z']], box_size=[20, 20, 20], spacing=0.375, exhaustiveness=exh)

                energies = vinaObject.energies()[:, 0]
                bestScore, meanScore, stdScore = np.round(energies[0], 4), np.round(
                    np.mean(energies), 4), np.round(np.std(energies), 4)
                vinaObject.write_poses(
                    f"{args.outputDir}/work/docking/ADV_{args.name}_pocket{i+1}_bestPose.pdbqt", overwrite=True, n_poses=1)
                vinaObject.write_poses(
                    f"{args.outputDir}/work/docking/ADV_{args.name}_pocket{i+1}.pdbqt", overwrite=True)
                print(row['Name'], bestScore, f"{meanScore} \u00B1 {stdScore}",
                      f"{np.round(row['Center X'],2)},{np.round(row['Center Y'],2)},{np.round(row['Center Z'],2)}", "20,20,20", 0.375, exh, sep='\t', file=f)

            logging.info("Blind docking...")
            center_x, box_x, center_y, box_y, center_z, box_z = wideBox(
                args.receptor)
            vinaObject = dock(receptor=args.receptor, ligand=args.ligand, center=[
                center_x, center_y, center_z], box_size=[box_x, box_y, box_z], spacing=1, exhaustiveness=exh)
            energies = vinaObject.energies()[:, 0]
            bestScore, meanScore, stdScore = np.round(energies[0], 4), np.round(
                np.mean(energies), 4), np.round(np.std(energies), 4)
            vinaObject.write_poses(
                f"{args.outputDir}/work/docking/ADV_{args.name}_blindDocking_bestPose.pdbqt", overwrite=True, n_poses=1)
            vinaObject.write_poses(
                f"{args.outputDir}/work/docking/ADV_{args.name}_blindDocking.pdbqt", overwrite=True)
            print('Wide Box', bestScore, f"{meanScore} \u00B1 {stdScore}", f"{np.round(center_x,2)},{np.round(center_y,2)},{np.round(center_z,2)}",
                  f"{int(box_x)},{int(box_y)},{int(box_z)}", 1, exh, sep='\t', file=f)


if __name__ == '__main__':
    main()
