'''
# File name: AFPAP_molecular_docking.py
# Description: Python script for molecular docking with AutoDock Vina
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import logging
import pathlib
import os.path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from vina import Vina


def dock(receptor=None, flexible=None, ligand=None, center=(0, 0, 0), box_size=(20, 20, 20), spacing=0.375, exhaustiveness=32, n_poses=20):
    '''
    Docks the ligand to the receptor using AutoDock Vina.

    Parameters:
        receptor: Rigid receptor pdbqt file.
        flexible: Flexible receptor pdbqt file.
        center: Docking center x,y,z coordinates.
        box_size: Docking grid box x,y,z points.
        spacing: Docking spacing between grid points.
        exhaustiveness: Numbers of Monte Carlo runs.
        n_poses: Number of poses to generate.
    Returns:
        Vina object with docked poses.
    '''
    vina = Vina(sf_name='vina')
    vina.set_receptor(rigid_pdbqt_filename=receptor, flex_pdbqt_filename=flexible)
    vina.set_ligand_from_file(ligand)
    vina.compute_vina_maps(center=center, box_size=box_size, spacing=spacing)
    vina.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    return vina


def calculate_wide_box(structure_file=None):
    '''
    Calculate wide box encompassing molecule.

    Parameters:
        structure_file: pdb or pdbqt file.
    Returns:
        x center coordinate, x axis points, y center coordinate, y axis points, z center coordinate, z axis points
    '''
    pdb_object = PDBParser(QUIET=True)
    structure = pdb_object.get_structure("APLHA", structure_file)
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
    '''
    Python script for molecular docking with AutoDock Vina
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")
    parser.add_argument("-r", "--receptor", required=True, help="Receptor")
    parser.add_argument("-l", "--ligand", required=True, help="Ligand")
    parser.add_argument("-n", "--name", default='ligand', help="Ligand name")
    parser.add_argument("-e", "--exhaustiveness", type=int, help="Exhaustiveness", default=8)
    parser.add_argument("--box_size", type=int, help="Box size", default=20)
    parser.add_argument("--spacing", type=float, help="Spacing", default=0.375)
    parser.add_argument('--dock_pockets', action="count", help="Dock ligand to predicted pockets")

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
    md_exhaustiveness = int(args.exhaustiveness)
    md_box_size = int(args.box_size)
    md_spacing = float(args.spacing)

    logging.info("Molecular docking ligand %s...", args.name)
    with open(f'{args.outputDir}/work/multiqc_files/molecular_docking_{args.name}_mqc.txt', 'w', encoding="utf8") as mqc_file:
        with open(f'{args.AFPAPpath}/config/molecularDocking_template.txt', 'r', encoding="utf8") as template:
            template_data = template.read()
            template_data = template_data.replace('--ligand--', args.name)
            print(template_data, file=mqc_file)

        with open(f'{args.outputDir}/work/docking/{args.name}/generate_complex.pml', 'a', encoding="utf8") as pml:
            print("set retain_order\nset pdb_retain_ids", file=pml)

            if args.dock_pockets and os.path.exists(f"{args.outputDir}/work/p2rank_predictions.csv"):
                logging.info("%s - Pocket docking...", args.name)
                p2rank_precitions = pd.read_csv(f"{args.outputDir}/work/p2rank_predictions.csv")
                for i, row in p2rank_precitions.iterrows():
                    logging.info("%s - Docking pocket %d/%d...", args.name, i+1, p2rank_precitions.shape[0])

                    vina_object = dock(receptor=args.receptor, ligand=args.ligand, center=[row['Center X'], row['Center Y'], row['Center Z']], box_size=[
                                       md_box_size, md_box_size, md_box_size], spacing=md_spacing, exhaustiveness=md_exhaustiveness)

                    energies = vina_object.energies()[:, 0]
                    best_score, mean_score, std_score = np.round(energies[0], 4), np.round(np.mean(energies), 4), np.round(np.std(energies), 4)
                    vina_object.write_poses(f"{args.outputDir}/work/docking/{args.name}/{args.name}_pocket{i+1}_best_pose.pdbqt", overwrite=True, n_poses=1)
                    vina_object.write_poses(f"{args.outputDir}/work/docking/{args.name}/{args.name}_pocket{i+1}.pdbqt", overwrite=True)
                    print(f'load "{args.name}_pocket{i+1}_best_pose.pdbqt", ligand', file=pml)
                    print('load "../receptor.pdbqt", receptor', file=pml)
                    print(f'save {args.name}_pocket{i+1}_complex.pdb, state=1', file=pml)
                    print('delete all', file=pml)

                    print(row['Name'], best_score, f"{mean_score} \u00B1 {std_score}",
                          f"{np.round(row['Center X'],2)},{np.round(row['Center Y'],2)},{np.round(row['Center Z'],2)}", f"{md_box_size},{md_box_size},{md_box_size}", md_spacing, md_exhaustiveness, sep='\t', file=mqc_file)

            logging.info("%s - Blind docking...", args.name)

            center_x, box_x, center_y, box_y, center_z, box_z = calculate_wide_box(args.receptor)
            vina_object = dock(receptor=args.receptor, ligand=args.ligand, center=[center_x, center_y, center_z], box_size=[box_x, box_y, box_z], spacing=1, exhaustiveness=md_exhaustiveness)

            energies = vina_object.energies()[:, 0]
            best_score, mean_score, std_score = np.round(energies[0], 4), np.round(np.mean(energies), 4), np.round(np.std(energies), 4)
            vina_object.write_poses(f"{args.outputDir}/work/docking/{args.name}/{args.name}_blind_docking_best_pose.pdbqt", overwrite=True, n_poses=1)
            vina_object.write_poses(f"{args.outputDir}/work/docking/{args.name}/{args.name}_blind_docking.pdbqt", overwrite=True)
            print(f'load "{args.name}_blind_docking_best_pose.pdbqt", ligand', file=pml)
            print('load "../receptor.pdbqt", receptor', file=pml)
            print(f'save {args.name}_blind_docking_complex.pdb, state=1', file=pml)
            print('delete all', file=pml)

            print('Blind docking', best_score, f"{mean_score} \u00B1 {std_score}", f"{np.round(center_x,2)},{np.round(center_y,2)},{np.round(center_z,2)}",
                  f"{int(box_x)},{int(box_y)},{int(box_z)}", 1, md_exhaustiveness, sep='\t', file=mqc_file)


if __name__ == '__main__':
    main()
