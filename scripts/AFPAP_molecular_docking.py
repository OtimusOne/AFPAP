from vina import Vina
import numpy as np
import argparse
import pandas as pd
from Bio.PDB import PDBParser

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--receptor", help="Receptor")
parser.add_argument("-l", "--ligand", help="Ligand")
parser.add_argument("-n", "--name", help="Ligand name")
args = parser.parse_args()


def dock(receptor=None, flexible=None, ligand=None, center=[0, 0, 0], box_size=[20, 20, 20], spacing=0.375, exhaustiveness=32, n_poses=20):
    v = Vina(sf_name='vina')
    v.set_receptor(rigid_pdbqt_filename=receptor, flex_pdbqt_filename=flexible)

    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=center, box_size=box_size, spacing=spacing)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    en = v.energies()[:, 0]
    return(np.round(en[0], 4), np.round(np.mean(en), 4), np.round(np.std(en), 4))


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


with open(f'./output/work/molecular_docking_{args.name}_mqc.txt', 'w') as f:
    with open('config/molecularDocking_template.txt', 'r') as template:
        exh = 32
        templateData = template.read()
        templateData = templateData.replace('--ligand--', args.name)
        print(templateData, file=f)
        df = pd.read_csv("./output/work/p2rank_predictions.csv")
        for i, row in df.iterrows():
            bestScore, meanScore, stdScore = dock(receptor=args.receptor, ligand=args.ligand, center=[
                                                  row['Center X'], row['Center Y'], row['Center Z']], box_size=[20, 20, 20], spacing=0.375, exhaustiveness=exh)
            print(row['Name'], bestScore, f"{meanScore} \u00B1 {stdScore}",
                  f"{np.round(row['Center X'],2)},{np.round(row['Center Y'],2)},{np.round(row['Center Z'],2)}", "20,20,20", 0.375, exh, sep='\t', file=f)
        center_x, box_x, center_y, box_y, center_z, box_z = wideBox(
            args.receptor)
        bestScore, meanScore, stdScore = dock(receptor=args.receptor, ligand=args.ligand, center=[
                                              center_x, center_y, center_z], box_size=[box_x, box_y, box_z], spacing=1, exhaustiveness=exh)
        print('Wide Box', bestScore, f"{meanScore} \u00B1 {stdScore}", f"{np.round(center_x,2)},{np.round(center_y,2)},{np.round(center_z,2)}",
              f"{int(box_x)},{int(box_y)},{int(box_z)}", 1, exh, sep='\t', file=f)
