import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pymol_script", help="PyMol script")
parser.add_argument("-c", "--csv", help="CSV")
args = parser.parse_args()

with open(args.pymol_script, 'a') as file:
    with open("config/pymol.pml", 'r') as pymol:
        pymolData = pymol.read()
        file.write(pymolData)

df = pd.read_csv(args.csv)
df = df.rename(columns=lambda x: x.strip())
df = df.drop(["residue_ids","surf_atom_ids"],axis=1)
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
df.columns=["Name","Rank","Score","Probability","SAS points", "Surface atoms","Center X", "Center Y", "Center Z"]
df.to_csv("./output/work/p2rank_predictions.csv", index=False)
