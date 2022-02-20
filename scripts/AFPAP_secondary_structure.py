import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import one_to_three
import argparse
import re 

def prependWhitespace(s,length):
  return ' '*(length-len(s))+s if len(s)<length else s

def get_dssp_df(model, pdb_file):
    dssp = DSSP(model=model, in_file=pdb_file)
    appender = []
    for k in dssp.property_keys:
        to_append = []
        x = dssp.property_dict[k]
        chain = k[0]
        residue = k[1]
        resnum = residue[1]
        to_append.extend([chain, resnum])
        to_append.extend(x[1:3])
        appender.append(to_append)

    cols = ['chain', 'resnum', 'aa', 'ss']

    df = pd.DataFrame.from_records(appender, columns=cols)
    #df = df[df['aa'].isin(list(aa1))]
    df['aa_three'] = df['aa'].apply(one_to_three)
    return df

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
args = parser.parse_args()

p = PDBParser()
structure = p.get_structure("APLHA", args.input)
model = structure[0]
df=get_dssp_df(model, args.input)
#print(df)
ss = "".join(df["ss"].tolist())
for r in (("G", "H"), ("I","H"), ("S","-"), ("T","-"), ("B","E")):
    ss = ss.replace(*r)
#print(ss)


kHelix=1
ssRows=[]
helixes = [m.span() for m in re.finditer('[H]{3,}', ss)]
for helixStart, helixEnd in helixes:
  helixStart,helixEnd = df.iloc[helixStart],df.iloc[helixEnd-1]
  helixRow = f"HELIX  {prependWhitespace(str(kHelix),3)} {prependWhitespace(str(kHelix),3)} {prependWhitespace(str(helixStart['aa_three']),3)} {helixStart['chain']} {prependWhitespace(str(helixStart['resnum']),4)}  {prependWhitespace(str(helixEnd['aa_three']),3)} {helixEnd['chain']} {prependWhitespace(str(helixEnd['resnum']),4)}  1                               {prependWhitespace(str(helixEnd['resnum']-helixStart['resnum']+1),5)}"
  #print(helixRow)
  ssRows.append(helixRow)
  kHelix+=1


kSheets=1
sheets = [m.span() for m in re.finditer('[E]{3,}', ss)]
for sheetStart, sheetEnd in sheets:
  sheetStart,sheetEnd = df.iloc[sheetStart],df.iloc[sheetEnd-1]
  sheetRow = f"SHEET  {prependWhitespace(str(kSheets),3)} {prependWhitespace(str(kSheets),3)} 1 {prependWhitespace(str(sheetStart['aa_three']),3)} {sheetStart['chain']}{prependWhitespace(str(sheetStart['resnum']),4)}  {prependWhitespace(str(sheetEnd['aa_three']),3)} {sheetEnd['chain']}{prependWhitespace(str(sheetEnd['resnum']),4)}  0"
  #print(sheetRow)
  ssRows.append(sheetRow)
  kSheets+=1

with open(args.input, 'r') as pdbFile:
  pdbData = pdbFile.read()
  with open(f'./output/work/{os.path.basename(args.input).split(".")[0]}_secondary_structure.pdb', 'w') as f:
    for ss in ssRows:
      print(ss,file=f)
    print(pdbData,file=f)