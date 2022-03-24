import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
parser.add_argument("--template", help="Template")
args = parser.parse_args()

print("Structure analysis...")

with open(args.input, 'r') as pdbFile:
    pdbData = pdbFile.read()
    with open(args.template, 'r') as template:
        templateData = template.read()
        templateData = templateData.replace("--pdb--", pdbData)
        with open('./output/work/struct_viewer_mqc.html', 'w') as f:
            f.write(templateData)
