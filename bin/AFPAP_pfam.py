import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--json", help="Pfam json file")

args = parser.parse_args()

with open(args.json, 'r') as file:
    pfam_data = json.load(file);
    for record in pfam_data:
        print(record["name"],record["acc"],record["clan"],record["desc"],record["evalue"])
