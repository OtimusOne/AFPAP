# File name: AFPAP_pfam.py
# Description: Python script for generating Pfam family alignment
# Author: Maghiar Octavian
# Date: 04-04-2022
import argparse
import logging
import pathlib
from ctypes import alignment
import json


def pfamClanLink(clan):
    if clan != "No_clan":
        return f"href=\"http://pfam.xfam.org/clan/{clan}\" target=\"_blank\""
    return ''


def prependWhitespace(s, length):
    return ' '*(length-len(s))+s if len(s) < length else s


def createHmmRow(hmm, match):
    hmmRow = ""
    for h, m in zip(hmm, match):
        hClass = 'class=\"' + ('hmmmatch\"' if h ==
                               m else 'hmmplus\"' if m == '+' else '\"')
        hmmRow += f'<span {hClass}>{h}</span>'
    return(hmmRow)


def createSeqRow(seq, pp):
    seqRow = ""
    for s, p in zip(seq, pp):
        sClass = 'class=\"' + \
            ('ppstar\"' if p ==
             '*' else f'pp{p}"' if p in '0123456789' else '\"')
        seqRow += f'<span {sClass}>{s}</span>'
    return(seqRow)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count",
                        help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output",
                        help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True,
                        help="Path to AFPAP home")
    parser.add_argument("-j", "--json", help="Pfam json file")

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

    logging.info("Pfam HMM alignment...")
    with open(args.json, 'r') as json_file:
        pfam_data = json.load(json_file)
        pfam_table = ""
        if len(pfam_data) == 0:
            pfam_table += '<tr class="pfam-row2"><td colspan="14">No match found!</td></tr>'
        else:
            for recordID, record in enumerate(pfam_data):
                alignment_table = ""
                alignment_match = []
                alignment_section = 150
                matchLength = len(record["align"][0][11:])
                matchMaxPosition = max(
                    int(record["hmm"]["to"]), int(record["seq"]["to"]))
                for i in range(0, matchLength, alignment_section):
                    alignment_match.append([record["align"][0][11+i:11+i+alignment_section], record["align"][1][11+i:11+i +
                                                                                                                alignment_section], record["align"][2][11+i:11+i+alignment_section], record["align"][3][11+i:11+i+alignment_section]])
                for i, match in enumerate(alignment_match):
                    alignment_table += f'#HMM    {prependWhitespace(str(int(record["hmm"]["from"])+i*alignment_section),len(str(matchMaxPosition)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmMODEL">{createHmmRow(match[0],match[1])}</span><br>\n'
                    alignment_table += f'#MATCH  {prependWhitespace("",len(str(matchMaxPosition)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmMATCH">{match[1].replace(" ","&nbsp")}</span><br>\n'
                    alignment_table += f'#PP     {prependWhitespace("",len(str(matchMaxPosition)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmPP">{match[2]}</span><br>\n'
                    alignment_table += f'#SEQ    {prependWhitespace(str(int(record["seq"]["from"])+i*alignment_section),len(str(matchMaxPosition)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmSEQ">{createSeqRow(match[3],match[2])}</span><br><br>\n'

                pfam_table += f"""<tr class="pfam-row1">
                    <td>{record["type"]}</td>
                    <td>
                        <a href="http://pfam.xfam.org/family/{record["name"]}" target="_blank">{record["name"]}</a>
                    </td>
                    <td><a href="http://pfam.xfam.org/family/{record["acc"]}" target="_blank">{record["acc"]}</a></td>
                    <td>
                        <a {pfamClanLink(record["clan"])}>{record["clan"]}</a>
                    </td>
                    <td>{record["desc"]}</td>
                    <td>{record["env"]["from"]}</td>
                    <td>{record["env"]["to"]}</td>
                    <td>{record["seq"]["from"]}</td>
                    <td>{record["seq"]["to"]}</td>
                    <td>{record["hmm"]["from"]}</td>
                    <td>{record["hmm"]["to"]}</td>
                    <td>{record["model_length"]}</td>
                    <td>{record["bits"]}</td>
                    <td>{record["evalue"]}</td>
                </tr>
                <tr class="pfam-row2">
                    <td colspan="14">
                        <p>
                            <a class="btn btn-primary alignmentButton" data-toggle="collapse" href="#pfamAlignment{recordID}" role="button"
                                aria-expanded="false" aria-controls="pfamAlignment{recordID}">
                                Show Alignment
                            </a>
                        </p>
                        <div class="collapse" id="pfamAlignment{recordID}">
                            <div class="card card-body alignmentContainer">
                                {alignment_table}
                            </div>
                        </div>
                    </td>
                </tr>"""
        with open(f"{args.AFPAPpath}/config/pfam_template.html", 'r') as template:
            templateData = template.read()
            templateData = templateData.replace("---pfam---", pfam_table)
            with open(f'{args.outputDir}/work/multiqc_files/pfam_mqc.html', 'w') as f:
                print(templateData, file=f)


if __name__ == '__main__':
    main()
