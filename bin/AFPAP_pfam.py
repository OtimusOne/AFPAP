'''
# File name: AFPAP_pfam.py
# Description: Python script for generating Pfam family alignment
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import json
import logging
import pathlib


def prepend_whitespace(string="", length=1):
    '''
    Prepend whitespace to string.

    Parameters:
        string(str)
        length(int)
    Returns:
        string with added whitespace
    '''
    return ' '*(length-len(string))+string if len(string) < length else string


def create_hmm_row(hmm_string="", match_string=""):
    '''
    Create HTML row for HMM alignment.

    Parameters:
        hmm_string: HMM model sequence.
        match_string: HMM model match.
    Returns:
        HTML HMM alignment row.
    '''
    hmm_row = ""
    for hmm, match in zip(hmm_string, match_string):
        h_class = 'class=\"' + ('hmmmatch\"' if hmm == match else 'hmmplus\"' if match == '+' else '\"')
        hmm_row += f'<span {h_class}>{hmm}</span>'
    return hmm_row


def create_seq_row(seq_string="", pp_string=""):
    '''
    Create HTML row for sequence alignment.

    Parameters:
        seq_string: Sequence alignment.
        pp_string: Sequence posterior probability.
    Returns:
        HTML sequence alignment row.
    '''
    seq_row = ""
    for seq, p_probability in zip(seq_string, pp_string):
        seq_class = 'class=\"' + ('ppstar\"' if p_probability == '*' else f'pp{p_probability}"' if p_probability in '0123456789' else '\"')
        seq_row += f'<span {seq_class}>{seq}</span>'
    return seq_row


def main():
    '''
    Python script for generating Pfam family alignment.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True,  help="Path to AFPAP home")
    parser.add_argument("-j", "--json", required=True, help="Pfam json file")

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

    logging.info("Pfam HMM alignment...")
    with open(args.json, 'r', encoding="utf8") as json_file:
        pfam_data = json.load(json_file)
        pfam_table = ""
        if len(pfam_data) == 0:
            pfam_table += '<tr class="pfam-row2"><td colspan="14">No match found!</td></tr>'
        else:
            for record_id, record in enumerate(pfam_data):
                alignment_table = ""
                alignment_match = []
                alignment_section = 150
                match_length = len(record["align"][0][11:])
                match_max_position = max(int(record["hmm"]["to"]), int(record["seq"]["to"]))
                for i in range(0, match_length, alignment_section):
                    alignment_match.append([record["align"][0][11+i:11+i+alignment_section], record["align"][1][11+i:11+i + alignment_section],
                                           record["align"][2][11+i:11+i+alignment_section], record["align"][3][11+i:11+i+alignment_section]])
                for i, match in enumerate(alignment_match):
                    alignment_table += f'#HMM    {prepend_whitespace(str(int(record["hmm"]["from"])+i*alignment_section),len(str(match_max_position)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmMODEL">{create_hmm_row(match[0],match[1])}</span><br>\n'
                    alignment_table += f'#MATCH  {prepend_whitespace("",len(str(match_max_position)))} '.replace(' ', "&nbsp")+f'<span class="hmmMATCH">{match[1].replace(" ","&nbsp")}</span><br>\n'
                    alignment_table += f'#PP     {prepend_whitespace("",len(str(match_max_position)))} '.replace(' ', "&nbsp")+f'<span class="hmmPP">{match[2]}</span><br>\n'
                    alignment_table += f'#SEQ    {prepend_whitespace(str(int(record["seq"]["from"])+i*alignment_section),len(str(match_max_position)))} '.replace(
                        ' ', "&nbsp")+f'<span class="hmmSEQ">{create_seq_row(match[3],match[2])}</span><br><br>\n'

                entry_clan = f"href=\"http://pfam.xfam.org/clan/{record['clan']}\" target=\"_blank\"" if record['clan'] != "No_clan" else ''
                pfam_table += f"""<tr class="pfam-row1">
                    <td>{record["type"]}</td>
                    <td>
                        <a href="http://pfam.xfam.org/family/{record["name"]}" target="_blank">{record["name"]}</a>
                    </td>
                    <td><a href="http://pfam.xfam.org/family/{record["acc"]}" target="_blank">{record["acc"]}</a></td>
                    <td>
                        <a {entry_clan}>{record["clan"]}</a>
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
                            <a class="btn btn-primary alignmentButton" data-toggle="collapse" href="#pfamAlignment{record_id}" role="button"
                                aria-expanded="false" aria-controls="pfamAlignment{record_id}">
                                Show Alignment
                            </a>
                        </p>
                        <div class="collapse" id="pfamAlignment{record_id}">
                            <div class="card card-body alignmentContainer">
                                {alignment_table}
                            </div>
                        </div>
                    </td>
                </tr>"""
        with open(f"{args.AFPAPpath}/config/pfam_template.html", 'r', encoding="utf8") as template:
            template_data = template.read()
            template_data = template_data.replace("---pfam---", pfam_table)
            with open(f'{args.outputDir}/work/multiqc_files/pfam_mqc.html', 'w', encoding="utf8") as mqc_file:
                print(template_data, file=mqc_file)


if __name__ == '__main__':
    main()
