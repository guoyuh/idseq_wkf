import sys
import os
import json
import argparse
import csv

"""
This script is used to convert CARD json file to GFF3 format.
The card.json file can be downloaded at https://card.mcmaster.ca/download under CARD data.
Usage:
	python convert_card_json_to_gff3.py -i card.json
"""


def main(args):
    if args.input_file == None:
        exit("Missing input card json file")

    with open(os.path.join(args.input_file), 'r') as jfile:
        data = json.load(jfile)
        try:
            version = data["_version"]
            print(version)
        except Exception as e:
            print("Error: missing version number")
            exit()
        with open(os.path.join("card_{}.gff3".format(version)), "w") as af:
            sequences = []
            headers = []
            body = []
            writer = csv.writer(af, delimiter='\t')
            headers.append(['##gff-version 3.2.1'])
            af.write("\t".join(["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "taxonomy_id","taxonomy_name"]) + "\n")
            for i in data:
                if i.isdigit():
                    if "model_sequences" in data[i].keys():
                        for k in data[i]["model_sequences"]["sequence"]:
                            _source = "{}".format("CARD")
                            _type = "{}".format("CDS")
                            _phase = "{}".format(".")
                            _score = "{}".format(".")
                            _seqid = "{gi}{aro}".format(
                                gi="gi|" + data[i]["model_sequences"]["sequence"][k]["dna_sequence"]["accession"],
                                aro="|ARO:" + data[i]["ARO_accession"])
                            _start = "{}".format(
                                int(data[i]["model_sequences"]["sequence"][k]["dna_sequence"]["fmin"]) + 1)
                            _end = "{}".format(data[i]["model_sequences"]["sequence"][k]["dna_sequence"]["fmax"])
                            _strand = "{}".format(data[i]["model_sequences"]["sequence"][k]["dna_sequence"]["strand"])
                            _attributes = "Name={name};Alias={aro}".format(name=data[i]["ARO_name"],aro="ARO:" + data[i]["ARO_accession"])
                            taxonomy_id = "{}".format(data[i]['model_sequences']['sequence'][k]['NCBI_taxonomy']['NCBI_taxonomy_id'])

                            taxonomy_name = "{}".format(data[i]['model_sequences']['sequence'][k]['NCBI_taxonomy']['NCBI_taxonomy_name'])

                            af.write("\t".join([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes, taxonomy_id,taxonomy_name]) + '\n')
                            headers.append(['##sequence-region ' + _seqid + ' ' + _start + ' ' + _end])
                            #body.append(c)
                            sequences.append(format_fasta(str(_seqid), "{}".format(data[i]["model_sequences"]["sequence"][k]["dna_sequence"]["sequence"])))
        # headers
        # for head_item in headers:
        #     af.write(head_item[0] + '\n')
        # body
        # for body_item in body:
        #     writer.writerow(body_item)


        # # footer
        # writer.writerow(["##FASTA"])
        # for sequence in sequences:
        # 	af.write(sequence)

    with open(os.path.join("card_{}.fa".format(version)),'w') as fh:
        for sequence in sequences:
            fh.write(sequence)

def format_fasta(name, sequence):
    fasta_string = '>' + name + '\n' + sequence + '\n'
    return fasta_string


def run():
    parser = argparse.ArgumentParser(description='convert card json to gff3')
    parser.add_argument('-i', '--input_file', dest="input_file", default=None, required=True,
                        help='card.json input file')
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()