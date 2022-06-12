import json,os

#card_json = "/mnt/home/huanggy/CARD/CARD-DATA/card.json"
card_json = "/mnt/data3/test/sample1_rgi.json"
with open(card_json,'r') as jsfile:
    data = json.load(jsfile)

    #print(data)
    # try:
    #     version = data["_version"]
    #     print(version)
    # except Exception as e:
    #     print("Error: missing version number")
    #     exit()

    print(type(data))
    ini = 1
    b = 10 
    for k,v in data.items():
        print(k,"==========>",v)
        ini +=1
        if ini > 3:
            break

"""
annotations = []
"""
#write card reference fasta (FASTA format)
"""
working_directory = "/mnt/home/huanggy/CARD/CARD-DATA"
with open(os.path.join(working_directory, "card_database_v{}.fasta".format(version)), 'w') as fout:

    for i in data:
        if i.isdigit():
            drug_class = []
            mechanism = []
            group = []
            if "ARO_category" in data[i]:
                for c in data[i]["ARO_category"]:
                    if "category_aro_class_name" in data[i]["ARO_category"][c]: 
                        if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
                            drug_class.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
                        if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
                            mechanism.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
                        if data[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
                            group.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))

            try:
                for seq in data[i]["model_sequences"]["sequence"]:
                    # header used to be able to validate CARD sequences with genbank sequences
                    header = ("gb|{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}|Taxonid:{taxonomy_id}|Taxoname:{taxonomy_name}".format(
                        ARO_accession=data[i]['ARO_accession'],
                        model_id=data[i]['model_id'],
                        ARO_name=(data[i]['ARO_name']).replace(" ", "_"),
                        ncbi=data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"],
                        taxonomy_id=data[i]['model_sequences']['sequence'][seq]["NCBI_taxonomy"]["NCBI_taxonomy_id"],
                        taxonomy_name=data[i]['model_sequences']['sequence'][seq]["NCBI_taxonomy"]["NCBI_taxonomy_name"]
                        ))
                    sequence = data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["sequence"]
                    fout.write(">{}\n".format(header))
                    fout.write("{}\n".format(sequence))
            except Exception as e:
                print("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))


        # print("drug_class:",drug_class)
        # print("mechanism:",mechanism)
        # print("group:",group)

"""



