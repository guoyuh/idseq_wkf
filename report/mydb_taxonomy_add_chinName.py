#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


name_dict = { }

with open('/mnt/data/kraken2_mask/taxonomy/种_拉丁2中文名.txt','r') as fh:
    for line in fh:
        array = line.strip().split("\t")
        if array[0] not in name_dict:
            name_dict.update({array[0]:array[1]})
        else:
            pass


out_file = "/mnt/data/kraken2_mask/taxonomy/mydb_taxonomy_with_chinName"
out = open(out_file,'w')
header = ['taxid','super_taxid','level','distance','Latin_name','chinName']
print("\t".join(header),file=out)
with open("/mnt/data/kraken2_mask/taxonomy/mydb_taxonomy_2.txt","r") as fh:
    for line in fh:
        if line.startswith("taxid"):
            continue

        line = line.strip()
        taxid,super_taxid,level,distance,Latin_name = line.split("\t")
        if level == "S" and Latin_name in name_dict:
            print(line,name_dict[Latin_name],sep="\t",file=out)
        else:
            print(line, "-", sep="\t", file=out)

out.close()




