#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd

input = "/mnt/data/NCBI/taxonomy/refseq_taxonomy/refseq_protozoa_taxid_lineage_db"
tb = pd.read_table(input,sep = '\t',header = 0,)
print(tb)

tb2 = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_protozoa.txt",sep='\t',header=0)
tb2.columns = ["taxid","species_name","accession_version" ,"refseq_description","refseq_status" ,"length"]
print('====>',tb2.head())
print('====>',tb2.shape)
mm=pd.merge(tb,tb2,how='right',on="taxid",left_index = False,right_index=False)
mm=mm[["taxid","species_name","accession_version","refseq_description","refseq_status","length","names_lineage","taxid_lineage","name","format_names_lineage"]]
mm.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_protozoa.lineage_db",index=False,sep="\t")

RefSeq_bac_fun_viral_db = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db",sep='\t',header=0)

new_db =RefSeq_bac_fun_viral_db.append(mm,ignore_index=True)
new_db.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral_protozoa.lineage_db",index=False,sep="\t")