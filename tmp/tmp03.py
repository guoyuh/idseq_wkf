

import pandas as pd 


# tb1 = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.txt", sep="\t",header=0)
# tb1.columns = ["taxid","species_name","accession_version" ,"refseq_description","refseq_status" ,"length"]
# print(tb1.head())

# tb2= pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/refseq_taxid_lineage_db",sep = "\t",header = 0)
# print(tb2.head())
# tb = pd.merge(tb1,tb2,how="left" ,on = ["taxid"])


# left = pd.DataFrame({'key': ['K0', 'K1', 'K2', 'K3'],
#                        'A': ['A0', 'A1', 'A2', 'A3'],
#                        'B': ['B0', 'B1', 'B2', 'B3']})
# right = pd.DataFrame({'key': ['K0', 'K1', 'K2', 'K3'],
#                         'C': ['C0', 'C1', 'C2', 'C3'],
#                         'D': ['D0', 'D1', 'D2', 'D3']})
# result = pd.merge(left, right, on='key')

# # on参数传递的key作为连接键
# print(left)
# print(right)
# print(result)

#tb.head().to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db.csv",sep=",",index=False)
#tb.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db",sep="\t",index=False)

# taxid_lineage_db = "/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db"
# errfile = open("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db.errline","w")
# cnt = 0 
# with open (taxid_lineage_db,'r') as  db_path_f:
#     #next(db_path_f)
#     for line in db_path_f:
#         line = line.strip()
#         if line.startswith("taxid"):
#             pass
#         else:
#             cnt +=1
#             if len(line.strip().split("\t")) < 10:
#                 print(line,sep="\t",file=errfile)

# errfile.close()


import pandas as pd
tb =pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db",sep="\t",header = 0)
print(tb)

tb.fillna("abc")
print(tb)
tb.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db.tmp",sep="\t",header = 0)
