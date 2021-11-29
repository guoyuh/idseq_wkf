#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
from Bio import SeqIO

# data =pd.read_csv("/mnt/home/huanggy/project/meta/result/align/top10.id.test",sep="\t",header=0)
# print(data)
# for k ,v in data.items():
#     print(v)

file_path = "/mnt/home/huanggy/project/meta/temp/qc/norg_13.R1_kneaddata_unmatched.priceseq.merge.fa"

# for ind ,row in  data.iterrows():
#     qseqid = row["qseqid"].lstrip("[").rstrip("]").split(",")
#     #qseqid = row["qseqid"]
#     print(type(qseqid),qseqid)
#     taxid = row["taxid"]
#     for id in qseqid:
#         id = id.replace("/2","")
#         print(id)
#         break

# with open(file_path, 'r') as handle :
#     record = SeqIO.parse(handle, "fasta")
#     for seq_record in record:
#         print(seq_record.id)


from pathlib import Path
import json

outdir = "/mnt/home/huanggy/project/meta/result/align/"
outdir = Path(outdir)

with open ("/mnt/home/huanggy/project/meta/result/align/top10.json") as f:

    a = json.load(f)
    print(a)
    print(type(a))

    for k , v in a.items():
        #print(k)
        if k == "qseqid":
            print(v)
            for specise ,seqid in v.items():
                for s in seqid:
                    s = s.replace("/2","")
                    s = s.replace("/1","")
                    save_file = outdir/(specise.replace(" ","_") + ".fa")
                    #print(save_file)
                    with open(file_path, 'r') as handle ,open(save_file,"a") as out :
                        record = SeqIO.parse(handle, "fasta")
                        for seq_record in record:
                            #print(seq_record.id)
                            if s == seq_record.id:
                                seq = str(seq_record.seq)
                                header = ">" + s
                                out.write(header + "\n" +seq + "\n" )

                    


