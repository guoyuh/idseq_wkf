#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import re
# out  = open("add.txt","w")
# with open ("/mnt/data/NCBI_Refseq/bacteria/PubMLST.bacteria.delspp.txt",'r') as fh:
#     with open("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt","r") as fh2:
#
#         for line1 in fh:
#             line1 = line1.strip()
#             for line2 in fh2:
#                 line2 = line2.strip()
#                 #print("=======>",line2)
#                 # id = str(line1)
#                 if re.search(line1,line2):
#                     print("===>",line2)
# out.close()

assembl_file = '/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt'
with open(assembl_file,'r') as f:
    with open("/mnt/data/NCBI_Refseq/bacteria/PubMLST.bacteria.delspp.txt",'r') as fh:
        for line in f :
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                refseq_category = line.strip().split("\t")[4]
                taxid = line.strip().split('\t')[5]
                species_taxid = line.strip().split('\t')[6]
                organism_name = line.strip().split("\t")[7]
                ftp_path = line.strip().split('\t')[19]

                for line2 in fh:
                    line2 = line2.strip()
                    if re.search(line2,organism_name):
                         print('=====>',line)
