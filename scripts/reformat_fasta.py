#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
import numpy as np
import os,re,gzip
from Bio import SeqIO
import glob
from multiprocessing.dummy import Pool as ThreadPool


def parser_file2strain2taxid(assembly_summary_file):
#def parser_file2strain2taxid(assembly_summary_file,assembly_genome_num):
    tb = pd.read_table(assembly_summary_file,header = 0)
    assembly_genome_dict = {}
    for ind,row in tb.iterrows():
        speciess_taxid = row["taxid"]
        organism_name = row["organism_name"]
        file_name_ftp = row["ftp_path"]
        file_name = re.search(r'(GCF_.+)',file_name_ftp).groups()[0]
        if file_name not in assembly_genome_dict:
            assembly_genome_dict.update({file_name : {"assembly_genome_num":file_name,"organism_name":organism_name,"speciess_taxid":speciess_taxid}})
        else:
            pass
    # try :
    #     return assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
    # except:
    #     raise ValueError(assembly_genome_num ,"no exist")
    return  assembly_genome_dict

#assembly_genome_dict = parser_file2strain2taxid(assembly_summary_file=assembly_summary)

def parser_genomic_fna_file(file_path,quire_taxid):
    assembly_genome_info = os.path.dirname(file_path) + "/" + "assembly_genome_info"
    assembly_genome_info_out = open(assembly_genome_info,"a+")
    header = ["seq_id","taxid","description","genome length"]
    assembly_genome_info_out.write("\t".join(header) + "\n")
    with open(file_path,'r') as handle:
        record = SeqIO.parse(handle,"fasta")
        #print("------->",record)

        """
        print(record.__next__())   迭代对象
        ID: NC_002695.2
        Name: NC_002695.2
        Description: NC_002695.2 Escherichia coli O157:H7 str. Sakai DNA, complete genome
        Number of features: 0
        Seq('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAG...TTC')
        """
        n=1
        for seq_record in record:
            #SeqIO.write()
            """
            rec : <class 'Bio.SeqRecord.SeqRecord'>
            """
            # print("=====>",rec,"+++++\n+++++",type(rec))
            # print(type(rec.id))
            #rec.id = "taxid|123456 "
            #quire_taxid =  parser_file2strain2taxid(assembly_summary_file=assembly_summary,assembly_genome_num=assembly_genome_num)
            seq_record.id = "taxid" + "|" + str(quire_taxid)
            assembly_genome_info_out.write(seq_record.id + "\t" + str(quire_taxid) + "\t" + seq_record.description + "\t" + str(len(seq_record.seq)) + "\n")
            if hasattr(seq_record,"name"):
                delattr(seq_record,"name")

            #print(rec.__dict__['name'])
            #out_file = file_path + "_"  + str(n)
            #out_file = os.path.dirname(file_path) + "/split/" + os.path.basename(file_path)  + "_"  + str(n)
            #SeqIO.write(seq_record,out_file,"fasta")
            # print('---->',len(seq_record))
            # print('=====>',len(seq_record.seq))
            # print(seq_record.id)
            # print(len(seq_record.id))
            #print(seq_record.__dict__)
            n = n+1
        assembly_genome_info_out.close()

#单个文件解压 每条序列拆分,最后合并为一个fasta(https://www.jb51.net/article/156992.htm) ，构建index ，用于比对，
## 构建taxid 物种分类谱系
#parser_genomic_fna_file("/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/GCF_000008865.2_ASM886v2_genomic.fna")
#parser_genomic_fna_file("/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/GCF_900149625.2_5-27b_chromosome_genomic.fna")
if __name__ == "__main__":
    assembly_summary = "/mnt/data/NCBI_Refseq/Bacteria/representative_genome.dir.txt"
    file_list = glob.glob("/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/*_genomic.fna")

    # for i in file_list:
    #     assembly_genome_num = os.path.basename(i).strip("_genomic.fna")
    #     quire_taxid = parser_file2strain2taxid(assembly_summary_file=assembly_summary,assembly_genome_num=assembly_genome_num)
    #     print(i,"||||",assembly_genome_num,"||||||",quire_taxid)
    assembly_genome_dict = parser_file2strain2taxid(assembly_summary_file=assembly_summary)

    for file_path in file_list:
        assembly_genome_num = os.path.basename(file_path).replace("_genomic.fna", "")
        quire_id = assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
        if assembly_genome_dict[assembly_genome_num]:
            print(assembly_genome_dict[assembly_genome_num])
            quire_id = assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
        else:
            print(assembly_genome_num ,"不存在")
            quire_id = "NA"
        # print(file_path,"=========>",assembly_genome_num)
        parser_genomic_fna_file(file_path=file_path,quire_taxid=quire_id)





