#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


# file1 = "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/GCF_014621655.1_ASM1462165v1_genomic.fna.gz"
# outfile = "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/nt.clin.fasta"
#
# with open(outfile,"w") as fh:
#     pass


import pandas as pd
import numpy as np
import os,re,gzip
from Bio import SeqIO
import glob
from multiprocessing import Pool
from subprocess import Popen

#assembly_summary = "/mnt/data/NCBI_Refseq/Bacteria/representative_genome.dir.txt"
#def parser_file2strain2taxid(assembly_summary_file):
def parser_file2strain2taxid(assembly_summary_file,assembly_genome_num):
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
    try :
        return assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
    except:
        raise ValueError(assembly_genome_num ,"no exist")

#assembly_genome_dict = parser_file2strain2taxid(assembly_summary_file=assembly_summary)

def parser_genomic_fna_file_max_id(file_path):
    assembly_genome_num = os.path.basename(file_path).replace("_genomic.fna", "")
    assembly_genome_info = os.path.dirname(file_path) + "/" + "assembly_genome_info"
    #assembly_genome_info_out = open(assembly_genome_info,"w")
    header = ["seq_id","taxid","description","genome length"]
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

        #rec =[ max(len(rec.seq)) for rec in record ]
        # init_seq_length = 0
        # for seq_record in record:
        #     #current_seq_length =len(seq_record.seq)
        #     init_seq_length = len(seq_record.seq)
        #     continue
        record_list = []
        for (ind,seq_record) in enumerate(record):
            print(ind,len(seq_record))
            record_list.append((ind,len(seq_record)))

        #print(max([ i[1] for i in record_list]),[ i[1] for i in record_list].index(max([ i[1] for i in record_list])))
        max_length = max([ i[1] for i in record_list])
        max_length_ind = [ i[1] for i in record_list].index(max([ i[1] for i in record_list]))
        return max_length_ind

def change_genomic_fna_file(file_path):
    assembly_genome_num = os.path.basename(file_path).replace("_genomic.fna", "")

    max_length_ind = parser_genomic_fna_file_max_id(file_path=file_path)
    with open(file_path,'r') as handle:
        record = SeqIO.parse(handle,"fasta")
        for (ind,seq_record) in enumerate(record):
            # print("--------->")
            # print(ind,len(seq_record))
            if ind == max_length_ind:
                quire_taxid =  parser_file2strain2taxid(assembly_summary_file=assembly_summary,assembly_genome_num=assembly_genome_num)
                seq_record.id = "kraken:taxid" + "|" + str(quire_taxid)
                if hasattr(seq_record,"name"):
                    delattr(seq_record,"name")
                out_file = os.path.dirname(file_path) + "/kraken2/" + os.path.basename(file_path)
                SeqIO.write(seq_record,out_file,"fasta")
                # print(seq_record.id)
                # print(seq_record.seq)
                return seq_record


#change_genomic_fna_file(file_path="/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/GCF_000008865.2_ASM886v2_genomic.fna")

def get_assembly_info(assembly_summary_file):
    tb = pd.read_table(assembly_summary_file,header = 0)
    assembly_genome_dict = {}
    for ind,row in tb.iterrows():
        speciess_taxid = row["taxid"]
        organism_name = row["organism_name"]
        assembly_accession = row["# assembly_accession"]
        asm_name = row["asm_name"]
        file_name_ftp = row["ftp_path"]
        file_name = re.search(r'(GCF_.+)',file_name_ftp).groups()[0]
        if file_name not in assembly_genome_dict.keys():
            print(file_name,file_name_ftp)
            assembly_genome_dict.update({file_name : {"assembly_genome_num":file_name,"organism_name":organism_name,"speciess_taxid":speciess_taxid,"assembly_accession":assembly_accession,"asm_name":asm_name}})
        else:
            pass

        return assembly_genome_dict



if __name__ == "__main__":
    assembly_summary = "/mnt/data/NCBI_Refseq/Bacteria/representative_genome.dir.txt"
    file_list = glob.glob("/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/*_genomic.fna")

    # for i in file_list:
    #     assembly_genome_num = os.path.basename(i).strip("_genomic.fna")
    #     quire_taxid = parser_file2strain2taxid(assembly_summary_file=assembly_summary,assembly_genome_num=assembly_genome_num)
    #     print(i,"||||",assembly_genome_num,"||||||",quire_taxid)
    #assembly_genome_dict = parser_file2strain2taxid(assembly_summary_file=assembly_summary)

    # for file_path in file_list:
    #     assembly_genome_num = os.path.basename(file_path).replace("_genomic.fna", "")
    #     quire_id = assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
    #     if assembly_genome_dict[assembly_genome_num]:
    #         print(assembly_genome_dict[assembly_genome_num])
    #         quire_id = assembly_genome_dict[assembly_genome_num]["speciess_taxid"]
    #     else:
    #         print(assembly_genome_num ,"不存在")
    #         quire_id = "NA"
    #     # print(file_path,"=========>",assembly_genome_num)
    #     parser_genomic_fna_file(file_path=file_path,quire_taxid=quire_id)

    # task_list = [f for f in file_list ]
    # print ("Total task: {task_num}\n".format(task_num=len(task_list)))
    # pool = Pool(10)
    # pool.map(change_genomic_fna_file, task_list)
    # pool.close()
    # pool.join()



    # for f in file_list:
    #     assembly_genome_num = os.path.basename(f).replace("_genomic.fna", "")
    #     assembly_genome_info = os.path.dirname(f) + "/" + "assembly_genome_info_20210824"
    #     assembly_genome_info_out = open(assembly_genome_info,"w")
    #     header = ["seq_id","taxid","description","genome length"]
    #     assembly_genome_info_out.write("\t".join(header) + "\n")
    #     seq_record = change_genomic_fna_file(f)
    #     quire_taxid = parser_file2strain2taxid(assembly_summary_file=assembly_summary,assembly_genome_num=assembly_genome_num)
    #     seq_record.id = "taxid" + "|" + str(quire_taxid)
    #     assembly_genome_info_out.write(seq_record.id + "\t" + str(quire_taxid) + "\t" + seq_record.description + "\t" + str(len(seq_record.seq)) + "\n")
    #     assembly_genome_info_out.close()

    assembly_genome_dict = get_assembly_info(assembly_summary_file=assembly_summary)
    print(assembly_genome_dict)
