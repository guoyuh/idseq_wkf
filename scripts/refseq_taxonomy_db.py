#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


from Bio import SeqIO
import  pandas as pd


##获取bac fungi nt 的seqid

def refseq_accession2taxid(file_path,out):
    with open(file_path, 'r') as handle,open(out,"w") as out:
        record = SeqIO.parse(handle, "fasta")
        for seq_record in record:
            #print(seq_record.id)
            out.write(seq_record.id + "\t" + seq_record.description + "\n")


def refseq_accession2taxid_2(file_path,out):
    with open(file_path, 'r') as handle,open(out,"w") as out:
        record = SeqIO.parse(handle, "fasta")
        out.write("taxid\taccession_id\tdescription\n")
        for seq_record in record:
            #print(seq_record.id)
            taxid = seq_record.id.split(":")[1].split("|")[1]
            accession_id = seq_record.id.split(":")[1].split("|")[2]
            out.write(taxid + "\t" + accession_id + "\t" + seq_record.description + "\n")




### 获取病毒的seqid 和 taxid
def viral_accession2taxid(file_path,out):
    with open(file_path, 'r') as handle,open(out,"w") as out:
        record = SeqIO.parse(handle, "fasta")
        out.write("taxid\taccession_id\tdescription\n")
        for seq_record in record:
            # print(seq_record.id)
            # print(seq_record.description)
            # print(dir(seq_record))

            taxid = seq_record.id.split(":")[1].split("|")[1]
            accession_id = seq_record.id.split(":")[1].split("|")[2]
            out.write(taxid + "\t" + accession_id + "\t" + seq_record.description + "\n")


#refseq_accession2taxid(file_path="/mnt/data/kraken2_db/library/bacteria/library.fna",out="/mnt/data/NCBI_Refseq/bacteria/refseq2taxid.txt")
#refseq_accession2taxid(file_path="/mnt/data/kraken2_db/library/fungi/library.fna",out="/mnt/data/NCBI_Refseq/fungi/refseq2taxid.txt")
#viral_accession2taxid(file_path="/mnt/data/kraken2_db/library/viral/library.fna",out="/mnt/data/kraken2_db/library/viral/refseq2taxid.txt")
#refseq_accession2taxid_2(file_path="/mnt/data/kraken2_mask/library/protozoa/library.fna",out="/mnt/data/NCBI_Refseq/protozoa/refseq2taxid.txt")
# tb1 = pd.read_table("/mnt/data/NCBI_Refseq/bacteria/refseq2taxid.txt",header=0,sep = "\t")
# tb2 = pd.read_table("/mnt/data/NCBI_Refseq/fungi/refseq2taxid.txt",header=0,sep = "\t")
# # print(tb1.head())
# # print(tb2.head())
# tb = tb1.append(tb2,ignore_index=True)
# # print(tb.head(),tb.info)
#
# tb3 = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq-release207.catalog",header=None,sep="\t")
# print(tb3.head())
# tb3.columns = ["taxonomy ID","species name","accession.version","refseq description","refseq status","length"]
#
# out=tb3[tb3["accession.version"].isin(tb["seqid"])]
# out.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun.txt",sep="\t",index=False)


# with open("/mnt/data/kraken2_db/library/viral/library.fna", 'r') as handle :
#     record = SeqIO.parse(handle, "fasta")
#     for seq_record in record:
#         print(seq_record.id)
#         print(seq_record.description)
#         print(dir(seq_record))
#
#         taxid = seq_record.id.split(":")[1].split("|")[1]
#         accession_id = seq_record.id.split(":")[1].split("|")[2]
#         break


#
# tb = pd.read_table("/mnt/data/kraken2_db/library/viral/refseq2taxid.txt",header=0,sep="\t")
# tb3 = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq-release207.catalog",header=None,sep="\t")
# tb3.columns = ["taxonomy ID","species name","accession.version","refseq description","refseq status","length"]
# out=tb3[tb3["accession.version"].isin(tb["accession_id"])]
# out.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_viral.txt",sep="\t",index=False)


#思路：
    # 解析protozoa/library.fna ，获取 taxid	accession_id
    # 从RefSeq-release207.catalog   根据accession_id    获取其余信息  "species name","accession.version","refseq description","refseq status","length"
    # 跟据taxid 从 taxonkit    获取格式化的lineage

#20211014
tb = pd.read_table("/mnt/data/NCBI_Refseq/protozoa/refseq2taxid.txt",header=0,sep="\t")
tb3 = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq-release207.catalog",header=None,sep="\t")
tb3.columns = ["taxonomy ID","species name","accession.version","refseq description","refseq status","length"]
out=tb3[tb3["accession.version"].isin(tb["accession_id"])]
out.to_csv("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_protozoa.txt",sep="\t",index=False)




