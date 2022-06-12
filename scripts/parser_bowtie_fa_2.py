

from Bio import SeqIO
import pandas as pd 
import os
import datetime


def parser_accessionNum2species_from_genomic_fna(fa_file):
    with open(fa_file,'r') as handle:
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
        for (ind,seq_record) in enumerate(record):
            #print(dir(seq_record))
            #print(ind,seq_record.id,seq_record.name,seq_record.description,len(seq_record))
            # print(seq_record.id)
            # print(seq_record.name)
            # print(seq_record.description)
            # print(len(seq_record))
            # print(len(seq_record.seq))
            # print('========>',seq_record.letter_annotations)
            
            yield seq_record


def parser_accessionNum2species_from_genomic_fna_2(fa_file):
    scaffold_length = {}
    with open(fa_file,'r') as handle,open(os.path.join(os.path.dirname(fa_file),'bacteria_fungi_protozoa_viral.scaffold.info'),'w') as f_out:
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
        print ("\t".join(["seq_id","organism","g_len"]),file=f_out)
        for (ind,seq_record) in enumerate(record):
            #print(dir(seq_record))
            #print(ind,seq_record.id,seq_record.name,seq_record.description,len(seq_record))
            # print(seq_record.id)
            # print(seq_record.name)
            # print(seq_record.description)
            # print(len(seq_record))
            # print(len(seq_record.seq))
            # print('========>',seq_record.letter_annotations)
            
            
            seq_id = seq_record.id
            seq_organism = seq_record.description.strip(seq_id).lstrip(" ")
            s_len = len(seq_record.seq)
            if seq_id not in scaffold_length:
                scaffold_length[seq_id]={"seq_organism":seq_organism,"seq_len":s_len}
                print(seq_id,seq_organism,s_len,sep="\t",file=f_out)
            else:
                pass
            




def get_species_meanDepth_breadth(bowtie2_fa, depth):
    start = datetime.datetime.now()
    print (start)
    tb = pd.read_table(depth, header=None, sep="\t", index_col=None)

    # print(tb.head)
    tb.columns = ["scaffold", "pos", "depth"]
    # print(tb.head)
    seq_record_generator = parser_accessionNum2species_from_genomic_fna(bowtie2_fa)

    for seq_record in seq_record_generator:
        seq_id = seq_record.id
        seq_organism = seq_record.description.strip(seq_id).lstrip(" ")
        s_len = len(seq_record.seq)
        #print(seq_id,seq_organism,s_len)


        #if tb["scaffold"].str.contains(seq_id):
        if seq_id in tb["scaffold"].to_list():
            scaffold_tb = tb[tb["scaffold"] == seq_id]
            
            pos_len = len(pd.Series(scaffold_tb["pos"].values))
            breadth  = "{:.2%}".format(pos_len /s_len )
            meandepth = "{:.2f}".format(pd.Series(scaffold_tb["depth"].values).values.sum()/pos_len)
            print(seq_id,seq_organism,s_len,meandepth,breadth,sep = "  ")

    end = datetime.datetime.now()
    print("total time is:%s"%(end-start))
        

def get_scaffold_length(bacteria_fungi_protozoa_viral_scaffold_info):
    scaffold_length = {}
    with open(bacteria_fungi_protozoa_viral_scaffold_info) as fh:
        for line in fh:
            if line.startswith("seq_id"):
                continue
            if not line:
                break
            line = line.strip()
            seq_id = line.split("\t")[0]
            s_len = line.split("\t")[2]
            if seq_id not in scaffold_length:
                scaffold_length[seq_id]={"seq_len":s_len}
            else:
                pass
    
    return scaffold_length



def get_s_meanDepth_breadth(scaffold_info,depth):
    start = datetime.datetime.now()
    scaffold_length = get_scaffold_length(scaffold_info)
    print(len(scaffold_length))
    tb = pd.read_table(depth, header=None, sep="\t", index_col=None)
    # print(tb.head)
    tb.columns = ["scaffold", "pos", "depth"]


    for scaffold in list(set(tb["scaffold"])):
        scaffold_tb = tb[tb["scaffold"] == scaffold]
        
        pos_len = len(pd.Series(scaffold_tb["pos"].values))
        meandepth = "{:.2f}".format(pd.Series(scaffold_tb["depth"].values).values.sum()/pos_len)
        breadth  = "{:.2%}".format(pos_len / int(scaffold_length[scaffold]["seq_len"]))

        print(scaffold,meandepth,breadth)
        del scaffold_tb

    end = datetime.datetime.now()
    print("total time is:%s"%(end-start))


fa_file="/mnt/data/kraken2_db/bacteria_fungi_protozoa_viral.mask.fna"
#fa_file="/mnt/data/kraken2_mask/library/GCF_000002855.3_ASM285v2_genomic.fna"
#fa_file="/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_001618505.2_ASM161850v2_genomic.fna"
#fa_file="/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_002277995.2_ASM227799v2_genomic.fna"
#parser_accessionNum2species_from_genomic_fna_2(fa_file)
# scaffold_length = parser_accessionNum2species_from_genomic_fna(fa_file)
scaffold_info = "/mnt/data/kraken2_db/bacteria_fungi_protozoa_viral.scaffold.info"
depth_out = "/mnt/home/huanggy/project/20220408/result/21JS944006/align/21JS944006.depth"

get_species_meanDepth_breadth(bowtie2_fa=fa_file,depth= depth_out)
#get_s_meanDepth_breadth(scaffold_info=scaffold_info,depth=depth_out)
"""
50bp /5M  =  0.001 
breadth_cutoff = 0.001

def plot(scaffold_tb):
    pass

if breadth_cutoff > 0.01 :
    plot(scaffold_tb)


#根据scaffold length 等分20 份 ，像上取整 得到20 份 区间右端点与区间的深度（平均深度？？？？？？？？？？？）字典
#然后画柱子图？？？？
"""
