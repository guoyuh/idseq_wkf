

from Bio import SeqIO
import pandas as pd 




def parser_accessionNum2species_from_genomic_fna(fa_file):
    scaffold_length = {}
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
            
            seq_id = seq_record.id
            seq_organism = seq_record.description.strip(seq_id).lstrip(" ")
            s_len = len(seq_record.seq)

            if seq_id not in scaffold_length:
                scaffold_length[seq_id]={"seq_organism":seq_organism,"seq_len":s_len}
            else:
                pass
    return scaffold_length


fa_file="/mnt/data/kraken2_db/bacteria_fungi_protozoa_viral.mask.fna"
#fa_file="/mnt/data/kraken2_mask/library/GCF_000002855.3_ASM285v2_genomic.fna"
#fa_file="/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_001618505.2_ASM161850v2_genomic.fna"
#fa_file="/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_002277995.2_ASM227799v2_genomic.fna"
#scaffold_length = parser_accessionNum2species_from_genomic_fna(fa_file)
#print(scaffold_length)
scaffold_length ={
    "NZ_WBWE01000003.1":100,
    "NZ_JMLH01000056.1":100
}
depth_out = "/mnt/home/huanggy/project/20220408/result/21JS944006/align/21JS944006.depth"

tb = pd.read_table(depth_out,header=None,sep="\t",index_col=None)

# print(tb.head)
tb.columns = ["scaffold","pos","depth"]
# print(tb.head)


scaffold_list = list(set(tb["scaffold"]))  ##大列表 耗内存
print(scaffold_list)
print(len(scaffold_list))

#print(tb[tb["scaffold"] == "NZ_VMQV01000054.1"])
sub_tb = tb[tb["scaffold"] == "NZ_VMQV01000054.1"]
print(sub_tb["scaffold"].to_list())

# for scaffold in scaffold_list:
#     scaffold_tb = tb[tb["scaffold"] == scaffold]

#     meandepth = pd.Series(scaffold_tb["depth"].values).values.mean()
#     pos_len = len(pd.Series(scaffold_tb["pos"].values))

#     if scaffold in scaffold_length:
#         breadth  = "{:.2%}".format(pos_len / scaffold_length[scaffold]["seq_len"])
#     else :
#         breadth = 0

#     print(scaffold,scaffold_length[scaffold]["seq_organism"],meandepth,breadth,sep = "  ")
#根据scaffold_tb  画深度覆盖度图
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