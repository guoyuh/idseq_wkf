
import pandas as pd


tb1 = pd.read_table("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt",sep="\t",header=0)

my_tb = pd.read_table("/mnt/data/NCBI_Refseq/my_bac_assembly_summary_refseq.info20211101",sep="\t",header=0)
print(tb1.head())

print(tb1.columns)
print(my_tb.head())
index =tb1["assembly_accession"].isin(my_tb["assembly_accession"])

add_tb = tb1["GCF_903181465.1"].is

["US Food and Drug Administration","FDA","NIH"]



def my_fun(s):
    if s.find("US Food and Drug Administration") != -1:
        return True
    elif s.find("FDA") != -1:
        return True
    
    else:
        return False



#print(tb1['submitter'].str.contains("US Food and Drug Administration")['submitter'].str.contains("FDA")['submitter'].str.contains("NIH"))

