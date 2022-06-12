
import pandas as pd
import os

# path_micro = {
#     {"Taxonid":562,"Taxoname":"Escherichia coli"},
#     {"Taxonid":573,"Taxoname":"Klebsiella pneumoniae"}
# }

card_gff = '/mnt/home/huanggy/idseq_wkf/card_3.2.2.gff3'

ARO = {}
with open(card_gff,'r') as fh:
    for line in fh :
        if line.startswith("seqid"):
            continue
        if not line:
            break
        line = line.strip()
        Taxoname = line.split("\t")[-1]
        Taxonid  = line.split("\t")[-2]
        aro = line.split("\t")[0].split(":")[-1]
        if aro not in ARO:
            ARO[aro] = {}
            ARO[aro]['Taxonid'] = Taxonid
            ARO[aro]['Taxoname'] = Taxoname
        else:
            pass

#print(ARO)
card_out = '/mnt/home/huanggy/biosoft/GenseqAMR/test/sample1_mNGS/sample1_rgi.txt'
tb = pd.read_table(card_out,header=0,sep = "\t")
out_file = os.path.join(os.path.dirname(card_out),"sample1_rgi_format.txt")
out = open(out_file,'w')
header=["Cut_Off","Pass_Bitscore","Best_Hit_Bitscore","Best_Hit_ARO","Best_Identities","Model_type","SNPs_in_Best_Hit_ARO","Other_SNPs","Resistance_Mechanism","AMR_Gene_Family","per_len_of_ref"]
print("Micro","\t".join(header),file=out,sep='\t')
for ind,row in tb.iterrows():
    aro_ = str(row["ARO"])
    Cut_Off = row["Cut_Off"]
    Pass_Bitscore = row["Pass_Bitscore"]
    Best_Hit_Bitscore = row["Best_Hit_Bitscore"]
    Best_Hit_ARO = row["Best_Hit_ARO"]
    Best_Identities = row["Best_Identities"]
    Model_type = row["Model_type"]
    SNPs_in_Best_Hit_ARO = row["SNPs_in_Best_Hit_ARO"]
    Other_SNPs = row["Other_SNPs"]
    Drug_Class = row["Drug Class"]
    Resistance_Mechanism = row["Resistance Mechanism"]
    AMR_Gene_Family = row["AMR Gene Family"]
    per_len_of_ref = row["Percentage Length of Reference Sequence"]

    if ARO.get(aro_):
        print(ARO[aro_]["Taxoname"],Cut_Off,Pass_Bitscore,Best_Hit_Bitscore,Best_Hit_ARO,Best_Identities,Model_type,
              SNPs_in_Best_Hit_ARO,Other_SNPs,Resistance_Mechanism,AMR_Gene_Family,per_len_of_ref,sep="\t",file=out)
    else:
        print("no")



# with  open(card_out,'r') as fh:
#     for line in fh:
#
#         if line.startswith("ORF_ID"):
#             continue
#         line = line.strip()
#         aro_ = line.split("\t")[10]
#         if aro_ in ARO:
#             ARO




