#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang





ref_id_file = "/mnt/data/NCBI_Refseq/viral/my_virus_species.id.txt"
ref_id = {}
with open( ref_id_file,'r') as handle:
    for line in handle:
         qid = line.strip()
         if qid not in ref_id:
             ref_id[qid] = 1
#print(ref_id)



assembl_file = '/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt'
out = open("/mnt/data/NCBI_Refseq/viral/my_viral_assembly_summary_refseq.20210830","w")
out2 = open("/mnt/data/NCBI_Refseq/viral/my_viral_assembly_summary_refseq.20210830_not_reference","w")
with open(assembl_file,'r') as f:
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

            if taxid in ref_id or species_taxid in ref_id:
                if refseq_category == "representative genome" or refseq_category == "reference genome":

                    #print(organism_name, taxid, species_taxid,refseq_category,ftp_path,sep="\t",file=out)
                    print(line, sep="\t", file=out)
                else:
                    #print(organism_name, taxid, species_taxid,refseq_category,ftp_path,sep="\t",file=out2)
                    print(line, sep="\t", file=out2)
out.close()
out2.close()