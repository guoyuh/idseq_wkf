import os 
refseq_dir_bac = "/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back"
refseq_dir_fungi = "/mnt/data/NCBI_Refseq/fungi/my_fungi_refseq"
refseq_dir_protozoa = "/mnt/data/NCBI_Refseq/protozoa/my_protozoa_refseq"
Karius_refseq = "/mnt/data/kraken2_mask/library/Karius_refseq"
#print(os.listdir(refseq_dir))
files_list_bac = os.listdir(refseq_dir_bac)
files_list_fungi =  os.listdir(refseq_dir_fungi)
files_list_protozoa = os.listdir(refseq_dir_protozoa)
files_list_Karius_refseq = os.listdir(Karius_refseq)

out = open("/mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20220406.txt","w")
bac_num=0
fungi_num = 0
protozoa_num = 0 
Karius_refseq_num = 0

header = ["assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date"
]
print("\t".join(header),file=out)
for file in files_list_bac:
    file = file.strip("_genomic.fna")
    with open ("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                if line.find(file) != -1:
                    bac_num +=1
                    print(line,file=out)



for file in files_list_fungi:
    file = file.strip("_genomic.fna")
    with open ("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                if line.find(file) != -1:
                    fungi_num += 1
                    print(line,file=out)



for file in files_list_protozoa:
    file = file.strip("_genomic.fna")
    with open ("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                if line.find(file) != -1:
                    protozoa_num += 1
                    print(line,file=out)



for file in files_list_Karius_refseq:
    file = file.strip("_genomic.fna")
    with open ("/mnt/data/NCBI_Refseq/assembly_summary_refseq.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                if line.find(file) != -1:
                    Karius_refseq_num += 1
                    print(line,file=out)



print("bac_num:",bac_num)
print("fungi_num:",fungi_num)
print("protozoa_num:",protozoa_num)
print("Karius_refseq_num:",Karius_refseq_num)

out.close()
