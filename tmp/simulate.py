#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

#Staphylococcus_aureus
import random
import numpy as np
import path,os


print(np.random.randint(1,35269,20))


def generate_random_seq(random_seqnum=20):
    strain = "Staphylococcus_aureus"
    with open(f"/mnt/home/huanggy/project/20211026_simulator/result/filter/{strain}_{random_seqnum}.fa","w") as fh1:
        with open("/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_000013425.1_ASM1342v1_genomic.fna","r") as fh2:
            a = fh2.readlines()
            print(a[0],a[1],a[2],len(a))
            random_ind = np.random.randint(1,35269,random_seqnum)
            for i,ind in enumerate(random_ind):
                fh1.write(">" + f"{strain}_0000{i}" + "\n")
                start = random.randint(0, 30)
                stop = start + 49
                print(start, stop)
                seq = a[ind][start:stop]
                fh1.write(seq+"\n")


# generate_random_seq(random_seqnum=20)
# generate_random_seq(random_seqnum=200)
# generate_random_seq(random_seqnum=2000)


strain2refseq = {
    "Staphylococcus_aureus":"/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_000013425.1_ASM1342v1_genomic.fna",
    "Klebsiella_pneumoniae":"/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_000240185.1_ASM24018v2_genomic.fna",
    "Aspergillus_fumigatus":"/mnt/data/NCBI_Refseq/fungi/my_fungi_refseq/GCF_000002655.1_ASM265v1_genomic.fna"
}

def generate_random_fq(strain,out_dir ,random_seqnum=20):
    #strain = "Staphylococcus_aureus"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir,mode=0o777)
    strain_genome = strain2refseq.get(strain,"/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_000013425.1_ASM1342v1_genomic.fna")
    #strain_genome = os.path.basename(strain_genome)
    out_file = path.Path(out_dir)/f"{strain}_{random_seqnum}.fq"
    print(strain_genome,"\n",out_file)
    with open(out_file,"w") as fh1:
        with open(strain_genome,"r") as fh2:
            a = fh2.readlines()
            #print(a[0],a[1],a[2],len(a))
            random_ind = np.random.randint(1,35269,random_seqnum)
            for i,ind in enumerate(random_ind):
                fh1.write("@" + f"{strain}_0000{i}" + "\n")
                start = random.randint(0, 30)
                stop = start + 50
                #print(start, stop)
                seq = a[ind][start:stop]
                fh1.write(seq+"\n")
                fh1.write("+" + "\n")
                fh1.write("GACB+ICCICHCIC=CCCIC9C0=IC9FHCHC(IC*E>CC*0ID8CIGIB" + "\n")


# generate_random_fq(strain="Staphylococcus_aureus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=10)
# generate_random_fq(strain="Klebsiella_pneumoniae",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=10)
# generate_random_fq(strain="Aspergillus_fumigatus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=10)

generate_random_fq(strain="Staphylococcus_aureus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=20)
generate_random_fq(strain="Klebsiella_pneumoniae",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=20)
generate_random_fq(strain="Aspergillus_fumigatus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=20)


generate_random_fq(strain="Staphylococcus_aureus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=100)
generate_random_fq(strain="Klebsiella_pneumoniae",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=100)
generate_random_fq(strain="Aspergillus_fumigatus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=100)

generate_random_fq(strain="Staphylococcus_aureus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=200)
generate_random_fq(strain="Klebsiella_pneumoniae",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=200)

generate_random_fq(strain="Staphylococcus_aureus",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=2000)
generate_random_fq(strain="Klebsiella_pneumoniae",out_dir='/mnt/home/huanggy/project/simulateor',random_seqnum=2000)

