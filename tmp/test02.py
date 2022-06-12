#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import datetime
import time

# def fun1():
#     start = datetime.datetime.now()
#     time.sleep(6)
#     end = datetime.datetime.now()
#     print("total time is:%s"%(end-start))


# fun1()

import pysam
# fa_file='/mnt/data/kraken2_db/bacteria_fungi_protozoa_viral.mask.fna'
# fa = pysam.FastaFile(fa_file)


# print(len(fa.fetch("NC_000913.3")))

# print(fa.fetch("NC_000913.3")[10])


# print("names of reference sequences: " + ",".join(fa.references))

total = 0 
with pysam.FastxFile("/mnt/data/NCBI_Refseq/bac_fungi_protozoa_viral/GCF_003205575.1_ASM320557v1_genomic.fna") as fh:
    for record in fh:
        print(record.name)
        print(record.comment)
        # print(type(record.sequence))
        print(len(record.sequence))
        total += len(record.sequence)
        
        #print(record.quality)
print(total)
    
