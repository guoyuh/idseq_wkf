#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

from Bio import SeqIO
import re

outdir = '/mnt/data/NCBI_Refseq/viral/kraken2_refseq'
kraken2_refseq_info = '/mnt/data/NCBI_Refseq/viral/kraken2_refseq_info'
out = open(kraken2_refseq_info,'w+')
with open('/mnt/data/kraken2_mask/library/viral/library.fna','r') as fh:
    record = SeqIO.parse(fh, "fasta")
    for seq_record in record:
        #print(dir(seq_record))
        # print(seq_record.id)
        # print(seq_record.description)

        genome_accseeion_name  = seq_record.id.split("|")[-1]+"_genomic.fna"
        genome_accseeion_name_file = outdir + '/' + genome_accseeion_name
        # with open(genome_accseeion_name_file,'w') as h:
        #     h.write(">"+seq_record.description)
        #     h.write(seq_record.seq)
        #SeqIO.write(seq_record,genome_accseeion_name_file,"fasta")
        organism = seq_record.description.strip(seq_record.id).lstrip(" ")
        organism = re.sub(pattern=',.*',repl='',string=organism)
        # print(re.sub(pattern=',.*',repl='',string=organism))
        # print('==================================================>')
        out.write(seq_record.id + "\t" + organism + "\n")

out.close()

