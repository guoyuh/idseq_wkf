#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


import pandas as pd


header = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
tb = pd.read_table("/mnt/home/huanggy/project/meta/result/align/my_gsnap.out",sep = '\t',header=None)
print(tb.head(10))
tb.columns = header
print(tb.columns)
print(tb["sseqid"])


def sub_fun(sseqid):
    species_taxid = sseqid.split("|")[1]
    sub_seqid = sseqid.split("|")[2]
    return str(species_taxid)

#tb['species_taxid'] = tb["sseqid"].map(sub_fun)
tb['species_taxid'] = tb["sseqid"].apply(sub_fun)
print(tb.head())

print(tb.groupby("species_taxid")["qseqid"].value_counts())
print(tb["species_taxid"].describe())
