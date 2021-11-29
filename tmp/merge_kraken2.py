#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
import glob
import os

sp = {
    "S200026380_L01_90":"LX2110582",
    "S200026521_L01_102":"LX2110650",
    "S200026521_L01_101":"LX2110665",
    "S200023771_L01_102":"LX2110966",
    "S200023771_L01_101":"LX2110965",
    "S200023771_L01_127":"LX2110967"
}


file_list = glob.glob("/mnt/home/huanggy/project/20210927_test/result/kraken2/*.s.bracken")
print(file_list)
tb = pd.DataFrame()
for f in file_list:
    id = os.path.basename(f).replace(".s.bracken","")
    if id in sp:
        id = sp[id]
    else:
        pass
    print(id)
    tmp = pd.read_table(f, sep="\t",header=0)
    tmp["id"] = [id]*len(tmp)
    tb=tb.append(tmp,ignore_index=False)

tb = tb[["id","name","taxonomy_id","taxonomy_lvl","kraken_assigned_reads","added_reads","new_est_reads","fraction_total_reads"]]
tb.to_csv("/mnt/home/huanggy/project/20210927_test/result/kraken2/merge.s.bracken.csv",index=False)


file_list = glob.glob("/mnt/home/huanggy/project/20210927_test/result/kraken2/*.g.bracken")
print(file_list)
tb = pd.DataFrame()
for f in file_list:
    id = os.path.basename(f).replace(".g.bracken","")
    if id in sp:
        id = sp[id]
    else:
        pass
    print(id)
    tmp = pd.read_table(f, sep="\t",header=0)
    tmp["id"] = [id]*len(tmp)
    tb=tb.append(tmp,ignore_index=False)

tb = tb[["id","name","taxonomy_id","taxonomy_lvl","kraken_assigned_reads","added_reads","new_est_reads","fraction_total_reads"]]
tb.to_csv("/mnt/home/huanggy/project/20210927_test/result/kraken2/merge.g.bracken.csv",index=False)
