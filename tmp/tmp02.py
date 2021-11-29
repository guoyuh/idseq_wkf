#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
rt =pd.read_table("/mnt/data/NCBI/taxonomy//RefSeq_bac_fun_viral.lineage_db",sep="\t",header = 0)
print(rt.shape)
