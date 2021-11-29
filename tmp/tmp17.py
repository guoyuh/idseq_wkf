#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
tb = pd.read_table('/mnt/home/huanggy/idseq_dag/report/samples.txt',sep = "\t")
print(tb)

print(tb['sample_name'].to_list())