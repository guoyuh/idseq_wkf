#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd

blast = {'s__Klebsiella pneumoniae': {'lineage_id': '131567;2;1224;1236;91347;543;570;573;72407;1125630', 'lineage_name': 'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella pneumoniae|t__Klebsiella pneumoniae subsp. pneumoniae HS11286', 'count': 7711, 's_abundance': '28.60%'}, 's__Klebsiella variicola': {'lineage_id': '131567;2;1224;1236;91347;543;570;244366', 'lineage_name': 'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella variicola|t__unclassified Klebsiella variicola subspecies/strain', 'count': 564, 's_abundance': '2.09%'}}

tb = pd.DataFrame.from_dict(blast).transpose()
print(tb)
print(type(tb))
print("======>",tb.columns)
tb2 =tb.sort_values(by = "s_abundance",axis=0,ascending=False,inplace=False)
print(tb2)

