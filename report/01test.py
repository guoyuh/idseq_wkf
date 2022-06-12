#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import os,glob
import pandas as pd
from path import  Path
from docx import Document



in_path = '/mnt/home/huanggy/project/report/LX2110650-吴伟双.docx'
# 按二维矩阵索引遍历表格
doc = Document(in_path)
tables = doc.tables

sample_info_num=1
DNA_tb_num=3 #DNA 病毒
bac_tb_num=4
fungi_tb_num=5
pro_tb_num=6
tolarent_gene_tb_num=7
likely_bacground_tb_num=8
Test_result_description = 9
mytable = tables[Test_result_description]


res = []
print(dir(mytable))
sample_type = tables[sample_info_num].row_cells(1)[0].text
from_part = tables[sample_info_num].row_cells(2)[1].text
print(sample_type,from_part)


print(f'第{Test_result_description}个表格: {len(mytable.rows)} 行 X {len(mytable.columns)}列')
for i, row in enumerate(tables[Test_result_description].rows):
    print(f'第 {i+1} 行有 {len(row.cells)} 列')
    tmp = []
    for j, cell in enumerate(row.cells):
        print(cell.text)
#
#
#
#
# for i, row in enumerate(tables[sample_info_num].rows):
#     for j, cell in enumerate(row.cells):
#         print(f"第{i+1}行,{j+1}列的value:{cell.text}")



