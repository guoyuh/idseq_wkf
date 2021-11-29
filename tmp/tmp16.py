import re,os

mpa_report = "/mnt/home/huanggy/project/20211101/result/LX2110582/kraken2/LX2110582.mpa.txt"


# with open(mpa_report,"r") as fh:
#     for line in fh:
#         line = line.strip()
#         Fungi_total_cnt = re.match("^k__Eukaryota\|k__Fungi\t(\d+)$",line)
#         if Fungi_total_cnt:
#             #Fungi_total_cnt.group(0) 相当于全词
#             print("=========================")
#             print(Fungi_total_cnt.group(1),"\n",line)


#         fungi_re2 = re.match('k__Eukaryota\|k__Fungi\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
#         if fungi_re2:
#             print("=======================")
#             print(fungi_re2.group(0)) #匹配全部
#             print(fungi_re2.group(1)) #匹配类别 kpcofgs
#             print(fungi_re2.group(2)) #匹配类别学名
#             print(fungi_re2.group(3)) #匹配该类别reads 数
#             break



# import subprocess
# #cmd = $(())"cat /mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq | wc -l "
# #cmd2 = "echo $(($(cat /mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq | wc -l)/4))"
# def count_reads(fastq):
#     cmd = "echo $(($(zcat {} | wc -l)/4))".format(fastq) if fastq.endswith("gz") else "echo $(($(cat {} | wc -l)/4))".format(fastq)
#     print(cmd)
#     (status,res) = subprocess.getstatusoutput(cmd)
#     print(status,res)
#     return res
#
#
# count_reads(fastq="/mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq")



import matplotlib.pyplot as plt
import pandas as pd
data = pd.read_csv(
    "https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/3_TwoNumOrdered.csv",
    delim_whitespace=True
)
print(data)