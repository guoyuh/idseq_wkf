#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import re
import os
import sys
import subprocess

def count_reads(fastq):
    #cmd = "echo $(($(cat /mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq | wc -l)/4))"
    cmd = "echo $(($(zcat {} | wc -l)/4))".format(fastq) if fastq.endswith("gz") else "echo $(($(cat {} | wc -l)/4))".format(fastq)
    print(cmd)
    (status,res) = subprocess.getstatusoutput(cmd)
    print(status,res)
    return int(res)

def run(fastq,mpa_report):
    fungi_amount = 0
    bacteria_amount = 0
    parasites_amount = 0
    viruses_amount = 0
    #total_counts = 0
    #去除接头，宿主之后的total reads
    total_counts = count_reads(fastq)
    #if os.path.exists(mpa_report):
    output_file = os.path.join(os.path.dirname(mpa_report),os.path.basename(mpa_report).strip(".txt") + ".rpm.csv")
    print(output_file)
    out = open(output_file,"w",encoding="gb2312")
    out.write('类别,类别reads数,分类级别,名称,lineage,reads数,丰度,RPM\n')
    with open(mpa_report, 'r')as L:

        for line in L:
            line = line.strip()
            fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            if fungi_re:
                fungi_amount = int(fungi_re.group(1))
        L.seek(0, 0)
        for line in L:
            line = line.strip()
            # 真菌 k__Eukaryota|k__Fungi
            fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            if fungi_re:
                fungi_amount = int(fungi_re.group(1))
                print('fungi_amount', fungi_amount)
            fungi_re2 = re.match('k__Eukaryota\|k__Fungi\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            # print(fungi_re2.group(0)) #匹配全部
            # print(fungi_re2.group(1)) #匹配类别 kpcofgs
            # print(fungi_re2.group(2)) #匹配类别学名
            # print(fungi_re2.group(3)) #匹配该类别reads 数
            if fungi_re2:
                out.write('fungi,{0},{1},{2},{3},{4},{5:.2%},{6:.2f}\n'.format(fungi_amount,
                                                                               fungi_re2.group(1),
                                                                               fungi_re2.group(2),
                                                                               line,
                                                                               fungi_re2.group(3),
                                                                               int(fungi_re2.group(3)) / fungi_amount,
                                                                               int(fungi_re2.group(3)) / total_counts * 1.0e6))
            #真核k__Eukaryota
            Eukaryota_re = re.match('k__Eukaryota\t(\d+)', line)
            if Eukaryota_re:
                parasites_amount = int(Eukaryota_re.group(1)) - fungi_amount
                print('parasites_amount', parasites_amount)

            # 寄生虫 (真核-真菌)
            parasites_re = re.match('k__Eukaryota\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            #真核里面包含寄生虫的parasites_re ，但是不包含真菌的 fungi_re2 is None  ，并集
            #if (parasites_re is not None) & (fungi_re2 is None):

            if (not parasites_re is None) & (fungi_re2 is None):
                out.write('parasites,{0},{1},{2},{3},{4},{5:.2%},{6:.2f}\n'.format(parasites_amount,
                                                                               parasites_re.group(1),
                                                                               parasites_re.group(2),
                                                                               line,
                                                                               parasites_re.group(3),
                                                                               int(parasites_re.group(3)) / parasites_amount,
                                                                               int(parasites_re.group(3)) / total_counts * 1.0e6))

            # 细菌 k__Bacteria
            bacteria_re = re.match('k__Bacteria\t(\d+)', line)
            if bacteria_re:
                bacteria_amount = int(bacteria_re.group(1))
                print('bacteria_amount', bacteria_amount)
            bacteria_re2 = re.match('k__Bacteria\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            if bacteria_re2:
                out.write('bacteria,{0},{1},{2},{3},{4},{5:.2%},{6:.2f}\n'.format(bacteria_amount,
                                                                              bacteria_re2.group(1),
                                                                              bacteria_re2.group(2),
                                                                              line,
                                                                              bacteria_re2.group(3),
                                                                              int(bacteria_re2.group(3)) / bacteria_amount,
                                                                              int(bacteria_re2.group(3)) / total_counts * 1.0e6))

            # 病毒 k__Viruses
            viruses_re = re.match('k__Viruses\t(\d+)', line)
            if viruses_re:
                viruses_amount = int(viruses_re.group(1))
                print('viruses_amount', viruses_amount)
            viruses_re2 = re.match('k__Viruses\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            if viruses_re2:
                out.write('viruses,{0},{1},{2},{3},{4},{5:.2%},{6:.2f}\n'.format(viruses_amount,
                                                                                 viruses_re2.group(1),
                                                                                 viruses_re2.group(2),
                                                                                 line,
                                                                                 viruses_re2.group(3),
                                                                                 int(viruses_re2.group(3)) / viruses_amount,
                                                                                 int(viruses_re2.group(3)) / total_counts * 1.0e6))
    out.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("usage: python {} remove.host.fastq  mpa_report ".format(sys.argv[0]))
    else:
        run(sys.argv[1],sys.argv[2])























