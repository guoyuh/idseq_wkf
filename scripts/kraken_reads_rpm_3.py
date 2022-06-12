#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import re
import os
import sys
import pandas as pd
import subprocess

### reads uniq  过滤之前 ，kraken_reads_rpm_2 计算mpa_report物种相对丰度，RPM,RPM20
####reads uniq  过滤之后 ，kraken_reads_rpm_3 计算过滤后mpa_report物种相对丰度，RPM,RPM20 

def count_reads(fastq):
    #cmd = "echo $(($(cat /mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq | wc -l)/4))"
    cmd = "echo $(($(zcat {} | wc -l)/4))".format(fastq) if fastq.endswith("gz") else "echo $(($(cat {} | wc -l)/4))".format(fastq)
    print(cmd)
    (status,res) = subprocess.getstatusoutput(cmd)
    print(status,res)
    return int(res)



def get_rawReads(single_kneaddata_sum):
    #/mnt/home/huanggy/project/20211203_mNGS_2/result/stats/kneaddata.sum.txt
    tb = pd.read_table(single_kneaddata_sum,sep="\t",header=0)
    raw_reads = int(tb['raw single'].values[0])
    return raw_reads

def get_category_counts(before_filter_mpa_report):
    """
    :param before_filter_mpa_report: kreport2mpa ---> before_filter_mpa_report
    :return:
    """
    fungi_amount = 0
    bacteria_amount = 0
    parasites_amount = 0
    viruses_amount = 0
    with open(before_filter_mpa_report,"r") as fh:
        for line in fh:
            line = line.strip()
            fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            if fungi_re:
                fungi_amount = int(fungi_re.group(1))

            #真核k__Eukaryota
            Eukaryota_re = re.match('k__Eukaryota\t(\d+)', line)
            if Eukaryota_re:
                # 寄生虫 (真核-真菌)
                parasites_amount = int(Eukaryota_re.group(1)) - fungi_amount
                print('parasites_amount', parasites_amount)

            # 细菌 k__Bacteria
            bacteria_re = re.match('k__Bacteria\t(\d+)', line)
            if bacteria_re:
                bacteria_amount = int(bacteria_re.group(1))
                print('bacteria_amount', bacteria_amount)

            viruses_re = re.match('k__Viruses\t(\d+)', line)
            if viruses_re:
                viruses_amount = int(viruses_re.group(1))
                print('viruses_amount', viruses_amount)
    return fungi_amount,bacteria_amount,parasites_amount,viruses_amount

"""
k__Eukaryota|k__Fungi   475
k__Eukaryota|k__Fungi|p__Ascomycota     373
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes  270
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales     270
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Phaffomycetaceae 145
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Phaffomycetaceae|g__Komagataella 145
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Phaffomycetaceae|g__Komagataella|s__Komagataella_phaffii 145
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae    61
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae|g__Yarrowia        61
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae|g__Yarrowia|s__Yarrowia_lipolytica 61
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae       50
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Eremothecium       23
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Eremothecium|s__Eremothecium_cymbalariae   21
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Eremothecium|s__Eremothecium_gossypii      2
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Torulaspora        7
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Torulaspora|s__Torulaspora_delbrueckii     6
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Torulaspora|s__Torulaspora_globosa 1
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Saccharomyces      5
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Saccharomyces|s__Saccharomyces_paradoxus   3
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Saccharomyces|s__Saccharomyces_cerevisiae  1
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Saccharomyces|s__Saccharomyces_eubayanus   1
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Kluyveromyces      4
k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Kluyveromyces|s__Kluyveromyces_marxianus   4
"""


def run(single_kneaddata_sum,before_filter_mpa_report,filter_mpa_report):
    fungi_amount, bacteria_amount, parasites_amount, viruses_amount = get_category_counts(before_filter_mpa_report)
    total_counts = get_rawReads(single_kneaddata_sum)
    #if os.path.exists(mpa_report):
    output_file = os.path.join(os.path.dirname(filter_mpa_report),os.path.basename(filter_mpa_report).strip(".txt") + ".rpm.xls")
    print(output_file)
    out = open(output_file,"w",encoding="gb2312")
    #out.write('类别\t类别reads数\t分类级别\t名称\tlineage\treads数\trelative_Abundance\tRPM\tRP20M\n')
    out.write('category\tcategory_reads\tlevel\tLatin_name\tlineage\tLatin_reads\trelative_Abundance\tRPM\tRP20M\n')
    with open(filter_mpa_report, 'r')as L:
        # for line in L:
        #     line = line.strip()
        #     fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
        #     if fungi_re:
        #         fungi_amount = int(fungi_re.group(1))
        # L.seek(0, 0)
        for line in L:
            line = line.strip()
            ##真菌 k__Eukaryota|k__Fungi
            # fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            # if fungi_re:
            #     fungi_amount = int(fungi_re.group(1))
            #     print('fungi_amount', fungi_amount)
            fungi_re2 = re.match('k__Eukaryota\|k__Fungi\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            # print(fungi_re2.group(0)) #匹配全部
            # print(fungi_re2.group(1)) #匹配类别 kpcofgs
            # print(fungi_re2.group(2)) #匹配类别学名
            # print(fungi_re2.group(3)) #匹配该类别reads 数
            if fungi_re2:
                if fungi_re2.group(1).upper() == "G" or fungi_re2.group(1).upper() == "S":
                    out.write('fungi\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}\t{6:.2f}\t{7:.2f}\n'.format(fungi_amount,
                                                                                   fungi_re2.group(1).upper(),
                                                                                   fungi_re2.group(2).replace("_"," ",2),
                                                                                   line.strip("\t" + fungi_re2.group(3)),
                                                                                   fungi_re2.group(3),
                                                                                   int(fungi_re2.group(3)) / fungi_amount,
                                                                                   int(fungi_re2.group(3)) / total_counts * 1.0e6,
                                                                                   int(fungi_re2.group(3)) / total_counts * 20*1.0e6))
                else:
                    pass

            ##真核k__Eukaryota
            # Eukaryota_re = re.match('k__Eukaryota\t(\d+)', line)
            # if Eukaryota_re:
            #     # 寄生虫 (真核-真菌)
            #     parasites_amount = int(Eukaryota_re.group(1)) - fungi_amount
            #     print('parasites_amount', parasites_amount)
            parasites_re = re.match('k__Eukaryota\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            #真核里面包含寄生虫的parasites_re ，但是不包含真菌的 fungi_re2 is None  ，并集
            #if (parasites_re is not None) & (fungi_re2 is None):
            if (not parasites_re is None) & (fungi_re2 is None):
                if parasites_re.group(1).upper() == "G" or parasites_re.group(1).upper() == "S":
                    out.write('parasites\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}\t{6:.2f}\t{7:.2f}\n'.format(parasites_amount,
                                                                                   parasites_re.group(1).upper(),
                                                                                   parasites_re.group(2).replace("_"," ",2),
                                                                                   line.strip("\t" + parasites_re.group(3)),
                                                                                   parasites_re.group(3),
                                                                                   int(parasites_re.group(3)) / parasites_amount,
                                                                                   int(parasites_re.group(3)) / total_counts * 1.0e6,
                                                                                   int(parasites_re.group(3)) / total_counts * 20 * 1.0e6))
                else:
                    pass

            ## 细菌 k__Bacteria
            # bacteria_re = re.match('k__Bacteria\t(\d+)', line)
            # if bacteria_re:
            #     bacteria_amount = int(bacteria_re.group(1))
            #     print('bacteria_amount', bacteria_amount)
            bacteria_re2 = re.match('k__Bacteria\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            if bacteria_re2:
                if bacteria_re2.group(1).upper() == "G" or bacteria_re2.group(1).upper() == "S":
                    out.write('bacteria\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}\t{6:.2f}\t{7:.2f}\n'.format(bacteria_amount,
                                                                                  bacteria_re2.group(1).upper(),
                                                                                  bacteria_re2.group(2).replace("_"," ",2),
                                                                                  line.strip("\t" + bacteria_re2.group(3)),
                                                                                  bacteria_re2.group(3),
                                                                                  int(bacteria_re2.group(3)) / bacteria_amount,
                                                                                  int(bacteria_re2.group(3)) / total_counts * 1.0e6,
                                                                                  int(bacteria_re2.group(3)) / total_counts * 20 * 1.0e6))
                else:
                    pass

            # # 病毒 k__Viruses
            # viruses_re = re.match('k__Viruses\t(\d+)', line)
            # if viruses_re:
            #     viruses_amount = int(viruses_re.group(1))
            #     print('viruses_amount', viruses_amount)
            viruses_re2 = re.match('k__Viruses\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            if viruses_re2:
                if viruses_re2.group(1).upper() == "G" or viruses_re2.group(1).upper() == "S" :
                    out.write('viruses\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}\t{6:.2f}\t{7:.2f}\n'.format(viruses_amount,
                                                                                     viruses_re2.group(1).upper(),
                                                                                     viruses_re2.group(2).replace("_"," ",2),
                                                                                     line.strip("\t" + viruses_re2.group(3)),
                                                                                     viruses_re2.group(3),
                                                                                     int(viruses_re2.group(3)) / viruses_amount,
                                                                                     int(viruses_re2.group(3)) / total_counts * 1.0e6,
                                                                                     int(viruses_re2.group(3)) / total_counts * 20 * 1.0e6))
                else:
                    pass
    out.close()

"""
#添加taxid
./tsv-join -z -H --filter-file /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy_2.txt  --key-fields level,Latin_name  --append-fields  taxid /mnt/home/huanggy/project/20211216/result/RUN1_LX003/kraken2/RUN1_LX003.uniq_kreport.mpa.rpm.xls
"""

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("usage: python {} single_kneaddata_sum  before_filter_mpa_report filter_mpa_report ".format(sys.argv[0]))
    else:
        run(sys.argv[1],sys.argv[2],sys.argv[3])























