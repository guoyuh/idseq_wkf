#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import re
import os
import sys
import pandas as pd
import subprocess
from path import Path

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
    if os.path.exists(single_kneaddata_sum):
        #/mnt/home/huanggy/project/20211203_mNGS_2/result/stats/kneaddata.sum.txt
        tb = pd.read_table(single_kneaddata_sum,sep="\t",header=0)
        raw_reads = int(tb['raw single'].values[0])
        return raw_reads
    else:
        pass



def get_category_counts(before_filter_mpa_report):
    """
    :param before_filter_mpa_report: kreport2mpa ---> before_filter_mpa_report
    :return:
    """
    Eukaryota_amount = 0
    fungi_amount = 0
    bacteria_amount = 0
    viruses_amount = 0
    fungi_organism_name = is_fungi(fungi_assemble_refseq="/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt")

    with open(before_filter_mpa_report,"r") as fh:
        for line in fh:
            line = line.strip()
            #fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            #Eukaryota_re = re.match('k__Eukaryota\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)

            #真核k__Eukaryota
            Eukaryota_re = re.match('k__Eukaryota\t(\d+)', line)
            if Eukaryota_re:
                Eukaryota_amount = int(Eukaryota_re.group(1))

            #fh.seek(0,0)
            ##所有种
            species_re = re.search(r'\|s__(\w+)\t(\d+)', line)
            if species_re and line.startswith("k__Eukaryota"):
                #寄生虫=真核-真菌
                species_name = species_re.group(1).replace("_", " ", 2)

                #真集中真菌判断
                for k in fungi_organism_name.keys():
                    if k.find(species_name) != -1:
                        fungi_amount += int(species_re.group(2))
                        ### 如果不break 掉，比如同种有多个菌株Aspergillus niger CBS 513.88，Aspergillus niger CBS 101883 ，这样Aspergillus niger 就计算多次
                        break

            # 细菌 k__Bacteria
            bacteria_re = re.match('k__Bacteria\t(\d+)', line)
            if bacteria_re:
                bacteria_amount = int(bacteria_re.group(1))
                print('bacteria_amount', bacteria_amount)
            ##病毒
            viruses_re = re.match('k__Viruses\t(\d+)', line)
            if viruses_re:
                viruses_amount = int(viruses_re.group(1))
                print('viruses_amount', viruses_amount)
    # 寄生虫=真核-真菌
    parasites_amount = Eukaryota_amount - fungi_amount
    print("fungi_amount:",fungi_amount,"bacteria_amount:",bacteria_amount,"parasites_amount:",parasites_amount,"viruses_amount:",viruses_amount)
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

"""
#添加taxid
./tsv-join -z -H --filter-file /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy_2.txt  --key-fields level,Latin_name  --append-fields  taxid /mnt/home/huanggy/project/20211216/result/RUN1_LX003/kraken2/RUN1_LX003.uniq_kreport.mpa.rpm.xls
"""

#before_filter_mpa_report  计算 kingdom_reads  为了计算基于kingdom 的相对丰度


# kingdom = {
#     Bacteria :0,
#     Fungi:0,
#     Parasite:0
#     viral:0
#
# }

# "Genus_name\tGenus_reads\tG_rpm\tGenus_relative_Abundance"  正则匹配
# "Latin_name\tLatin_reads\tlineage\trelative_Abundance\tRPM\tRP20M"  正则匹配
# "Coverage(ratio)\tDepth\tBSD-conclusion\tFoldChange(BSD)"

#Genus 建立字典 ，遍历species
def get_Genus_counts(filter_mpa_report):
    genus_cnt = {}
    with open(filter_mpa_report, "r") as fh:
        for line in fh:
            if not line :
                break
            line = line.strip()
            genus_re = re.search(r'\|g__(\w+)\t(\d+)', line)

            ##针对病毒，世系可能没有genus 这一node 它的super—node 就是family
            family_re = re.search(r'\|f__(\w+)\t(\d+)', line)
            if genus_re:
                genus_name = genus_re.group(1).replace("_"," ",2)
                genus_counts = int(genus_re.group(2))

                if genus_name not in genus_cnt:
                    ##一个genus 可以有多个speciess
                    genus_cnt[genus_name] = genus_counts
                else:
                    pass
            elif family_re:
                viral_family_name = family_re.group(1).replace("_"," ",2)
                viral_family_counts = int(family_re.group(2))
                if viral_family_name not in genus_cnt:
                    genus_cnt[viral_family_name] = viral_family_counts
                else:
                    pass

            else:
                pass
    return genus_cnt

def is_fungi(fungi_assemble_refseq="/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt"):
    fungi_organism_name = {}
    with open(fungi_assemble_refseq) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            organism_name = line.strip().split("\t")[7]
            if organism_name not in fungi_organism_name:
                fungi_organism_name[organism_name] = 1
            else:
                pass
    return fungi_organism_name


def judge_fungi(oganism):
    fungi_assemble_refseq = "/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt"
    with open(fungi_assemble_refseq) as fh:
        for line in fh:
            if line.find(oganism) != -1:
                return True


def negative_control(workdir,single_kneaddata_sum):
    tmp = Path(workdir).glob('result/*/kraken2/*.report.mpa.txt')
    negative_Species_set = {}
    for f in tmp:
        #print(f.basename())
        if f.basename().find("yinkong") != -1:
            print("find ======>:", f.basename())

            negative_filter_mpa_report = f
            with open(negative_filter_mpa_report, "r") as fh:
                for line in fh:
                    if not line :
                        break
                    line = line.strip()
                    species_re = re.search(r'\|s__(\w+)\t(\d+)', line)

                    if species_re:
                        species_name = species_re.group(1).replace("_", " ", 2)
                        Species_reads = int(species_re.group(2))
                        raw_single = get_rawReads(single_kneaddata_sum=single_kneaddata_sum)
                        Species_rpm = "{:.2f}".format(Species_reads / raw_single * 20 * 1.0e6)

                        if species_name not in negative_Species_set:
                            negative_Species_set[species_name] = {"Species_reads":Species_reads,"Species_rpm":Species_rpm}
                        else:
                            pass

            break
    return negative_Species_set

#os.path.join(os.path.dirname(filter_mpa_report),os.path.basename(filter_mpa_report).strip(".txt") + ".rpm.xls")
def read_mpa_report(single_kneaddata_sum,before_filter_mpa_report,filter_mpa_report,workdir,prefix):
    #fungi_organism_name = is_fungi(fungi_assemble_refseq="/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt")
    #print(fungi_organism_name)
    fungi_amount, bacteria_amount, parasites_amount, viruses_amount = get_category_counts(before_filter_mpa_report)
    genus_cnt = get_Genus_counts(filter_mpa_report)

    ##阴控样本的标准化reads 必须与阴控样本的single_kneaddata_sum 绑定在一起,##或者与阴控样本的clean_fastq
    tmp_kneaddata = Path(workdir).glob('result/*/filter/kneaddata.sum.txt')
    if len(tmp_kneaddata) >=1:
        for f in tmp_kneaddata:
            #print(f.basename())
            if f.find("yinkong") != -1:
                print("find yinkong_single_kneaddata_sum======>:", f)
                negative_Species_set = negative_control(workdir=workdir,single_kneaddata_sum=f)
    #print(genus_cnt)
    total_counts = get_rawReads(single_kneaddata_sum)
    #output_file = os.path.join(os.path.dirname(filter_mpa_report),prefix + "uniq_kreport.mpa.rpm.xls2")
    output_file = os.path.join(os.path.dirname(filter_mpa_report),os.path.basename(filter_mpa_report).strip(".txt") + ".rpm.xls")
    out = open(output_file, "w", encoding="gb2312")
    #sample_name = os.path.basename(filter_mpa_report).split(".")[0]
    sample_name=prefix
    header = ["#sample","category","category_reads","Genus_name","Genus_reads","G_rpm","G_relative_Abundance","Latin_name","Latin_reads","lineage","S_relative_Abundance","RP20M","BSD-conclusion","FoldChange(BSD)"]
    print("\t".join(header),file=out)
    with open(filter_mpa_report, "r") as fh:
        for line in fh:
            if not line :
                break
            line = line.strip()
            # fungi_re = re.match('k__Eukaryota\|k__Fungi\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)
            # if fungi_re:
            #     category = "Fungi"
            #     category_amount = fungi_amount
            #     lineage = fungi_re.group()
            #     Genus_name =fungi_re.group(1).replace("_"," ",2)
            #     Species_name = fungi_re.group(2).replace("_"," ",2)
            #     Species_reads = int(fungi_re.group(3))
            #
            #     Genus_reads = genus_cnt.get(Genus_name,0)
            #     Genus_rpm = Genus_reads / total_counts * 20*1.0e6   ##( (Genus_reads * 20*1.0e6 )/total_counts )
            #     Genus_abundance = "{:.2f}".format(Genus_reads / category_amount)
            #     Species_abundance = "{:.2f}".format(Species_reads / category_amount )
            #     Species_rpm = Species_reads / total_counts * 20*1.0e6
            #     print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,"-","-","-","-",sep="\t",file=out)


            # parasites_re = re.match('k__Eukaryota\|.*?\|([kpcofgs])__(\w+)\t(\d+)', line)
            # if (not parasites_re is None) & (fungi_re is None):
            #     category = "Parasite"
            #     category_amount = parasites_amount
            #     parasites_re2 = re.match('k__Eukaryota\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)
            #     if parasites_re2:
            #         lineage = parasites_re2.group()
            #         Genus_name = parasites_re2.group(1).replace("_", " ", 2)
            #         Species_name = parasites_re2.group(2).replace("_", " ", 2)
            #         Species_reads = int(parasites_re2.group(3))
            #
            #         Genus_reads = genus.get(Genus_name, 0)
            #         Genus_rpm = Genus_reads / total_counts * 20 * 1.0e6  ##( (Genus_reads * 20*1.0e6 )/total_counts )
            #         Genus_abundance = "{:.2f}".format(Genus_reads / category_amount)
            #         Species_abundance = "{:.2f}".format(Species_reads / category_amount)
            #         Species_rpm = Species_reads / total_counts * 20 * 1.0e6
            #         print(sample_name, category, category_amount, Genus_name, Genus_reads, Genus_rpm, Genus_abundance,Species_name, Species_reads, lineage, Species_abundance, Species_rpm, "-", "-", "-", "-",sep="\t", file=out)

            ##匹配到species 的行
            Eukaryota_re = re.match('k__Eukaryota\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)
            if Eukaryota_re:
                #真核生物
                species_re = re.search(r'\|s__(\w+)\t(\d+)', line)
                if species_re:
                    species_name = species_re.group(1).replace("_", " ", 2)
                    if judge_fungi(species_name):
                        ###真核生物里面的 真菌fungi
                        category = "Fungi"
                        category_amount = fungi_amount
                        lineage = Eukaryota_re.group().strip("\t" + Eukaryota_re.group(3))
                        Genus_name =Eukaryota_re.group(1).replace("_"," ",2)
                        Species_name = Eukaryota_re.group(2).replace("_"," ",2)
                        Species_reads = int(Eukaryota_re.group(3))

                        Genus_reads = genus_cnt.get(Genus_name,0)
                        Genus_rpm = round(Genus_reads / total_counts * 20*1.0e6)   ##( (Genus_reads * 20*1.0e6 )/total_counts )
                        Genus_abundance = "{:.3f}".format(Genus_reads / category_amount)
                        Species_abundance = "{:.3f}".format(Species_reads / category_amount )
                        Species_rpm = "{:.2f}".format(Species_reads / total_counts * 20*1.0e6)

                        if Species_name in negative_Species_set:
                            BSD_conclusion = "Y"
                            BSD_FoldChange = "{:.2f}".format(float(Species_rpm)/float(negative_Species_set[Species_name]["Species_rpm"]))
                        else:
                            BSD_conclusion = "N"
                            BSD_FoldChange = "-"
                        print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,BSD_conclusion,BSD_FoldChange,sep="\t",file=out)

                    else:
                        #parasites = 真核生物-真菌fungi
                        category = "Parasite"
                        category_amount = parasites_amount

                        lineage = Eukaryota_re.group().strip("\t" + Eukaryota_re.group(3))
                        Genus_name =Eukaryota_re.group(1).replace("_"," ",2)
                        Species_name = Eukaryota_re.group(2).replace("_"," ",2)
                        Species_reads = int(Eukaryota_re.group(3))

                        Genus_reads = genus_cnt.get(Genus_name,0)
                        Genus_rpm = round(Genus_reads / total_counts * 20*1.0e6)   ##( (Genus_reads * 20*1.0e6 )/total_counts )
                        Genus_abundance = "{:.3f}".format(Genus_reads / category_amount)
                        Species_abundance = "{:.3f}".format(Species_reads / category_amount )
                        Species_rpm = "{:.2f}".format(Species_reads / total_counts * 20*1.0e6)
                        if Species_name in negative_Species_set:
                            BSD_conclusion = "Y"
                            BSD_FoldChange = "{:.2f}".format(float(Species_rpm)/float(negative_Species_set[Species_name]["Species_rpm"]))
                        else:
                            BSD_conclusion = "N"
                            BSD_FoldChange = "-"
                        print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,BSD_conclusion,BSD_FoldChange,sep="\t",file=out)

            bacteria_re = re.match('k__Bacteria\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)
            if bacteria_re:
                category = "Bacteria"
                category_amount = bacteria_amount
                lineage = bacteria_re.group().strip("\t" +bacteria_re.group(3))
                Genus_name =bacteria_re.group(1).replace("_"," ",2)
                Species_name = bacteria_re.group(2).replace("_"," ",2)
                Species_reads = int(bacteria_re.group(3))

                Genus_reads = genus_cnt.get(Genus_name,0)
                Genus_rpm = round(Genus_reads / total_counts * 20*1.0e6)   ##( (Genus_reads * 20*1.0e6 )/total_counts )
                Genus_abundance = "{:.3f}".format(Genus_reads / category_amount)
                Species_abundance = "{:.3f}".format(Species_reads / category_amount )
                Species_rpm = "{:.2f}".format(Species_reads / total_counts * 20*1.0e6)
                if Species_name in negative_Species_set:
                    BSD_conclusion = "Y"
                    BSD_FoldChange = "{:.2f}".format(float(Species_rpm) / float(negative_Species_set[Species_name]["Species_rpm"]))
                else:
                    BSD_conclusion = "N"
                    BSD_FoldChange = "-"
                print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,BSD_conclusion,BSD_FoldChange,sep="\t",file=out)

            viruses_re = re.match('k__Viruses\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)
            if viruses_re:
                category = "viral"
                category_amount = viruses_amount
                lineage = viruses_re.group().strip("\t" +viruses_re.group(3))
                Genus_name =viruses_re.group(1).replace("_"," ",2)
                Species_name = viruses_re.group(2).replace("_"," ",2)
                Species_reads = int(viruses_re.group(3))

                Genus_reads = genus_cnt.get(Genus_name,0)
                Genus_rpm = round(Genus_reads / total_counts * 20*1.0e6)   ##( (Genus_reads * 20*1.0e6 )/total_counts )
                Genus_abundance = "{:.3f}".format(Genus_reads / category_amount)
                Species_abundance = "{:.3f}".format(Species_reads / category_amount )
                Species_rpm = "{:.2f}".format(Species_reads / total_counts * 20*1.0e6)
                if Species_name in negative_Species_set:
                    BSD_conclusion = "Y"
                    BSD_FoldChange = "{:.2f}".format(float(Species_rpm)/float(negative_Species_set[Species_name]["Species_rpm"]))
                else:
                    BSD_conclusion = "N"
                    BSD_FoldChange = "-"
                print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,BSD_conclusion,BSD_FoldChange,sep="\t",file=out)

            viruses_re2 = re.match('k__Viruses\|.*?\|f__(\w+)\|s__(\w+)\t(\d+)', line)
            if viruses_re2:
                category = "viral"
                category_amount = viruses_amount
                lineage = viruses_re2.group().strip("\t" +viruses_re2.group(3))
                Genus_name =viruses_re2.group(1).replace("_"," ",2)
                Species_name = viruses_re2.group(2).replace("_"," ",2)
                Species_reads = int(viruses_re2.group(3))

                Genus_reads = genus_cnt.get(Genus_name,0)
                Genus_rpm = Genus_reads / total_counts * 20*1.0e6   ##( (Genus_reads * 20*1.0e6 )/total_counts )
                Genus_abundance = "{:.3f}".format(Genus_reads / category_amount)
                Species_abundance = "{:.3f}".format(Species_reads / category_amount )
                Species_rpm = "{:.2f}".format(Species_reads / total_counts * 20*1.0e6)
                if Species_name in negative_Species_set:
                    BSD_conclusion = "Y"
                    BSD_FoldChange = "{:.2f}".format(float(Species_rpm) / float(negative_Species_set[Species_name]["Species_rpm"]))
                else:
                    BSD_conclusion = "N"
                    BSD_FoldChange = "-"
                print(sample_name,category,category_amount,Genus_name,Genus_reads,Genus_rpm,Genus_abundance,Species_name,Species_reads,lineage,Species_abundance,Species_rpm,BSD_conclusion,BSD_FoldChange,sep="\t",file=out)


    out.close()


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("usage: python {} single_kneaddata_sum  before_filter_mpa_report filter_mpa_report workdir prefix".format(sys.argv[0]))
        print(
            """
            example:\n
            python  {} /mnt/home/huanggy/project/20220408/result/21JS944006/filter/kneaddata.sum.txt /mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.kreport.mpa.txt /mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.uniq80_kreport.mpa.txt /mnt/home/huanggy/project/20220325 21JS944006_2
            """.format(sys.argv[0])
              )

    else:
        read_mpa_report(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
























