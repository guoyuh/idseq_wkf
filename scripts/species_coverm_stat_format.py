#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import os,sys
import pandas as pd

def is_organism(cove_d,species):
    if isinstance(cove_d,dict):
        if species in cove_d:
            return cove_d.get(species)

def run(workdir,prefix):
    my_refseq_Read_Count_cut_off = 1
    coverm_file = f'{workdir}/result/{prefix}/align/{prefix}.refseq.cove.xls'
    cove_stat = {}
    if os.path.exists(coverm_file):
        with open(coverm_file,'r') as fh:
            for line in fh:
                if line.startswith("Genome"):
                    continue
                if not line:
                    break
                else:
                    genome_name = line.strip().split("\t")[0]
                    my_refseq_Read_Count = int(line.strip().split("\t")[6])
                    if my_refseq_Read_Count < my_refseq_Read_Count_cut_off:
                        continue
                    if genome_name not in cove_stat:
                        ###相同的微生物种名 取其中一个，因亚种存在，有些微生物存在多个scaffold
                        cove_stat[genome_name] = line.strip().split("\t")[1:]
                    else:
                        pass

    uniq_kreport_mpa = f'{workdir}/result/{prefix}/kraken2/{prefix}.uniq_kreport.mpa.rpm.xls'
    out = open(f'{workdir}/result/{prefix}/kraken2/{prefix}.uniq_kreport.mpa.rpm.stat.xls','w')
    tb = pd.read_table(uniq_kreport_mpa,sep='\t',header=0)
    header = tb.columns.to_list()
    print("\t".join(header),'Coverage(ratio)','mean depth',file=out ,sep='\t')
    tb2 = tb.drop_duplicates(keep='first')
    del tb
    for ind,r in tb2.iterrows():
        ks = r['Latin_name']
        row = [r['#sample'],r['category'],str(r['category_reads']),r['Genus_name'],str(r['Genus_reads']),str(r['G_rpm']),
               str(r['G_relative_Abundance']),r['Latin_name'],str(r['Latin_reads']),r['lineage'],str(r['S_relative_Abundance']),
               str(r['RP20M']),r['BSD-conclusion'],r['FoldChange(BSD)']]

        if ks in cove_stat:
            print("\t".join(row), cove_stat[ks][2], cove_stat[ks][0], sep='\t', file=out)
        else:
            print("\t".join(row), '-', '-', sep='\t', file=out)
    out.close()


def run2(k_report,refseq_cove):
    my_refseq_Read_Count_cut_off = 1
    cove_stat = {}
    with open(refseq_cove,'r') as fh:
        for line in fh:
            if line.startswith("Genome"):
                continue
            if not line:
                break
            else:
                genome_name = line.strip().split("\t")[0]
                my_refseq_Read_Count = int(line.strip().split("\t")[6])
                if my_refseq_Read_Count < my_refseq_Read_Count_cut_off:
                    continue
                if genome_name not in cove_stat:
                    ###相同的微生物种名 取其中一个，因亚种存在，有些微生物存在多个scaffold
                    cove_stat[genome_name] = line.strip().split("\t")[1:]
                else:
                    pass


    out_file = k_report.rstrip(".xls") + ".stat.xls"
    out = open(out_file,'w')
    tb = pd.read_table(k_report,sep='\t',header=0)
    header = tb.columns.to_list()
    print("\t".join(header),'Coverage(ratio)','mean depth',file=out ,sep='\t')
    tb2 = tb.drop_duplicates(keep='first')
    del tb
    for ind,r in tb2.iterrows():
        ks = r['Latin_name']
        row = [r['#sample'],r['category'],str(r['category_reads']),r['Genus_name'],str(r['Genus_reads']),str(r['G_rpm']),
               str(r['G_relative_Abundance']),r['Latin_name'],str(r['Latin_reads']),r['lineage'],str(r['S_relative_Abundance']),
               str(r['RP20M']),r['BSD-conclusion'],r['FoldChange(BSD)']]

        if ks in cove_stat:
            print("\t".join(row), cove_stat[ks][2], cove_stat[ks][0], sep='\t', file=out)
        else:
            print("\t".join(row), '-', '-', sep='\t', file=out)
    out.close()




if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage example: python /mnt/home/huanggy/project/20220428/result/21JS944006/kraken2/21JS944006.uniq_kreport.mpa.rpm.xls  /mnt/home/huanggy/project/20220428/result/21JS944006/align/21JS944006.refseq.cove.xls ".format(sys.argv[0]))
        exit(1)
    else:
        run2(k_report=sys.argv[1],refseq_cove=sys.argv[2])