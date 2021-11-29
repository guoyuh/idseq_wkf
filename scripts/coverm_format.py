#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import os,sys
import pandas as pd
import numpy as np

header = ["assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date"]
assemble_summary = "/mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20211101.txt"
def assembleFile2strain(assemble_summary):
    assemble_genome2organism_name = {}
    with open(assemble_summary) as fh:
        for line in fh:
            if line.startswith("assembly_accession"):
                continue
            if not line:
                break
            else:
                organism_name = line.strip().split("\t")[7]
                ftp_path = line.strip().split('\t')[19].split('/')[-1] + "_genomic"
                #print(organism_name + "------>" + ftp_path)
                if ftp_path not in organism_name:
                    assemble_genome2organism_name[ftp_path] = organism_name
                else:
                    print("assemble_genome_name 有重复")
    return assemble_genome2organism_name


def assembleFile2strain_update_viral(assemble_summary,kraken2_viral_refseq_info):
    assemble_genome2organism_name = assembleFile2strain(assemble_summary)
    with open(kraken2_viral_refseq_info) as fh:
        for line in fh:
            if not line :
                break
            else:
                organism_name = line.strip().split("\t")[1]
                file_name = line.strip().split("\t")[0].split("|")[-1] + "_genomic"
                if file_name not in assemble_genome2organism_name :
                    assemble_genome2organism_name[file_name] = organism_name
                else:
                    pass

    return  assemble_genome2organism_name


def run(coverm_file,assemble_summary):
    assemble_genome2organism_name = assembleFile2strain(assemble_summary)
    tb = pd.read_table(coverm_file,sep='\t',header=0)
    out = os.path.dirname(coverm_file) + '/' + os.path.basename(coverm_file).strip(".cove.tsv") + '.refseq.cove.tsv'
    tb.columns = ['Genome', 'Covered Fraction', 'Variance', 'Length', 'my_refseq_Read_Count', 'RPKM']
    # for i in range(len(tb)):
    #     print(tb.loc[i, 'Genome'])
    #     #if assemble_genome2organism_name.get(tb.loc[i, 'Genome']):
    #         tb.loc[i, 'Genome'] = assemble_genome2organism_name.get(tb.loc[i, 'Genome'])
    tb['Genome'] = tb['Genome'].map(assemble_genome2organism_name,na_action='[None]')

    ######去除没有refseq 即为NaN 的行
    ##方法 一
    # tmp = tb[['Covered Fraction', 'Variance', 'Length', 'my_refseq_Read_Count', 'RPKM']]
    # tb = tb[tmp.apply(np.sum, axis=1) != 0]
    ##方法 二
    tb.dropna(axis=0,how='any') ####how = 'any' ,该轴方向所有为NaN 的行删除；how = 'all' ,该轴方向存在空值，即删除该行
    tb.to_csv(out,sep='\t',index=False)

def run2(coverm_file,assemble_summary,kraken2_viral_refseq_info):
    my_refseq_Read_Count_cut_off = 1
    assemble_genome2organism_name = assembleFile2strain_update_viral(assemble_summary,kraken2_viral_refseq_info)
    header = ['Genome', 'Covered Fraction', 'Variance', 'Length', 'my_refseq_Read_Count', 'RPKM']
    out = os.path.dirname(coverm_file) + '/' + os.path.basename(coverm_file).strip(".cove.tsv") + '.refseq.cove.tsv'
    myout = open(out,'w')
    print("\t".join(header),sep='\t',file=myout)
    with open(coverm_file,'r') as fh:
        for line in fh:
            if line.startswith("Genome"):
                continue
            if not line:
                break
            else:
                genome = line.strip().split("\t")[0]
                my_refseq_Read_Count = int(line.strip().split("\t")[4])
                if my_refseq_Read_Count < my_refseq_Read_Count_cut_off:
                    continue
                if genome in assemble_genome2organism_name:
                    print(assemble_genome2organism_name[genome],'\t'.join(line.strip().split("\t")[1:]),sep='\t',file=myout)
                else:
                    print(f'{genome} not exist in my refseq db，please add it ')
    myout.close()

# run2(coverm_file = '/mnt/home/huanggy/project/20211111/result/LX2110582/align/LX2110582.cove.tsv',
#      assemble_summary = '/mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20211101.txt',
#      kraken2_viral_refseq_info = '/mnt/data/NCBI_Refseq/viral/kraken2_viral_refseq_info')
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage example: python {} LX2110582.cove.tsv my_assembly_summary_refseq.info20211101.txt kraken2_viral_refseq_info ".format(sys.argv[0]))
        exit(1)
    else:
        run2(coverm_file=sys.argv[1],assemble_summary = sys.argv[2],kraken2_viral_refseq_info=sys.argv[3])