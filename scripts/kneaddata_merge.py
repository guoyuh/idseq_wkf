#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang
import pandas as pd
from path import Path
import sys
def run(single_kneaddata_sum):
    tb = pd.read_table(single_kneaddata_sum,sep='\t',header=0)
    sample = tb.loc[0,"Sample"]
    raw_reads = tb.loc[0,"raw single"]
    #trimmed_single = tb.loc[0,'trimmed single']
    #filter_reads = tb.loc[0,'decontaminated yh_add_hg19_GRCh38_latest_rna single']
    final_single = tb.loc[0,'final single']
    suzhu_per = "{:.2f}".format((raw_reads - final_single)/raw_reads)
    #print("{:.2%}".format((raw_reads - filter_reads) / raw_reads))
    return sample,raw_reads,final_single,suzhu_per

def main(analysis_dir):
    p = Path(analysis_dir).iglob("result/*/filter/kneaddata.sum.txt")
    out_file = Path(analysis_dir)/'stat.txt'
    out = open(out_file,'w')
    print("Sample","raw_reads","filter_reads","suzhu_per",sep='\t',file=out)
    for f in p:
        sample, raw_reads, final_single, suzhu_per = run(f)
        print(sample, raw_reads, final_single, suzhu_per,sep='\t',file=out)
    out.close()

def run2(single_kneaddata_sum):
    tb = pd.read_table(single_kneaddata_sum, sep='\t', header=0)
    return pd.Series(tb.iloc[0,:])

def main2(analysis_dir):
    mm = pd.DataFrame()
    p = Path(analysis_dir).iglob("result/*/filter/kneaddata.sum.txt")
    out_file = Path(analysis_dir)/'stat.txt'
    for f in p:
        mm = mm.append(run2(f),ignore_index=False)
    print(mm)
    mm.sort_values(by='Sample',inplace=True,ascending=True)
    mm['suzhu_per'] = (mm['raw single'] - mm['final single'])/mm['raw single']
    mm['suzhu_per'] = mm['suzhu_per'].apply(lambda x:format(x,".4f"))
    mm = mm[["Sample","raw single","final single","suzhu_per"]]
    mm.columns = ["Sample","raw_reads","filter_reads","suzhu_per"]
    mm.to_csv(out_file,sep="\t",index=False)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: python3 %s analysis_dir"%(sys.argv[0]))
    else:
        main2(sys.argv[1])








