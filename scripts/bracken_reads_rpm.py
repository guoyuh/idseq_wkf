#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang
import pandas as pd
import argparse
import os

def count_clena_reads(clean_fastq):
    cnt = 0
    with open(clean_fastq,'r') as fh:
        for line in fh:
            if line.startswith("@"):
                cnt += 1
            else:
                pass
    return  cnt

# def npm(reads_num):
#     a = "%.2f" % (float(reads_num / incount * 1.0e6))
#     return  a

def float2(rpm):
    a = "%.2f" %(rpm)
    return  a


def run(workdir,prefix,cleanfq):
    #clean_fastq = f'{workdir}/result/{prefix}/filter/{prefix}.filter.fq'
    kraken_out =  f'{workdir}/result/{prefix}/kraken2/{prefix}.s.bracken'
    if os.path.exists(kraken_out):
        out = kraken_out + ".csv"
        tb = pd.read_table(kraken_out,sep = "\t")
        in_count = count_clena_reads(cleanfq)
        print(in_count)
        # def npm(reads_num):
        #     a = "%.2f" % (float(reads_num / in_count * 1.0e6))
        #     return a
        #tb['RPM'] = tb['kraken_assigned_reads'].apply(npm)
        tb['RPM_tmp'] = tb['kraken_assigned_reads'] / in_count * 1.0e6
        tb['RPM'] = tb['RPM_tmp'].map(float2)
        tb.drop('RPM_tmp',inplace=True,axis=1)
        tb.to_csv(out,sep=',',index=False)
    else:
        print("err: kraken out report no exist")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-w','--workdir',type=str,help='')
    parser.add_argument('-p','--prefix',type = str ,help='')
    parser.add_argument('-c', '--cleanfq', type=str, help='')
    args = parser.parse_args()
    run(workdir=args.workdir,prefix=args.prefix,cleanfq=args.cleanfq)


