#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang
import os,sys
def uniqreads_filter(kraken2_out,kmer_length = 34):
    output_file_1 = os.path.join(os.path.dirname(kraken2_out), os.path.basename(kraken2_out).strip(".output") + ".uniq_output")
    output_file_2 = os.path.join(os.path.dirname(kraken2_out),os.path.basename(kraken2_out).strip(".output") + ".uniq80_output")
    out100 = open(output_file_1,'w')
    out80 = open(output_file_2,'w')
    #uniq_match_kmer_counts = reads_length - kmer_length  #41
    with open(kraken2_out, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("U"):
                # 未分类
                continue
            else:
                taxid = line.split("\t")[2]
                reads_length = int(line.split("\t")[3])
                kmer_cly = line.split("\t")[4]
                if len(kmer_cly.split(" ")) > 1:
                    id_cnt_list = kmer_cly.split(" ")
                    for i in range(len(kmer_cly.split(" "))):
                        id = id_cnt_list[i].split(":")[0]
                        cnt = int(id_cnt_list[i].split(":")[1])
                        if id == taxid and cnt >= (reads_length - kmer_length) * 0.8:
                            ## %80 匹配
                            print(line,sep='\t',file=out80)
                else:
                    id = kmer_cly.split(":")[0]
                    cnt = int(kmer_cly.split(":")[1])
                    if id == taxid and cnt == (reads_length - kmer_length):
                        ## reads 唯一匹配
                        print(line, sep='\t', file=out100)
                        ## reads 唯一匹配显然也符合80% 匹配
                        print(line, sep='\t', file=out80)

    out80.close()
    out100.close()

##借用工具
#python ~/idseq_wkf/KrakenTools/make_ktaxonomy.py --nodes /mnt/data/kraken2_mask/taxonomy/nodes.dmp --names /mnt/data/kraken2_mask/taxonomy/names.dmp --seqid2taxid /mnt/data/kraken2_mask/seqid2taxid.map  -o  /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy.txt
#KrakenTools/make_kreport.py

###自己解析  RUN1_LX001.report.mpa.rpm.xls   纠正kraken_output 里面的reasd 数目（唯一匹配或者80%匹配）
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: python {} kraken2_out ".format(sys.argv[0]))
    else:
        uniqreads_filter(sys.argv[1])


