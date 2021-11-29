#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang
import re,os,sys
import datetime
import pandas as pd
from pathlib import Path


def extrat_fa_from_qseqid(output_path,qsseqid,top10idfile):
    pass




def dbasdict(db):
    db_dict = { }
    with open(db,'r') as fh:
        for line in fh:
            if not line:
                break
            if line.startswith("taxid"):
                continue
            taxid,species_name,accession_version,refseq_description,refseq_status,length,names_lineage,taxid_lineage,name,format_names_lineage = line.strip().split("\t")
            if accession_version not in db_dict:
                #db_dict[accession_version] =
                db_dict.update({accession_version:{"taxid":taxid,"taxid_lineage":taxid_lineage,"format_names_lineage":format_names_lineage}})
            else:
                pass

    return db_dict




def parser_blast4(blast_outm8,workdir,prefix):

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    start = datetime.datetime.now()
    print(start)
    mydb = dbasdict(db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")

    blast_out_s = {}
    blast_out_g = {}
    total_counts = 0
    with open(blast_outm8, "r") as blast_outm8_f:

        evalue_cutoff = float(1.0e-10)
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            if not line :
                break
            qseqid, subseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")

            #evalue 过滤
            if float(evalue) > evalue_cutoff:
                #print("true:",evalue)
                continue

            #subseqid is seq accession_num
            if subseqid.startswith("kraken"):
                subseqid = subseqid.split("|")[2]
            if subseqid in mydb:
                taxid = mydb[subseqid]["taxid"]
                lineage_id = mydb[subseqid]["taxid_lineage"]
                #lineage_name = mydb[sseqid]["format_names_lineage"]
                format_names_lineage = mydb[subseqid]["format_names_lineage"]
                s_name = re.match('.*\|(s__.*)\|.*', format_names_lineage).groups()[0]
                if s_name not in blast_out_s:
                    blast_out_s.update({s_name: {"taxid":taxid,"lineage_id": lineage_id, "lineage_name": format_names_lineage, "count": 1,"qseqid": [qseqid]}})
                    # blast_out_s.update({s_name: { "lineage_id": lineage_id,
                    #                              "lineage_name": format_names_lineage, "count": 1, "qseqid": [qseqid]}})
                else:
                    counts = blast_out_s[s_name]["count"] + 1
                    # blast_out.update({s_level}.update({"count":counts}))
                    blast_out_s[s_name].update({"count": counts})

                    ##同种 不同的reads append
                    #if qseqid not in blast_out_s[s_name]["qseqid"] and len(blast_out_s[s_name]["qseqid"])<3:
                    if qseqid not in blast_out_s[s_name]["qseqid"]:
                        blast_out_s[s_name]["qseqid"].append(qseqid)


                #### parser gensu#####
                for ind, i in enumerate(format_names_lineage.split("|")):
                    if i.startswith("g__"):
                        g_name = i
                        g_taxid = lineage_id.split(";")[ind - len(format_names_lineage.split("|"))]
                    else:
                        g_name = re.match('.*\|(g__.*)\|.*', format_names_lineage).groups()[0]
                        g_taxid = taxid
                    #my_g_name = re.match('.*\|(g__.*)\|.*', format_names_lineage).groups()[0]
                    if g_name not in blast_out_g:
                        blast_out_g.update({g_name: {"taxid":g_taxid,"lineage_id": lineage_id, "lineage_name": format_names_lineage, "count": 1, "qseqid": [qseqid]}})
                    else:
                        ##不同seqid 但是在同一个属  ,此时该属 +1 条hits
                        ##同一seqid 比对到不同的种水平 ，但是是同一个属 ，此时属于该属水平的 不能加1 ，属于重复reads
                        # blast_out.update({s_level}.update({"count":counts}))

                        # 同属 total count +1
                        #    同属不同种 ，该属水平 reads + 1
                        #    同属同种，该属水平 reads 不加1
                        if qseqid not in blast_out_g[g_name]["qseqid"]:
                            counts = blast_out_g[g_name]["count"] + 1
                            blast_out_g[g_name].update({"count": counts})
                            blast_out_g[g_name]["qseqid"].append(qseqid)
                        else:
                            # 此reads 已存在该属里面，不用append ，counts也不用 +1
                            pass

            total_counts += 1   ###文件的行数


    print("total counts is %s"%(total_counts))
    #############compute species abudance###################
    for k, v in blast_out_s.items():
        #s_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        s_abundance = "%.4f" % (float(v["count"]) / float(total_counts))
        #del blast_out_s[k]["qseqid"]
        if s_abundance not in blast_out_s[k].keys():
            blast_out_s[k].update({"s_abundance": s_abundance})
        else:
            pass

    ##################保存种id 并排序######################################
    tmp_tb = pd.DataFrame.from_dict(blast_out_s).transpose()
    tb= tmp_tb.sort_values(by="s_abundance", axis=0, ascending=False,inplace=False)


    ###删除seqid  保存blast_s 结果
    s_tb = tb[['taxid', 'lineage_id', 'lineage_name', 'count','s_abundance']]
    outblast_s = Path(workdir)/(f'{prefix}.s.csv')
    s_tb.to_csv(outblast_s,sep=',',index=False)
    del s_tb


    ### 保存序列top10 物种id  到json #####
    tb_top10 = tb.head(n=10)
    # print("11111111111111111")
    # print(tb_top10)
    tb_top10 = tb_top10.set_index("taxid",drop=False) ##taxid 作为后续的fa 命名前缀
    #tb_top10.to_csv("/mnt/home/huanggy/project/meta/result/align/top10.id",sep="\t",index=False)
    outblast_j = f'{workdir}/{prefix}.top10.json'
    #print(outblast_j,type(outblast_j))
    tb_top10.to_json(outblast_j)


    #############################保存属id ，排序############################
    for k, v in blast_out_g.items():
        #s_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        g_abundance = "%.4f" % (float(v["count"]) / float(total_counts))
        #del blast_out_s[k]["qseqid"]
        if g_abundance not in blast_out_g[k].keys():
            blast_out_g[k].update({"g_abundance": g_abundance})
        else:
            pass


    g_tmp_tb = pd.DataFrame.from_dict(blast_out_g).transpose()
    g_tb= g_tmp_tb.sort_values(by="g_abundance", axis=0, ascending=False,inplace=False)
    ###删除seqid  保存blast_g 结果#######
    g_tb = g_tb[['taxid', 'lineage_id', 'lineage_name', 'count','g_abundance']]
    outblast_g = Path(workdir)/(f'{prefix}.g.csv')
    g_tb.to_csv(outblast_g,sep=',',index=False)
    del g_tb


    #### 根据json 结果  seqid  从unmap 数据中提取taxid.fa

    #####gmap 比对####

    ###samtools  sort index  depth 统计深度#####


    ###可视化 覆盖度######


    #species_name = tb
    return  blast_out_s



#print(mydb)
# print(blast_out_s)
# print(end)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage")
    else:
        #blast_out_s = parser_blast4(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out")
        blast_out_s = parser_blast4(blast_outm8=sys.argv[1],workdir=sys.argv[2],prefix=sys.argv[3])
        end = datetime.datetime.now()
        print(end)


