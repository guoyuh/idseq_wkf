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

def count_clena_reads(clean_fastq):
    cnt = 0
    with open(clean_fastq,'r') as fh:
        for line in fh:
            if line.startswith("@"):
                cnt += 1
            else:
                pass

    return  cnt




def parser_blast(blast_outm8,workdir,prefix,topn = 10):
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    #计算去人源，去低质量后的总reads ，并计算RPM 每百万序列数（reads per million，RPM）#################################
    clean_fq = Path(workdir)/(f'{prefix}/filter/{prefix}.f.fastq')
    if os.path.exists(clean_fq):
        in_counts = count_clena_reads(clean_fastq = clean_fq)
    else:
        in_counts = 1000000

    ####################################################################################################################
    mydb = dbasdict(db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral_protozoa.lineage_db")
    #total_counts_bac, total_counts_fungi, total_counts_viral, total_counts_protozoa, mydb = parser_blast4_bac_fungi_viral_protozoa_counts(blast_outm8=blast_outm8)
    total_counts_bac = 0
    total_counts_fungi = 0
    total_counts_viral = 0
    total_counts_protozoa = 0
    blast_out_s = {}
    blast_out_g = {}
    total_counts_valid = []

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


                if re.search("k__Bacteria",format_names_lineage):
                    total_counts_bac +=  1
                elif re.search("k__Eukaryota",format_names_lineage):
                    total_counts_fungi += 1

                elif re.search("k__Viruses",format_names_lineage):
                    total_counts_viral +=1
                elif re.search("k__protozoa",format_names_lineage):
                    total_counts_protozoa += 1
                else:
                    print("其它类别,maybe human")
                ###################
                if s_name not in blast_out_s:
                    ##初始化新的物种信息
                    blast_out_s.update({s_name: {"taxid":taxid,"lineage_id": lineage_id, "lineage_name": format_names_lineage, "count": 1,"qseqid": [qseqid]}})

                else:
                    ##同种 不同的reads append(在基因组不同的位置)
                    #if qseqid not in blast_out_s[s_name]["qseqid"] and len(blast_out_s[s_name]["qseqid"])<3:
                    if qseqid not in blast_out_s[s_name]["qseqid"]:
                        counts = blast_out_s[s_name]["count"] + 1
                        blast_out_s[s_name].update({"count": counts})
                        blast_out_s[s_name]["qseqid"].append(qseqid)
                    else:
                        ## 相同的reads 比对不同的亚种 （或同种同源区？？？？）
                        pass


                #### parser gensu#####
                for ind, i in enumerate(format_names_lineage.split("|")):
                    if i.startswith("g__"):
                        g_name = i
                        g_taxid = lineage_id.split(";")[ind - len(format_names_lineage.split("|"))]
                    # else:
                    #     g_name = re.match('.*\|(g__.*)\|.*', format_names_lineage).groups()[0]
                    #     g_taxid = taxid
                    #my_g_name = re.match('.*\|(g__.*)\|.*', format_names_lineage).groups()[0]
                        if g_name not in blast_out_g:
                            blast_out_g.update({g_name: {"taxid":g_taxid,"lineage_id": lineage_id, "lineage_name": format_names_lineage, "count": 1, "qseqid": [qseqid]}})
                        else:
                            ##不同seqid 但是在同一个属  ,此时该属 +1 条hits
                            # blast_out.update({s_level}.update({"count":counts}))

                            # 同属 total count +1
                            #    同属不同种 ，该属水平 reads + 1
                            #    同属同种，该属水平 reads 不加1
                            if qseqid not in blast_out_g[g_name]["qseqid"]:
                                counts = blast_out_g[g_name]["count"] + 1
                                blast_out_g[g_name].update({"count": counts})
                                blast_out_g[g_name]["qseqid"].append(qseqid)
                            else:
                                ##同一seqid 比对到不同的种水平 ，但是是同一个属 ，此时属于该属水平的 不能加1 ，属于重复reads
                                # 此reads 已存在该属里面，不用append ，counts也不用 +1
                                pass

            total_counts_valid.append(qseqid)

    total_counts = len(list(set(total_counts_valid)))  ##total counts 是通过过滤的且不重复的reads 数
    print("total valid counts is %s"%(total_counts))
    #############compute species abudance###################
    for k, v in blast_out_s.items():
        ######计算相对丰度#####
        if re.search("k__Bacteria", v["lineage_name"]):
            s_abundance = "%.4f" % (float(v["count"]) / float(total_counts_bac))
        elif re.search("k__Eukaryota", v["lineage_name"]):
            s_abundance = "%.4f" % (float(v["count"]) / float(total_counts_fungi))

        elif re.search("k__Viruses", v["lineage_name"]):
            s_abundance = "%.4f" % (float(v["count"]) / float(total_counts_viral))
        elif re.search("k__protozoa", v["lineage_name"]):
            s_abundance = "%.4f" % (float(v["count"]) / float(total_counts_protozoa))
        else:
            print("不属于任何一个k__Bacteria，k__Eukaryota，k__Viruses，k__protozoa大类")
            s_abundance = 0

        #####计算相对丰度####
        #s_abundance = "%.4f" % (float(v["count"]) / float(total_counts))
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

    s_tb['RPM'] = "%.2f" % (float(tb['count']/in_counts *1.0e6))
    outblast_s = Path(workdir)/(f'{prefix}.s.csv')
    s_tb.to_csv(outblast_s,sep=',',index=False)
    del s_tb


    ### 保存序列top10 物种id  到json #####
    tb_top10 = tb.head(n=topn)
    # print("11111111111111111")
    # print(tb_top10)
    tb_top10 = tb_top10.set_index("taxid",drop=False) ##taxid 作为后续的fa 命名前缀
    #tb_top10.to_csv("/mnt/home/huanggy/project/meta/result/align/top10.id",sep="\t",index=False)
    outblast_j = f'{workdir}/{prefix}.top10.json'
    #print(outblast_j,type(outblast_j))
    tb_top10.to_json(outblast_j)


    #############################保存属id ，排序############################
    for k, v in blast_out_g.items():
        ######计算相对丰度#####
        if re.search("k__Bacteria", v["lineage_name"]):
            g_abundance = "%.4f" % (float(v["count"]) / float(total_counts_bac))
        elif re.search("k__Eukaryota", v["lineage_name"]):
            g_abundance = "%.4f" % (float(v["count"]) / float(total_counts_fungi))

        elif re.search("k__Viruses", v["lineage_name"]):
            g_abundance = "%.4f" % (float(v["count"]) / float(total_counts_viral))
        elif re.search("k__protozoa", v["lineage_name"]):
            g_abundance = "%.4f" % (float(v["count"]) / float(total_counts_protozoa))
        else:
            print("不属于任何一个k__Bacteria，k__Eukaryota，k__Viruses，k__protozoa大类")
            g_abundance = 0

        #######计算绝对丰度######
        #g_abundance = "%.4f" % (float(v["count"]) / float(total_counts))
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


def parser_blast4_bac_fungi_viral_protozoa_counts(blast_outm8):
    total_counts_bac = 0
    total_counts_fungi = 0
    total_counts_viral = 0
    total_counts_protozoa = 0
    mydb = dbasdict(db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral_protozoa.lineage_db")
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

            #subseqid is seq accession_num  如： >kraken:taxid|2698686|NZ_CP048031.1
            if subseqid.startswith("kraken"):
                subseqid = subseqid.split("|")[2]
            if subseqid in mydb:
                taxid = mydb[subseqid]["taxid"]
                lineage_id = mydb[subseqid]["taxid_lineage"]
                #lineage_name = mydb[sseqid]["format_names_lineage"]
                format_names_lineage = mydb[subseqid]["format_names_lineage"]
                s_name = re.match('.*\|(s__.*)\|.*', format_names_lineage).groups()[0]

                if re.search("k__Bacteria",format_names_lineage):
                    total_counts_bac +=  1
                elif re.search("k__Eukaryota",format_names_lineage):
                    total_counts_fungi += 1

                elif re.search("k__Viruses",format_names_lineage):
                    total_counts_viral +=1
                elif re.search("k__protozoa",format_names_lineage):
                    total_counts_protozoa += 1
                else:
                    print("其它类别,maybe human")

        return total_counts_bac,total_counts_fungi ,total_counts_viral,total_counts_protozoa,mydb

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage")
    else:
        start = datetime.datetime.now()
        print(start)
        #blast_out_s = parser_blast4(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out")
        blast_out_s = parser_blast(blast_outm8=sys.argv[1],workdir=sys.argv[2],prefix=sys.argv[3])
        end = datetime.datetime.now()
        print(end)


