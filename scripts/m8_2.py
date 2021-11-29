#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


from  database_2 import taxon_id2lineage
import re


###
##report   taxid_lineage_db  鉴定到种水平
##db nt  鉴定到 亚种 水平
# taxid :{"taxid":taxid,
#           "level":lineage,
#           "count";int,
#          "family_cnt":int,
#           "genus cnt ":int,
#             "species_cnt":int
#             "secies_name",
#             "genus_name"
# }

# {taxid :{"taxid":taxid,
#           "level":species,
#           "secies_name",: str
#           "count";int,}
#}


#blast_outm8 = '/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out'


#{('Klebsiella pneumoniae', '574'):1}
# blast_out = {('Klebsiella pneumoniae', '573'):0}
# blast_out.update({('Klebsiella pneumoniae', '574'):1})
# print(blast_out)

# with open(blast_outm8,"r") as blast_outm8_f:
#     total_counts = 0
#     for line in blast_outm8_f:
#         qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.split("\t")
#         s_level = get_accession_id2species_name(sseqid,taxid_lineage_db = "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db")
#         if s_level not in blast_out_s:
#             # qseqid_list = []
#             # qseqid_list.append(qseqid)
#             blast_out_s.update({s_level:{"species_id":s_level[1],"species_name":s_level[0],"count":1,"qseqid":[qseqid]}})
#         else:
#             counts = blast_out_s[s_level]["count"] + 1
#             #blast_out.update({s_level}.update({"count":counts}))
#             blast_out_s[s_level].update({"count":counts})
#             ###同一seqid 比对到不同的种水平 ，但是是同一个属
#             # if qseqid not in blast_out_s[s_level]["qseqid"]:
#             #     blast_out_s[s_level]["qseqid"].append(qseqid)
#
#         total_counts += 1
#
# print(blast_out_s)
# print(total_counts)
#
# for k , v in blast_out_s.items():
#     s_abundance = "{:.2%}".format(float(v["count"])/float(total_counts))
#     if s_abundance not in blast_out_s[k].keys():
#         blast_out_s[k].update({"s_abundance":s_abundance})
#     else:
#         pass
# print(blast_out_s)

def parser_blast_s():
    with open(blast_outm8, "r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split(
                "\t")
            lineage = taxon_id2lineage(sseqid)
            # print("------->",s_level)
            # ('k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella pneumoniae|t__Klebsiella pneumoniae subsp. pneumoniae HS11286', '131567;2;1224;1236;91347;543;570;573;72407;1125630')
            if lineage:
                s_name = re.match('.*\|(s__.*)\|.*', lineage[0]).groups()[0]
                if s_name not in blast_out_s:
                    blast_out_s.update({s_name: {"lineage_id": lineage[1], "lineage_name": lineage[0], "count": 1,
                                                 "qseqid": [qseqid]}})
                    # print(s_name)
                else:
                    counts = blast_out_s[s_name]["count"] + 1
                    # blast_out.update({s_level}.update({"count":counts}))
                    blast_out_s[s_name].update({"count": counts})

                    ##同种 不同的reads append
                    if qseqid not in blast_out_s[s_name]["qseqid"]:
                        blast_out_s[s_name]["qseqid"].append(qseqid)
                total_counts += 1

    for k, v in blast_out_s.items():
        s_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        if s_abundance not in blast_out_s[k].keys():
            blast_out_s[k].update({"s_abundance": s_abundance})
        else:
            pass

def parser_blast_g():
    with open(blast_outm8,"r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.split("\t")
            lineage = taxon_id2lineage(sseqid)
            if lineage:
                g_name = re.match('.*\|(g__.*)\|.*',lineage[0]).groups()[0]

                if g_name not in blast_out_g:
                    blast_out_g.update({g_name:{"lineage_id":lineage[1],"lineage_name":lineage[0],"count":1,"qseqid":[qseqid]}})
                else:
                    ##不同seqid 但是在同一个属  ,此时该属 +1 条hits
                    ##同一seqid 比对到不同的种水平 ，但是是同一个属 ，此时属于该属水平的 不能加1 ，属于重复reads
                    #blast_out.update({s_level}.update({"count":counts}))

                    #同属 total count +1
                    #    同属不同种 ，该属水平 reads + 1
                    #    同属同种，该属水平 reads 不加1
                    if qseqid not in blast_out_g[g_name]["qseqid"]:
                        counts = blast_out_g[g_name]["count"] + 1
                        blast_out_g[g_name].update({"count":counts})
                        blast_out_g[g_name]["qseqid"].append(qseqid)
                    else:
                        #此reads 已存在该属里面，不用append ，counts也不用 +1
                        pass

                total_counts += 1




def parser_blast(blast_outm8):
    blast_out_s = {}
    blast_out_g = {}
    with open(blast_outm8, "r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split(
                "\t")
            lineage = taxon_id2lineage(sseqid)
            # print("------->",s_level)
            # ('k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella pneumoniae|t__Klebsiella pneumoniae subsp. pneumoniae HS11286', '131567;2;1224;1236;91347;543;570;573;72407;1125630')
            if lineage:
                #### parser species#####
                s_name = re.match('.*\|(s__.*)\|.*', lineage[0]).groups()[0]
                if s_name not in blast_out_s:
                    blast_out_s.update({s_name: {"lineage_id": lineage[1], "lineage_name": lineage[0], "count": 1,
                                                 "qseqid": [qseqid]}})
                    # print(s_name)
                else:
                    counts = blast_out_s[s_name]["count"] + 1
                    # blast_out.update({s_level}.update({"count":counts}))
                    blast_out_s[s_name].update({"count": counts})

                    ##同种 不同的reads append
                    if qseqid not in blast_out_s[s_name]["qseqid"]:
                        blast_out_s[s_name]["qseqid"].append(qseqid)

                #### parser gensu#####
                g_name = re.match('.*\|(g__.*)\|.*', lineage[0]).groups()[0]
                if g_name not in blast_out_g:
                    blast_out_g.update({g_name: {"lineage_id": lineage[1], "lineage_name": lineage[0], "count": 1,"qseqid": [qseqid]}})
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

                total_counts += 1

    ### compute species abudance
    for k, v in blast_out_s.items():
        s_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        del blast_out_s[k]["qseqid"]
        if s_abundance not in blast_out_s[k].keys():
            blast_out_s[k].update({"s_abundance": s_abundance})
        else:
            pass

    ###computer genus abudance
    for k, v in blast_out_g.items():
        g_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        del blast_out_g[k]["qseqid"]
        if g_abundance not in blast_out_g[k].keys():
            blast_out_g[k].update({"g_abundance": g_abundance})
        else:
            pass

    return blast_out_s,blast_out_g

if __name__ == "__main__":
    blast_out_s, blast_out_g = parser_blast(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out")
    with open("/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.parser_s.out","w") as out1,open("/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.parser_g.out","w") as out2 :
        for k, v in blast_out_s.items():
            out1.write(k + "\t" + str(v) + "\n" )

        for k,v in blast_out_g.items():
            out2.write(k + "\t" + str(v) + "\n")
    #print(blast_out_s,"\n",blast_out_g)









