#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang
import re
import datetime
import pandas as pd
class Record(object):
    """
    One line information in accession2taxid db
    """

    def __init__(self,line):
        self.line = line
        info = self.line.split("\t")
        #####Accession2taxid_db######
        self.accession_version = info[2]
        self.taxid = info[0]
        self.scientific_name = info[1]
        self.taxid_lineage = info[9]
        self.format_names_lineage = info[9]





class Taxid_lineage_db(object):
    def __init__(self, db):
        self.reader = open(db, 'r')
        self.line = self.reader.readline().strip()
        while self.line.startswith('taxid'):
            continue
        self.record = Record(self.line)
        #self.foreign_key =Lineage_db(db="").record.taxid   #####通过taxid 建立外键

    def __iter__(self):
        return self

    def __next__(self):
        self.line = self.reader.readline().strip()
        if self.line != "":
            self.record = Record(self.line)
            return self.record
        else:
            self.reader.close()
            raise StopIteration()

    def reader_close(self):
        self.reader.close()




def dbasdict(db):
    db_dict = {

    }
    with open(db,'r') as fh:
        for line in fh:
            if not line:
                break
            if line.startswith("taxid"):
                continue
            taxid,species_name,accession_version,refseq_description,refseq_status,length,names_lineage,taxid_lineage,name,format_names_lineage = line.strip().split("\t")
            if accession_version not in db_dict:
                #db_dict[accession_version] =
                #print('----->',accession_version)
                db_dict.update({accession_version:{"taxid":taxid,"taxid_lineage":taxid_lineage,"format_names_lineage":format_names_lineage}})
            else:
                pass

    return db_dict




def iterator_taxid_lineage_db(taxid_lineage_db):
    "/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db"
    with open (taxid_lineage_db,'r') as  db_path_f:
        #next(db_path_f)
        for line in db_path_f:
            if line.startswith("taxid"):
                pass
            else:
                if len(line.strip().split("\t")) < 10:
                    print("数据格式？？？？？？？？？？？？？",line)
                else:
                    taxid, species_name, accession_version, refseq_description, refseq_status, length, names_lineage, taxid_lineage, name, format_names_lineage = line.strip().split("\t")
                    #yield (taxid,accession_version,taxid_lineage,format_names_lineage)
                    yield {
                        accession_version: {"taxid":taxid,
                                            "taxid_lineage":taxid_lineage,
                                            "format_names_lineage":format_names_lineage}
                    }


#iterator_taxid_lineage_db("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")
# for d in iterator_taxid_lineage_db("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db"):
#     #print ("=======>",d)
#     a = "NC_016845.1"
#     if a in d.keys():
#         print(d[a]["taxid"],"\n",d[a]["format_names_lineage"])


def get_accession_id2species_name(sseqid,taxid_lineage_db):
    for d in iterator_taxid_lineage_db(taxid_lineage_db):
        if sseqid in d.keys():
            return  (d[sseqid]["taxid"],d[sseqid]["format_names_lineage"],d[sseqid]["taxid_lineage"])


#get_accession_id2species_name("NC_016845.1","/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")

def parser_blast2(blast_outm8,taxid_lineage_db):
    blast_out_s = {}
    blast_out_g = {}
    with open(blast_outm8, "r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            if not line :
                break
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")

            ###由于病毒的fa 序列 是采用 kraken2 标准库的，故sseid 以kraken 开头 kraken:taxid|1922926|NC_033306.1  在refseq 和世系数据库不存在，故另行处理
            if sseqid.startswith("kraken"):
                sseqid = sseqid.split("|")[2]
            if sseqid:

                # 'NoneType' object is not iterable
                res = get_accession_id2species_name(sseqid, taxid_lineage_db)
                if res:
                    #print("=========>", qseqid,sseqid)
                    #taxid, format_names_lineage,taxid_lineage= get_accession_id2species_name(sseqid,taxid_lineage_db)
                    taxid, format_names_lineage, taxid_lineage =  res
                    if format_names_lineage :
                        #### parser species#####
                        s_name = re.match('.*\|(s__.*)\|.*', format_names_lineage).groups()[0]
                        if s_name not in blast_out_s:
                            blast_out_s.update({s_name: {"lineage_id": taxid_lineage, "lineage_name": format_names_lineage, "count": 1,
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

    ### compute species abudance
    for k, v in blast_out_s.items():
        s_abundance = "{:.2%}".format(float(v["count"]) / float(total_counts))
        del blast_out_s[k]["qseqid"]
        if s_abundance not in blast_out_s[k].keys():
            blast_out_s[k].update({"s_abundance": s_abundance})
        else:
            pass
    return blast_out_s



# import pandas as pd
# pd.merge

# blast_out_s = parser_blast2(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out",taxid_lineage_db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")
# print(blast_out_s)

#taxid_lineage_db_iteror = Taxid_lineage_db(db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")
def parser_blast3(blast_outm8):
    blast_out_s = {}
    blast_out_g = {}
    with open(blast_outm8, "r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            if not line :
                break
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")
            if sseqid.startswith("kraken"):
                sseqid = sseqid.split("|")[2]
            if sseqid:
                for itor in taxid_lineage_db_iteror:
                    if sseqid == itor.accession_version:
                        s_name = re.match('.*\|(s__.*)\|.*', itor.format_names_lineage).groups()[0]
                        lineage_id = itor.taxid_lineage
                        lineage_name = itor.format_names_lineage
                        if s_name not in blast_out_s:
                            blast_out_s.update({s_name: {"lineage_id": lineage_id, "lineage_name": lineage_name, "count": 1,"qseqid": [qseqid]}})
                        else:
                            counts = blast_out_s[s_name]["count"] + 1
                            # blast_out.update({s_level}.update({"count":counts}))
                            blast_out_s[s_name].update({"count": counts})

                            ##同种 不同的reads append
                            if qseqid not in blast_out_s[s_name]["qseqid"]:
                                blast_out_s[s_name]["qseqid"].append(qseqid)
                        total_counts += 1

    return  blast_out_s

# blast_out_s = parser_blast3(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out")
# print(blast_out_s)



def parser_blast4(blast_outm8):
    start = datetime.datetime.now()
    print(start)
    mydb = dbasdict(db="/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.lineage_db")

    blast_out_s = {}
    blast_out_g = {}
    with open(blast_outm8, "r") as blast_outm8_f:
        total_counts = 0
        evalue_cutoff = ''
        bitscore_cutoff = ''
        for line in blast_outm8_f:
            if not line :
                break
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")
            if sseqid.startswith("kraken"):
                sseqid = sseqid.split("|")[2]
            if sseqid in mydb:
                taxid = mydb[sseqid]["taxid"]
                lineage_id = mydb[sseqid]["taxid_lineage"]
                #lineage_name = mydb[sseqid]["format_names_lineage"]
                format_names_lineage = mydb[sseqid]["format_names_lineage"]
                s_name = re.match('.*\|(s__.*)\|.*', format_names_lineage).groups()[0]
                if s_name not in blast_out_s:
                    blast_out_s.update({s_name: {"lineage_id": lineage_id, "lineage_name": format_names_lineage, "count": 1,
                                                 "qseqid": [qseqid]}})
                else:
                    counts = blast_out_s[s_name]["count"] + 1
                    # blast_out.update({s_level}.update({"count":counts}))
                    blast_out_s[s_name].update({"count": counts})

                    ##同种 不同的reads append
                    if qseqid not in blast_out_s[s_name]["qseqid"]:
                        blast_out_s[s_name]["qseqid"].append(qseqid)
                total_counts += 1


                #
                # for itor in taxid_lineage_db_iteror:
                #     if sseqid == itor.accession_version:
                #         s_name = re.match('.*\|(s__.*)\|.*', itor.format_names_lineage).groups()[0]
                #         lineage_id = itor.taxid_lineage
                #         lineage_name = itor.format_names_lineage
                #         if s_name not in blast_out_s:
                #             blast_out_s.update({s_name: {"lineage_id": lineage_id, "lineage_name": lineage_name, "count": 1,"qseqid": [qseqid]}})
                #         else:
                #             counts = blast_out_s[s_name]["count"] + 1
                #             # blast_out.update({s_level}.update({"count":counts}))
                #             blast_out_s[s_name].update({"count": counts})
                #
                #             ##同种 不同的reads append
                #             if qseqid not in blast_out_s[s_name]["qseqid"]:
                #                 blast_out_s[s_name]["qseqid"].append(qseqid)
                #         total_counts += 1


    ### compute species abudance
    for k, v in blast_out_s.items():
        #s_abundance = "{:.2}".format(float(v["count"]) / float(total_counts))
        s_abundance = "%.4f"%(float(v["count"]) / float(total_counts))
        #del blast_out_s[k]["qseqid"]
        if s_abundance not in blast_out_s[k].keys():
            blast_out_s[k].update({"s_abundance": s_abundance})
        else:
            pass

    ####top 10 物种####
    tmp_tb = pd.DataFrame.from_dict(blast_out_s).transpose()
    tb= tmp_tb.sort_values(by="s_abundance", axis=0, ascending=False,inplace=False)
    del tmp_tb
    tb.to_csv("/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.parser",sep="\t")

    ### 保存序列top10 物种id  到json #####


    return  blast_out_s


blast_out_s = parser_blast4(blast_outm8="/mnt/home/huanggy/project/meta/result/align/my_gsnap20210916.out")
end = datetime.datetime.now()

#print(mydb)
#print(blast_out_s)
print(end)



















