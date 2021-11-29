#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

class Record(object):
    """
    One line information in accession2taxid db
    """

    def __init__(self,line):
        self.line = line
        info = self.line.split("\t")
        if len (info) == 6:
            #####Accession2taxid_db######
            self.accession = info[2]
            self.taxid = info[0]
            self.scientific_name = info[1]

        if len (info) == 3:
            ######Lineage_db#####
            self.taxid = info[0]
            self.scientific_lineage = info[1]
            self.taxid_lineage = info [2]


class Accession2taxid_db(object):
    def __init__(self, db):
        self.reader = open(db, 'r')
        self.line = self.reader.readline().strip()
        while self.line.startswith('#'):
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

    def quire_accessino(self,accessino):
        """
        通过accession 查找taxid
        :return:  taxid
        """
        if self.record.accession == accessino:
            return self.record.taxid

    def quire_taxid(self,taxid):
        """
        通过taxid 查找 accession
        :param taxid:
        :return:
        """
        if self.record.taxid == taxid:
            return self.record.accession

    def quire_taxid_lineage(self):
        """
        通过taxid 查找 accession
        :param taxid:
        :return:
        """
        if self.record.taxid == self.foreign_key:
            ###同过本记录的taxid 所对应的外键 去找外键对应Lineage_db  filed info
            return self.foreign_key.scientific_lineage,self.foreign_key.taxid_lineage


class Lineage_db(object):
    def __init__(self, db):
        self.reader = open(db, 'r')
        self.line = self.reader.readline().strip()
        while self.line.startswith('#'):
            continue
        self.record = Record(self.line)

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



accession2taxid_db_iteror = Accession2taxid_db(db="/mnt/data/NCBI/taxonomy//RefSeq_bac_fun.txt.back")
lineage_db_iteror = Lineage_db(db="/mnt/home/huanggy/project/meta/result/align/select_assembly_info_20210826_taxid_lineage_db")

accession_id = "NZ_WLVG01000001.1"

re_taxid = accession2taxid_db_iteror.quire_accessino(accession_id)
print(re_taxid)
for lineage_itor in lineage_db_iteror:
    for accession_itor in accession2taxid_db_iteror:

        if lineage_itor.taxid == accession_itor.taxid:
                print("_____")
                print(accession_id,lineage_itor.taxid,lineage_itor.scientific_lineage,lineage_itor.taxid_lineage)