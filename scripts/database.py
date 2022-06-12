#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


import pandas as pd
import pymysql


db = pymysql.connect(
    host="127.0.0.1",
    port=3306,
    user = "huanggy",
    password="654321",
    database = "taxdb",
    charset="utf8"
)
print(db)

# 创建表`lineage_db`
create_sql_lineage_db = """
    DROP table IF EXISTS lineage_db ;
    CREATE TABLE `lineage_db` (
      `id` int(16) NOT NULL AUTO_INCREMENT,
      `s_taxid` int(24) DEFAULT NULL,
      `names_lineage` varchar(100) DEFAULT NULL,
      `taxid_lineage` varchar(100) DEFAULT NULL,
    ) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
"""


def creat_lineage_db():
    sql = """
        CREATE TABLE lineage_db (
          id int(16) NOT NULL AUTO_INCREMENT,
          s_taxid int(24) DEFAULT NULL,
          names_lineage varchar(300) DEFAULT NULL,
          taxid_lineage varchar(100) DEFAULT NULL,
          primary key (id),
          constraint taxids foreign key (s_taxid) references  refseq_taxdb (id)
        )ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=utf8
        """
    cursor = db.cursor()

    cursor.execute("DROP TABLE IF EXISTS lineage_db")
    cursor.execute(sql)
    print("CREATE TABLE lineage_db   ok   ")
    db.close()

##foreign key (s_taxid) references refseq_taxdb(taxonomy_ID),

def insert_lineage_db():
    cursor = db.cursor()
    # cursor.execute("DROP TABLE IF EXISTS lineage_db")
    # cursor.execute(sql)
    # db.close()

    data = pd.read_table("/mnt/home/huanggy/project/meta/result/align/select_assembly_info_20210826_taxid_lineage_db",sep="\t", header=0)

    query = """insert into lineage_db (s_taxid, names_lineage, taxid_lineage) values (%s, %s, %s)"""
    for ind, row in data.iterrows():
        # print(ind,"\n",row)
        # id = ind
        s_taxid = row["speciess_taxid"]
        names_lineage = row["names_lineage"]
        taxid_lineage = row["taxid_lineage"]
        values = (s_taxid, names_lineage, taxid_lineage)
        # print(id)
        cursor.execute(query, values)

    db.commit()
    cursor.close()
    db.close()

def creat_refseq_taxdb():
    sql = """
        CREATE TABLE refseq_taxdb (
          id int(16) NOT NULL AUTO_INCREMENT,
          taxonomy_ID int(24) DEFAULT NULL,
          species_name varchar(100) DEFAULT NULL,
          accession_version varchar(100) DEFAULT NULL,
          refseq_description varchar (100) DEFAULT NULL,
          refseq_status varchar (16) DEFAULT NULL,
          genome_length int (32) DEFAULT NULL,
          primary key (id)
        )ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=utf8
        """

    cursor = db.cursor()

    cursor.execute("DROP TABLE IF EXISTS refseq_taxdb")
    cursor.execute(sql)
    db.close()



def insert_refseq_taxdb():
    cursor = db.cursor()
    data = pd.read_table("/mnt/data/NCBI/taxonomy//RefSeq_bac_fun.txt",sep="\t", header=0)

    query = """insert into refseq_taxdb (taxonomy_ID, species_name, accession_version,refseq_description,refseq_status,genome_length) values (%s,%s,%s,%s,%s,%s)"""
    for ind, row in data.iterrows():
        # print(ind,"\n",row)
        # id = ind
        taxonomy_ID = row["taxonomy ID"]
        species_name = row["species name"]
        accession_version = row["accession.version"]
        refseq_description = row["refseq description"]
        refseq_status = row["refseq status"]
        genome_length = row["length"]
        values = (taxonomy_ID, species_name, accession_version,refseq_description,refseq_status,genome_length)
        # print(id)
        cursor.execute(query, values)

    db.commit()
    cursor.close()
    db.close()


#creat_refseq_taxdb()
#insert_refseq_taxdb()


creat_lineage_db()
#insert_lineage_db()




