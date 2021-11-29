#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pandas as pd
import pymysql
# temp = {'a':[1,2,3,4,5,6],'b':['Alice','Bob','Cindy','Eric','Helen','Grace '],'c':[90,89,99,78,97,93],'d':[89,94,80,94,94,90]}
# data = pd.DataFrame(temp)
# print(data)


# for i in range(len(data)):
#     a = data['a'].iloc[i]
#     b = data['b'].iloc[i]
#     c = data['c'].iloc[i]
#     d = data['d'].iloc[i]
#     print(a,b,c,d)


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



sql = """
    CREATE TABLE lineage_db (
      id int(16) NOT NULL AUTO_INCREMENT,
      s_taxid int(24) DEFAULT NULL,
      names_lineage varchar(300) DEFAULT NULL,
      taxid_lineage varchar(100) DEFAULT NULL,
      primary key (id)
    )ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=utf8
    """




cursor = db.cursor()
# cursor.execute("DROP TABLE IF EXISTS lineage_db")
# cursor.execute(sql)
# db.close()

data = pd.read_table("/mnt/home/huanggy/project/meta/result/align/select_assembly_info_20210826_taxid_lineage_db",sep="\t",header=0)


query = """insert into lineage_db (s_taxid, names_lineage, taxid_lineage) values (%s, %s, %s)"""
for ind ,row in data.iterrows():
    #print(ind,"\n",row)
    #id = ind
    s_taxid = row["speciess_taxid"]
    names_lineage = row["names_lineage"]
    taxid_lineage = row["taxid_lineage"]
    values = (s_taxid, names_lineage, taxid_lineage)
    #print(id)
    cursor.execute(query,values)

db.commit()
cursor.close()
db.close()
