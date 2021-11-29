#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

HOST = 'localhost'
PORT = 3306
USERNAME = 'huanggy'
PASSWORD = '654321'
DB = 'taxdb'

# dialect + driver://username:passwor@host:port/database
DB_URI = f'mysql+pymysql://{USERNAME}:{PASSWORD}@{HOST}:{PORT}/{DB}'



from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, MetaData,Column, Integer, String
from sqlalchemy.orm import sessionmaker
import pandas as pd


engine = create_engine(DB_URI,encoding='utf-8')
Base = declarative_base(engine)  # SQLORM基类
meta = MetaData()
session = sessionmaker(engine)()  # 构建session对象


#
class Lineage_db(Base):
    """
    https://www.cnblogs.com/goldsunshine/p/10124859.html
    """
    __tablename__ = 'lineage_dbs'  # 表名
    id = Column(Integer, primary_key=True, autoincrement=True)
    taxid = Column(Integer)
    names_lineage = Column(String(600))
    taxid_lineage = Column(String(400))
    name = Column(String(320))
    f_names_lineage = Column(String(800))



class Refseq_taxdb(Base):
    __tablename__ = "refseq_taxdb"
    id = Column(Integer, primary_key=True, autoincrement=True)
    taxonomy_ID = Column(Integer)
    species_name = Column(String(128))
    accession_version = Column(String(128))
    refseq_description = Column(String(128))
    refseq_status = Column(String(12))
    genome_length = Column(String(32))
#Base.metadata.create_all()  # 将模型映射到数据库中
#在数据库中创建表
#Lineage_db.metadata.create_all(bind=engine)
#Refseq_taxdb.metadata.create_all(bind=engine)

def insert_Lineage_db():
    #data = pd.read_table("/mnt/home/huanggy/project/meta/result/align/select_assembly_info_20210826_taxid_lineage_db",sep="\t", header=0)
    data = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/refseq_taxid_lineage_db",sep="\t", header=0)
    for ind, row in data.iterrows():
        lineage = Lineage_db(
            taxid=row["taxid"],
            names_lineage=row["names_lineage"],
            taxid_lineage=row["taxid_lineage"],
            name = row["name"],
            f_names_lineage = row["format_names_lineage"]
        )
        session.add(lineage)  # 添加到session
        session.commit()  # 提交到数据库

    # student = LineageDb( s_taxid=280145, names_lineage='cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter;Acinetobacter colistiniresistens',taxid_lineage = '131567;2;1224;1236;72274;468;469;280145')  # 创建一个student对象
    # session.add(student)  # 添加到session
    # session.commit()  # 提交到数据库


def insert_Refseq_taxdb():
    #data = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun.txt", sep="\t", header=0)
    data = pd.read_table("/mnt/data/NCBI/taxonomy/refseq_taxonomy/RefSeq_bac_fun_viral.txt", sep="\t", header=0)
    for ind ,row in data.iterrows():
        refseq_taxdb = Refseq_taxdb(
            taxonomy_ID=row["taxonomy_ID"],
            species_name = row["species_name"],
            accession_version = row["accession_version"],
            refseq_description = row["refseq_description"],
            refseq_status = row["refseq_status"],
            genome_length = row["length"]
        )
        session.add(refseq_taxdb)  # 添加到session
        session.commit()  # 提交到数据库



#insert_Lineage_db()
#insert_Refseq_taxdb()

# item = session.query(LineageDb.s_taxid).first()
# print(item)


# item2 = session.query(Refseq_taxdb.taxonomy_ID).first()
# print(item2)
#
#
# #inputid = input_accessionid = input("NZ_JODL01000187.1")
# inputid = "NZ_JODL01000187.1"
#
# refobj3 = session.query(Refseq_taxdb.taxonomy_ID).filter(Refseq_taxdb.accession_version == inputid).all()
#
# print('333333',refobj3[0][0],type(refobj3[0]))
# print(refobj3[0][0])



#result = session.query(Refseq_taxdb).join(LineageDb).filter(Family.member>6)



def accession_id2taxon_id(accession_id):

    if session.query(Refseq_taxdb).filter(Refseq_taxdb.accession_version == accession_id).all():
        quire_taxonid = session.query(Refseq_taxdb.taxonomy_ID).filter(Refseq_taxdb.accession_version == accession_id).all()
        return  quire_taxonid[0][0]

    else:
        return "不存在数据库"


def taxon_id2lineage(accession_id):
    taxid = accession_id2taxon_id(accession_id)
    #print('---->',taxid)
    res = session.query(Lineage_db.f_names_lineage,Lineage_db.taxid_lineage).filter(Lineage_db.taxid == taxid).all()
    #print(res)
    if res:
        res = res[0]
    return res


if __name__ == "__main__":
    print(taxon_id2lineage("NZ_FNVE01000001.1"))
