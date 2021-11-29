#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

from csv import DictReader, DictWriter
from typing import Any, Dict, Iterable, Iterator, Sequence, Text, Tuple
import re


MAX_EVALUE_THRESHOLD = 1
class _TypedDictTSVReader(Iterator[Dict[str, Any]]):
    """
    Similar to DictReader but instead of a sequence of fieldnames it takes a sequence of
    tuples, the first element being the field name and the second being the field's type.
    After reading a row, this class will convert each element to it's associated type.
    """
    def __init__(self, f: Iterable[Text], schema: Sequence[Tuple[str, type]]) -> None:
        fieldnames = [field for field, _ in schema]
        self._types = {field: _type for field, _type in schema}
        # This is a class member instead of using  inheritance so we can return Dict[str, Any] instead of Dict[str, str]
        # Since str is a subtype of Any, if we add the Dict[str, Any] to our __next__ method signatures python will
        # use the more restrictive str type in place of Any. This makes it impossible to use this class as a
        # Iterator[Dict[str, Any]] as long as we inherit from DictReader, so it is made a class member instead.
        self._reader = DictReader(f, fieldnames=fieldnames, delimiter="\t")

    def __next__(self):
        row = next(self._reader)
        assert len(row) <= len(self._types), f"row {row} contains fields not in schema {self._types}"
        for key, value in row.items():
            if value:
                row[key] = self._types[key](value)
        return row

class _TypedDictTSVWriter(DictWriter):
    """
    This is just a convenience class so you don't need to pull out the field names
    """
    def __init__(self, f: Any, schema: Sequence[Tuple[str, type]]) -> None:
        fieldnames = [field for field, _ in schema]
        super().__init__(f, fieldnames, delimiter="\t")


class _BlastnOutput6ReaderBase(_TypedDictTSVReader):
    """
    This class is a bit of an oddball due to some compatibility concerns. In addition
    to the normal parsing stuff it also:
    1. Ignores comments (lines starting with '#') in tsv files, rapsearch2 adds them
    2. Supports filtering rows that we consider invalid
    """
    def __init__(self, f: Iterable[Text], schema: Sequence[Tuple[str, type]], filter_invalid: bool = False, min_alignment_length: int = 0):
        self._filter_invalid = filter_invalid
        self._min_alignment_length = min_alignment_length

        # The output of rapsearch2 contains comments that start with '#', these should be skipped
        filtered_stream = (line for line in f if not line.startswith("#"))
        super().__init__(filtered_stream, schema,)

    def __next__(self):
        if not self._filter_invalid:
            return super().__next__()
        row = super().__next__()
        while row and not self._row_is_valid(row):
            row = super().__next__()
        return row

    def _row_is_valid(self, row) -> bool:
        # GSNAP outputs bogus alignments (non-positive length /
        # impossible percent identity / NaN e-value) sometimes,
        # and usually they are not the only assignment, so rather than
        # killing the job, we just skip them. If we don't filter these
        # out here, they will override the good data when computing min(
        # evalue), pollute averages computed in the json, and cause the
        # webapp loader to crash as the Rails JSON parser cannot handle
        # NaNs. Test if e_value != e_value to test if e_value is NaN
        # because NaN != NaN.
        # *** E-value Filter ***
        # Alignments with e-value > 1 are low-quality and associated with false-positives in
        # all alignments steps (NT and NR). When the e-value is greater than 1, ignore the
        # alignment
        ###
        return all([
            row["length"] >= self._min_alignment_length,
            -0.25 < row["pident"] < 100.25,
            row["evalue"] == row["evalue"],
            row["evalue"] <= MAX_EVALUE_THRESHOLD,
        ])


class _BlastnOutput6Schema:
    """
    blastn output format 6 as documented in
    http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    it's also the format of our GSNAP and RAPSEARCH2 output
    """
    SCHEMA = [
        ("qseqid", str),
        ("sseqid", str),
        ("pident", float),
        ("length", int),
        ("mismatch", int),
        ("gapopen", int),
        ("qstart", int),
        ("qend", int),
        ("sstart", int),
        ("send", int),
        ("evalue", float),
        ("bitscore", float),
    ]

class _BlastnOutput6NTSchema:
    """
    Additional blastn output columns.
    """
    SCHEMA = _BlastnOutput6Schema.SCHEMA + [
        ("qlen", int),      # query sequence length, helpful for computing qcov
        ("slen", int),      # subject sequence length, so far unused in IDseq
    ]

class _BlastnOutput6NTRerankedSchema:
    """
    Re-ranked output of blastn.  One row per query.  Two additional columns.
    """
    SCHEMA = _BlastnOutput6NTSchema.SCHEMA + [
        ("qcov", float),     # fraction of query covered by the optimal set of HSPs
        ("hsp_count", int),   # cardihnality of optimal fragment cover;  see BlastCandidate
    ]

class BlastnOutput6NTRerankedReader(_BlastnOutput6NTRerankedSchema, _BlastnOutput6ReaderBase):
    def __init__(self,f:Iterable[Text],filter_invalid:bool=True,min_alignment_length :int=0):
        super().__init__(f,self.SCHEMA,filter_invalid,min_alignment_length)

def parser_blastn_m8(input_blastn_m8_path):
    with open(input_blastn_m8_path) as input_blastn_m8_f:
        for row in BlastnOutput6NTRerankedReader(input_blastn_m8_f,filter_invalid=True,min_alignment_length=50):
           read_id,accession_id,evalue =row['qseqid'],row["sseqid"],row["evalue"]




def parser_blast_outm8(gsnap_out_file_m8):

    """
    length
    pident
    evalue
    hit_summary_writer.writerow({
        "read_id": read_id,
        "level": hit_level,
        "taxid": taxid,
        "accession_id": best_accession_id,
        "species_taxid": species_taxid,
        "genus_taxid": genus_taxid,
        "family_taxid": family_taxid,
    })
    """
    pass
    with open (gsnap_out_file_m8,'r') as h1:
        for line in h1:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore=line.strip().split("\t")

def iterator_taxid_lineage_db(taxid_lineage_db):
    "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210825_taxid_lineage_db"
    with open (taxid_lineage_db,'r') as  db_path_f:
        #next(db_path_f)
        for line in db_path_f:
            if line.startswith("speciess_taxid"):
                pass
            else:
                taxid ,name_lineage,id_lineage = line.strip().split("\t")
                one_taid_d = dict(zip([name for name in name_lineage.split(";")],[id for id in id_lineage.split(";")]))
                yield {taxid:one_taid_d}

def get_lineage(accession_id,taxid_lineage_db):
    """
    :param accession_id:   自定义构建的数据库， kraken:taxid|1125630|NC_016845.1
    :return:

    NC_016845.1   ----》 taxid|1125630  ———》 taxid_lineage_db
    """
    #taxid = re.search(r'\|(^\d+\d$)\|',accession_id).groups()[0]
    taxid = accession_id.split("|")[1]
    for d in iterator_taxid_lineage_db(taxid_lineage_db):
        if taxid in d.keys():
            print(taxid,d.get(taxid))


def get_accession_id2species_name(accession_id,taxid_lineage_db):
    taxid = accession_id.split("|")[1]
    for d in iterator_taxid_lineage_db(taxid_lineage_db):
        if taxid in d.keys():
            species_name,species_name_taxid = [(k,v) for k,v in d.get(taxid).items()][-1]  ###暂且 -1
            return (species_name,species_name_taxid)


def get_accession_id2genus_name(accession_id,taxid_lineage_db):
    taxid = accession_id.split("|")[1]
    for d in iterator_taxid_lineage_db(taxid_lineage_db):
        if taxid in d.keys():
            genus_name,genus_name_taxid = [(k,v) for k,v in d.get(taxid).items()][-2]  ###暂且 -2
            return (genus_name,genus_name_taxid)



lineage_db_iterator =iterator_taxid_lineage_db(taxid_lineage_db="/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db")
# print(lineage_db_iterator.__next__())
# print(lineage_db_iterator.__next__())
# print(lineage_db_iterator.__next__())

# for d in lineage_db_iterator:
#     print(d)

get_lineage(accession_id="kraken:taxid|573|NC_016845.1",taxid_lineage_db="/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db")

print(get_accession_id2species_name(accession_id="kraken:taxid|573|NC_016845.1",taxid_lineage_db="/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db"))
print(get_accession_id2genus_name(accession_id="kraken:taxid|573|NC_016845.1",taxid_lineage_db="/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db"))
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


blast_outm8 = '/mnt/home/huanggy/project/meta/result/align/my_gsnap.out'
#{('Klebsiella pneumoniae', '574'):1}
# blast_out = {('Klebsiella pneumoniae', '573'):0}
# blast_out.update({('Klebsiella pneumoniae', '574'):1})
# print(blast_out)
blast_out_s = {}
blast_out_g = {}
with open(blast_outm8,"r") as blast_outm8_f:
    total_counts = 0
    for line in blast_outm8_f:
        qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.split("\t")
        s_level = get_accession_id2species_name(sseqid,taxid_lineage_db = "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db")
        if s_level not in blast_out_s:
            # qseqid_list = []
            # qseqid_list.append(qseqid)
            blast_out_s.update({s_level:{"species_id":s_level[1],"species_name":s_level[0],"count":1,"qseqid":[qseqid]}})
        else:
            counts = blast_out_s[s_level]["count"] + 1
            #blast_out.update({s_level}.update({"count":counts}))
            blast_out_s[s_level].update({"count":counts})
            ###同一seqid 比对到不同的种水平 ，但是是同一个属
            # if qseqid not in blast_out_s[s_level]["qseqid"]:
            #     blast_out_s[s_level]["qseqid"].append(qseqid)

        total_counts += 1

print(blast_out_s)
print(total_counts)

for k , v in blast_out_s.items():
    s_abundance = "{:.2%}".format(float(v["count"])/float(total_counts))
    if s_abundance not in blast_out_s[k].keys():
        blast_out_s[k].update({"s_abundance":s_abundance})
    else:
        pass
print(blast_out_s)



# with open(blast_outm8,"r") as blast_outm8_f:
#     total_counts = 0
#     for line in blast_outm8_f:
#         qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.split("\t")
#         g_level = get_accession_id2genus_name(sseqid,taxid_lineage_db = "/mnt/data/NCBI_Refseq/Bacteria/refseq/unzip/select_assembly_info_20210826_taxid_lineage_db")
#         g_level_name = g_level[0]
#         g_level_id = g_level[1]
#         if g_level_name not in blast_out_g:
#             blast_out_g.update({g_level_name:{"g_level_id":g_level_id,"genus_name":g_level_name,"count":1,"qseqid":[qseqid]}})
#         else:
#             ##不同seqid 但是在同一个属  ,此时该属 +1 条hits
#             ##同一seqid 比对到不同的种水平 ，但是是同一个属 ，此时属于该属水平的 不能加1 ，属于重复reads
#             #blast_out.update({s_level}.update({"count":counts}))
#             if qseqid not in blast_out_g[g_level_name]["qseqid"]:
#                 counts = blast_out_g[g_level_name]["count"] + 1
#                 blast_out_g[g_level_name].update({"count":counts})
#                 blast_out_g[g_level_name]["qseqid"].append(qseqid)
#             else:
#                 #此reads 已存在该属里面，不用append ，counts也不用 +1
#                 pass
#
#         total_counts += 1
#
#
#
# #print(blast_out_g)
# print(total_counts)
# for k , v in blast_out_g.items():
#     g_abundance = "{:.2%}".format(float(v["count"])/float(total_counts))
#     del blast_out_g[k]["qseqid"]
#     if g_abundance not in blast_out_g[k].keys():
#         blast_out_g[k].update({"g_abundance":g_abundance})
#     else:
#         pass
#
# print(blast_out_g)
#
# import json
# blast_out_g_json = json.dumps(blast_out_g)
# print(blast_out_g_json)


# with open('/mnt/home/huanggy/project/meta/blast_out_g.json','w',encoding='utf-8') as fObj:
#     json.dump(blast_out_g_json, fObj, ensure_ascii=False)

# out = open ('/mnt/home/huanggy/project/meta/blast_out_g.json','w')
# for k,v in blast_out_g.items():
#     pass










