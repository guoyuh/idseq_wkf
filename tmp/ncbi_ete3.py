from ete3 import NCBITaxa                       # 导入此模块
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()                 # 升级
print(ncbi.get_rank([93061,9443,1282,190485]))
print(ncbi.get_taxid_translator([9606]))
print(ncbi.get_name_translator(['Lactococcus lactis subsp. lactis']))



def get_taxon_path(taxid):
    try:
        taxid_list = ncbi.get_lineage(taxid)
    except ValueError:
        raise
    kept_levels = ["superkingdom", "phylum", "class",
                   "order", "family", "genus", "species", "strain"]

    rank_dict = ncbi.get_rank(taxid_list)
    kept_taxids = []

    for level in kept_levels:
        for k, v in rank_dict.items():
            if v == level:
                kept_taxids.append(k)
    taxsn_dict = ncbi.get_taxid_translator(kept_taxids)
    taxid_path = "|".join(map(str, kept_taxids))
    taxsn_path = "|".join([taxsn_dict[tax] for tax in kept_taxids])
    return [taxid_path, taxsn_path]


print(get_taxon_path(487214))

print(get_taxon_path(5061))
print(get_taxon_path(425011))

print(get_taxon_path(1566990))