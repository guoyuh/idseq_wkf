
from collections import OrderedDict

db1="/mnt/data/kraken2_mask/taxonomy/mydb_taxonomy.txt"
db2="/mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220407.txt"

uni = OrderedDict()
out = open("/mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt","w")

with open (db1,"r") as fh:
    for line in fh:
        line = line.strip()
        if line not in uni:
            uni.update({line:1})
        else:
            pass

with open (db2,"r") as fh:
    for line in fh:
        line = line.strip()
        if line not in uni:
            uni.update({line:1})
        else:
            pass


for k,v in uni.items():
    print(k,file=out)

out.close()
