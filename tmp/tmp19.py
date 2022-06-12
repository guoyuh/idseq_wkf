

import pandas as pd
uniq_kreport_mpa="/mnt/home/huanggy/project/PRJNA558701/result/S57_SRR9903807/kraken2/S57_SRR9903807.kreport.mpa.rpm.xls"
tb = pd.read_table(uniq_kreport_mpa,sep='\t',header=0)
header = tb.columns.to_list()

print(tb)
tb2 = tb.drop_duplicates(keep='first')
print("=======================")
print(tb2)