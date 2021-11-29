

import pandas as pd
import glob


#glob.glob("/mnt/home/huanggy/project/20210926_test/result/kraken2")

tb = pd .read_table("/mnt/home/huanggy/project/20210926_test/result/kraken2/S200023771_L01_101.s.bracken",sep="\t",header = 0 )




tb["id"] = ["S200023771_L01_101"]*len(tb)
print(tb)
print("=========")

tb2 = pd .read_table("/mnt/home/huanggy/project/20210926_test/result/kraken2/S200023771_L01_102.s.bracken",sep="\t",header = 0 )
tb2["id"] = ["S200023771_L01_102"]*len(tb2)
print(tb2)
print("=========")
tb3 = tb.append(tb2,ignore_index=False)
print(tb3)