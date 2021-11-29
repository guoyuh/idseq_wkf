import pandas as pd 
df = pd.read_csv('/mnt/home/huanggy/idseq_dag/report/samples.txt',sep = "\t",header=0)
print(df)
SAMPLES = df["sample_name"].tolist()
print(SAMPLES)