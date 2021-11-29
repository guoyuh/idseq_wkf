# idseq_dag

# install


# example
snakemake -s ~/idseq_dag/report/mNGS.snakemake.py --cores 2 -p --directory ~/project/20211124_test  --configfile ~/idseq_dag/report/config.yaml
# 流程图
snakemake -s ~/idseq_dag/report/mNGS.snakemake.py  --dag --cores 1 | dot -Tpdf > test.pdf
# idseq_dag
