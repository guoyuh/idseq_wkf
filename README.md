# idseq_wkf

# install
	##conda env create -f snakemake_py36.yml
# example
	##snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py --cores 2 -p --directory ~/project/20211124_test  --configfile ~/idseq_dag/runs/config.yaml
# 流程图
	##snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py  --dag --cores 1 | dot -Tpdf > test.pdf
# idseq_wkf
