# idseq_wkf

# install
	##conda env create -f snakemake_py36.yml
# example
	##snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py --cores 2 -p --directory ~/project/20211124_test  --configfile ~/idseq_dag/runs/config.yaml
# 流程图
	##snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py  --dag --cores 1 | dot -Tpdf > test.pdf
# idseq_wkf
	##you also run by workflow.py like this :
	python /path/to/idseq_wkf/runs/workflow_3.py /path/to/outdir  sample_prefix  /path/to/single.fq.gz &
