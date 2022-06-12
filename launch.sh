#!/bin/bash
out_put=$1
config_file=$2
thread_num=$3
if [ $thread_num -lt 4 ];then
	thread_num=$3
else
	thread_num=4
fi

help(){
	echo 'USAGE: snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py --cores 4 -p --directory path/to/outdir --configfile ~/idseq_wkf/runs/config.yaml -r'
}
run(){
    snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py --cores ${thread_num} -p --directory ${out_put} --configfile  ${config_file} -r
}


if [ $# -lt 2 ];then
	help
	exit
else
	echo "snakemake -s ~/idseq_wkf/runs/mNGS.snakemake.py --cores ${thread_num} -p --directory ${out_put} --configfile  ${config_file} -r"
	run
fi

if [ $? -eq 0 ];then
	echo "merge zhikong xinxi "
fi
