#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


import os,sys

from DAGflow.dagflow.dag import DAG,Task,ParallelTask

from DAGflow.dagflow import do_dag

PriceSeqFilter="/usr/local/bin/PriceSeqFilter"
kneaddata="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata"
kneaddata_read_count_table="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata_read_count_table"
bowtie2_path="/mnt/project/tools/bowtie2-2.4.4"
bowtie2_bin="/mnt/project/tools/bowtie2-2.4.4/bowtie2"
KRAKEN="/mnt/project/tools/kraken2/kraken2"
Bracken="/mnt/project/tools/Bracken/bracken"
samtools="/mnt/db/software/samtools-1.9/samtools"

#bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral"
refsesq_dir="/mnt/data/NCBI_Refseq/bac_fungi_protozoa_viral"
host_yanhaung="/mnt/data/yanhuang" ### 由炎黄1号 + somatic-b37_Homo_sapiens_assembly19.fasta + GRCh38_latest_rna
host_hg37="/mnt/db/kneaddata/human_genome" ### 由hg37 + human_contamination  组成from 官网
kraken2_db="/mnt/data/kraken2_mask"
Coverm="/mnt/home/huanggy/bin/coverm"

out_dir = sys.argv[1]
sample = sys.argv[2]
fq = sys.argv[3]



# create a DAG object
my_dag = DAG("mNGS")
# create the first task 'make_db'
mkdir = Task(
    id="mkdir",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 1",  # the option of "sge" or "local"
    script = "mkdir -p %s/result/%s/filter && mkdir -p %s/result/%s/align && mkdir -p %s/result/%s/kraken2 "% (out_dir,sample,out_dir,sample,out_dir,sample)
)

# host_remove = Task(
#     id="kneaddata",  # your task id, should be unique
#     work_dir=out_dir + "/bin/" + sample,  # you task work directory
#     type="local",  # the way your task run. if "sge", task will submit with qsub
#     option="-step 2",  # the option of "sge" or "local"
#     script = "echo %s -i %s -o %s/result/%s/filter -v -t 16 \
#     --remove-intermediate-output \
#     --trimmomatic /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2 \
#     --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10' \
#     --bowtie2 %s --bowtie2-options '--very-sensitive --dovetail' -db %s \
#     --output-prefix %s" %(kneaddata,fq,out_dir,sample,bowtie2_path,host_hg37,sample)
# )

host_remove = Task(
    id="kneaddata",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 2",  # the option of "sge" or "local"
    script = "{kneaddata} \
        -i {input_single} \
        -o {outdir}/result/{sample}/filter  \
        -v -t 16 --remove-intermediate-output \
        --trimmomatic /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10' \
        --bowtie2 {bowtie2_path} --bowtie2-options '--very-sensitive --dovetail' -db {host_hg37} \
        --output-prefix {sample}".format(kneaddata=kneaddata,input_single = fq,outdir= out_dir,sample=sample,bowtie2_path=bowtie2_path,host_hg37=host_hg37,)
)

priceSeqFilter = Task(
    id="PriceSeqFilter",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 3",  # the option of "sge" or "local"
    script = "PriceSeqFilter -a 12 -rnf 90 -log c -f {0}/result/{1}/filter/{2}.fastq -o {3}/result/{4}/filter/{5}.f.fastq \
     -rqf 85 0.98 -lenf 36".format(out_dir,sample,sample,out_dir,sample,sample)
)
"""
input:
    "{outdir}/result/{sample}/filter/{sample}.fastq"
output:
    "{outdir}/result/{sample}/filter/{sample}.f.fastq"
"""

kneaddata_stat = Task(
    id="kneaddata_stat",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 4",  # the option of "sge" or "local"
    script = "{kneaddata_read_count_table} --input {outdir}/result/{sample}/filter --output {outdir}/result/{sample}/filter/kneaddata.sum.txt  && rm -f {outdir}/result/{sample}/filter/{sample}_hg37dec_v0.1_bowtie2_contam.fastq"
        .format(kneaddata_read_count_table=kneaddata_read_count_table,outdir=out_dir,sample=sample)
)


kraken2_classify = Task(
    id="kraken2_classify",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "{KRAKEN} --db {db} {outdir}/result/{sample}/filter/{sample}.f.fastq \
        --threads 16 \
        --report {outdir}/result/{sample}/kraken2/{sample}.report \
        --output {outdir}/result/{sample}/kraken2/{sample}.output"
        .format(KRAKEN=KRAKEN,db=kraken2_db,outdir = out_dir,sample=sample)
)


bracken_Abundance_estimation = Task(
    id="bracken_Abundance_estimation",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "{Bracken} -d {db} -i {outdir}/result/{sample}/kraken2/{sample}.report \
        -o {outdir}/result/{sample}/kraken2/{sample}.s.bracken \
        -w {outdir}/result/{sample}/kraken2/{sample}.s.bracken.report -r 75 -l S \
        && {Bracken} -d {db} -i {outdir}/result/{sample}/kraken2/{sample}.report \
        -o {outdir}/result/{sample}/kraken2/{sample}.g.bracken \
        -w {outdir}/result/{sample}/kraken2/{sample}.g.bracken.report -r 75 -l G"
        .format(Bracken=Bracken,db=kraken2_db,outdir=out_dir,sample=sample)
)


bracken_format_rpm = Task(
    id="bracken_format_rpm",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python  ~/idseq_wkf/scripts/bracken_reads_rpm.py -w {outdir}  -p {sample}".format(outdir=out_dir,sample=sample)
)

kout2kreport = Task(
    id="kout2kreport",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python ~/idseq_wkf/KrakenTools/make_kreport.py \
    -i {outdir}/result/{sample}/kraken2/{sample}.output -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt \
    -o {outdir}/result/{sample}/kraken2/{sample}.kreport ".format(outdir=out_dir,sample=sample)
)

kreport2mpa = Task(
    id="kreport2mpa",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  \
    -r {outdir}/result/{sample}/kraken2/{sample}.kreport \
    -o {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt".format(outdir=out_dir,sample=sample)
)

kreport2krona = Task(
    id="kreport2krona",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python  ~/idseq_wkf/KrakenTools/kreport2krona.py --report-file  {outdir}/result/{sample}/kraken2/{sample}.kreport \
    -o {outdir}/result/{sample}/kraken2/{sample}.krona \
    && perl ~/idseq_wkf/Krona/bin/ktImportText {outdir}/result/{sample}/kraken2/{sample}.krona  \
    -o {outdir}/result/{sample}/kraken2/{sample}.krona.html".format(outdir=out_dir,sample=sample)
)

kraken2_filter = Task(
    id="kraken2_filter",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python ~/idseq_wkf/scripts/kraken2_uniq_reads_filter.py {outdir}/result/{sample}/kraken2/{sample}.output".format(outdir=out_dir,sample=sample)
)

kfiltout2kreport = Task(
    id="kfiltout2kreport",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "python ~/idseq_wkf/KrakenTools/make_kreport.py -i {outdir}/result/{sample}/kraken2/{sample}.uniq_output -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt -o {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport && \
    python ~/idseq_wkf/KrakenTools/make_kreport.py -i {outdir}/result/{sample}/kraken2/{sample}.uniq80_output -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt -o {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport".format(outdir=out_dir,sample=sample)
)

kreport2mpa2rpm = Task(
    id="kreport2mpa2rpm",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = """
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport -o {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt 
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport -o {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {outdir}/result/{sample}/filter/kneaddata.sum.txt {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {outdir}/result/{sample}/filter/kneaddata.sum.txt {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt
    """.format(outdir = out_dir,sample=sample)
)



# when you create a task, then add it to DAG object
my_dag.add_task(mkdir)
my_dag.add_task(host_remove)
host_remove.set_downstream(priceSeqFilter)
my_dag.add_task(priceSeqFilter)
priceSeqFilter.set_downstream(kneaddata_stat)
my_dag.add_task(kneaddata_stat)
kneaddata_stat.set_downstream(kraken2_classify)
my_dag.add_task(kraken2_classify)
kraken2_classify.set_downstream(bracken_Abundance_estimation)
my_dag.add_task(bracken_Abundance_estimation)
bracken_Abundance_estimation.set_downstream(bracken_format_rpm)
my_dag.add_task(bracken_format_rpm)
bracken_format_rpm.set_downstream(kout2kreport)
my_dag.add_task(kout2kreport)
kout2kreport.set_downstream(kreport2mpa)
my_dag.add_task(kreport2mpa)
kreport2mpa.set_downstream(kreport2krona)
my_dag.add_task(kreport2krona)
kreport2krona.set_downstream(kraken2_filter)
my_dag.add_task(kraken2_filter)
kraken2_filter.set_downstream(kfiltout2kreport)
my_dag.add_task(kfiltout2kreport)
kfiltout2kreport.set_downstream(kreport2mpa2rpm)
my_dag.add_task(kreport2mpa2rpm)
# all of you tasks were added to you workflow, you can run it
do_dag(my_dag)


"""
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 test01 /mnt/home/huanggy/raw_fq/2022-03/HZ200A014_20220316/21JS944006-1DL-DH-UDB-136_UDB-136.fq.gz
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 test02 /mnt/home/huanggy/raw_fq/2022-03/HZ200A014_20220316/21JS944009-1DL-DH-UDB-135_UDB-135.fq.gz

test
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 S200023771_L01_101 /mnt/data/NGSDATA/meta_data/S200023771_L01_101.fq.gz
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 S200023771_L01_102 /mnt/data/NGSDATA/meta_data/S200023771_L01_102.fq.gz

"""