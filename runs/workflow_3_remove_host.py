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
Coverm="/mnt/home/huanggy/bin/coverm"
amrfinder="/mnt/project/tools/amrfinder-3.10.18/amrfinder"

bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral"
refsesq_dir="/mnt/data/NCBI_Refseq/bac_fungi_protozoa_viral"
host_yanhaung="/mnt/data/yanhuang" ### 由炎黄1号 + somatic-b37_Homo_sapiens_assembly19.fasta + GRCh38_latest_rna
host_hg37="/mnt/db/kneaddata/human_genome" ### 由hg37 + human_contamination  组成from 官网
kraken2_db="/mnt/data/kraken2_mask"




out_dir = sys.argv[1]
sample = sys.argv[2]
fq = sys.argv[3]


# create a DAG object
wkf1 = DAG("workflow_1")
wkf2 = DAG("workflow_2")
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

##################################begin with nohost##########################################
kraken2_classify = Task(
    id="kraken2_classify",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "{KRAKEN} --db {db} {input_clean_fq} \
        --threads 16 \
        --report {outdir}/result/{sample}/kraken2/{sample}.report \
        --output {outdir}/result/{sample}/kraken2/{sample}.output"
        .format(KRAKEN=KRAKEN,db=kraken2_db,outdir = out_dir,sample=sample,input_clean_fq = fq)
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
    script = "python  ~/idseq_wkf/scripts/bracken_reads_rpm.py -w {outdir}  -p {sample} -c {cleanfq}".format(outdir=out_dir,sample=sample,cleanfq = fq)
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
    script = """
    python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {outdir}/result/{sample}/kraken2/{sample}.kreport -o {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt
    """.format(outdir=out_dir,sample=sample)
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

kreport2mpa2rpm_BSD = Task(
    #kreport2mpa2rpm_BSD  取代kreport2mpa2rpm  添加BSD
    id="kreport2mpa2rpm_BSD",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = """
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport -o {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt 
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport -o {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_bsd_cleanfq.py {clean_fq} {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir} {sample} &
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_bsd_cleanfq.py {clean_fq} {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt {outdir} {sample} &
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_bsd_cleanfq.py {clean_fq} {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt {outdir} {sample} &
    """.format(outdir = out_dir,sample=sample,clean_fq = fq)
)

clean_reads_bowtie2_bac_fungi_protozoa_viral = Task(
    id="clean_reads_bowtie2_bac_fungi_protozoa_viral",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = "{bowtie2_bin} -q --quiet --sensitive --threads {threads} -x {bac_fungi_protozoa_viral_index} \
    -U {clean_fq} --al-gz {outdir}/result/{sample}/align/{sample}.clean.refseq.fq.gz -S {outdir}/result/{sample}/align/{sample}.refseq.sam \
    1>{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log 2>&1".format(bowtie2_bin = bowtie2_bin,threads=64,bac_fungi_protozoa_viral_index=bac_fungi_protozoa_viral_index,clean_fq = fq,outdir=out_dir,sample=sample)
)

sam2bam = Task(
    id="sam2bam",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = """
    {samtools} view -bS {outdir}/result/{sample}/align/{sample}.refseq.sam > {outdir}/result/{sample}/align/{sample}.bam
    {samtools} sort -o {outdir}/result/{sample}/align/{sample}.sorted.bam {outdir}/result/{sample}/align/{sample}.bam
    rm -f {outdir}/result/{sample}/align/{sample}.refseq.sam &
    rm -f {outdir}/result/{sample}/align/{sample}.bam &
    """.format(samtools=samtools,outdir=out_dir,sample=sample)
)

Covem = Task(
    id="Covem",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = """
    {Coverm} genome -m mean covered_bases covered_fraction variance length count rpkm --min-covered-fraction 0 --bam-files {outdir}/result/{sample}/align/{sample}.sorted.bam --genome-fasta-directory {refsesq_dir} --genome-fasta-extension fna > {outdir}/result/{sample}/align/{sample}.cove.tsv
    python ~/idseq_wkf/scripts/coverm_format.py {outdir}/result/{sample}/align/{sample}.cove.tsv /mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20220406.txt /mnt/data/NCBI_Refseq/viral/kraken2_viral_refseq_info
    python ~/idseq_wkf/scripts/species_coverm_stat_format.py {outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.rpm.xls  {outdir}/result/{sample}/align/{sample}.refseq.cove.xls &
    python ~/idseq_wkf/scripts/species_coverm_stat_format.py {outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.rpm.xls  {outdir}/result/{sample}/align/{sample}.refseq.cove.xls &
    python ~/idseq_wkf/scripts/species_coverm_stat_format.py {outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.rpm.xls  {outdir}/result/{sample}/align/{sample}.refseq.cove.xls &
    """.format(Coverm=Coverm,outdir=out_dir,sample=sample,refsesq_dir= refsesq_dir)
)

SPAdes = Task(
    id="SPAdes",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    #SPAdes
    # if --meta  giveme error Sorry, current version of metaSPAdes can work either with single library (paired-end only) or in hybrid paired-end + (TSLR or PacBio or Nanopore) mode.
    #"/mnt/project/tools/SPAdes-3.15.3/bin/spades.py -s {input[0]} -o {params.out_dir}"
    script = """
        mkdir -p {outdir}/result/{sample}/filter/spades_out/
        /mnt/project/tools/SPAdes-3.15.3/bin/spades.py -s {outdir}/result/{sample}/filter/{sample}.f.fastq -o {outdir}/result/{sample}/filter/spades_out/contigs.fasta""".format(outdir=out_dir,sample=sample)
)
amrfinder = Task(
    id="amrfinder",  # your task id, should be unique
    work_dir=out_dir + "/bin/" + sample,  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option="-step 5",  # the option of "sge" or "local"
    script = """
        {amrfinder} -n {outdir}/result/{sample}/filter/spades_out/contigs.fasta -o {outdir}/result/{sample}/filter/spades_out/{sample}.amr""".format(amrfinder=amrfinder,outdir=out_dir,sample=sample)
)


def run_negative_control():
    # when you create a task, then add it to DAG object
    wkf1.add_task(mkdir)
    # wkf1.add_task(host_remove)
    # host_remove.set_downstream(priceSeqFilter)
    # wkf1.add_task(priceSeqFilter)
    # priceSeqFilter.set_downstream(kneaddata_stat)
    # wkf1.add_task(kneaddata_stat)
    mkdir.set_downstream(kraken2_classify)
    wkf1.add_task(kraken2_classify)
    kraken2_classify.set_downstream(bracken_Abundance_estimation)
    wkf1.add_task(bracken_Abundance_estimation)
    bracken_Abundance_estimation.set_downstream(bracken_format_rpm)
    wkf1.add_task(bracken_format_rpm)
    bracken_format_rpm.set_downstream(kout2kreport)
    wkf1.add_task(kout2kreport)
    kout2kreport.set_downstream(kreport2mpa)
    wkf1.add_task(kreport2mpa)
    kreport2mpa.set_downstream(kreport2krona)
    wkf1.add_task(kreport2krona)
    # kreport2krona.set_downstream(kraken2_filter)
    # wkf1.add_task(kraken2_filter)
    # kraken2_filter.set_downstream(kfiltout2kreport)
    # wkf1.add_task(kfiltout2kreport)
    # kfiltout2kreport.set_downstream(kreport2mpa2rpm)
    # wkf1.add_task(kreport2mpa2rpm)
    # all of you tasks were added to you workflow, you can run it
    do_dag(wkf1)


def run_clin():
    # when you create a task, then add it to DAG object
    wkf2.add_task(mkdir)
    # wkf2.add_task(host_remove)
    # host_remove.set_downstream(priceSeqFilter)
    # wkf2.add_task(priceSeqFilter)
    # priceSeqFilter.set_downstream(kneaddata_stat)
    # wkf2.add_task(kneaddata_stat)
    mkdir.set_downstream(kraken2_classify)
    wkf2.add_task(kraken2_classify)
    kraken2_classify.set_downstream(bracken_Abundance_estimation)
    wkf2.add_task(bracken_Abundance_estimation)
    bracken_Abundance_estimation.set_downstream(bracken_format_rpm)
    wkf2.add_task(bracken_format_rpm)
    bracken_format_rpm.set_downstream(kout2kreport)
    wkf2.add_task(kout2kreport)
    kout2kreport.set_downstream(kreport2mpa)
    wkf2.add_task(kreport2mpa)
    kreport2mpa.set_downstream(kreport2krona)
    wkf2.add_task(kreport2krona)

    print("======================only clinsample run task========================")
    kreport2krona.set_downstream(kraken2_filter)
    wkf2.add_task(kraken2_filter)
    kraken2_filter.set_downstream(kfiltout2kreport)
    wkf2.add_task(kfiltout2kreport)
    kfiltout2kreport.set_downstream(kreport2mpa2rpm_BSD)
    wkf2.add_task(kreport2mpa2rpm_BSD)
    kreport2mpa2rpm_BSD.set_downstream(clean_reads_bowtie2_bac_fungi_protozoa_viral)
    wkf2.add_task(clean_reads_bowtie2_bac_fungi_protozoa_viral)
    clean_reads_bowtie2_bac_fungi_protozoa_viral.set_downstream(sam2bam)
    wkf2.add_task(sam2bam)
    sam2bam.set_downstream(Covem)
    wkf2.add_task(Covem)


    # wkf2.add_task(SPAdes)
    # SPAdes.set_downstream(amrfinder)
    # wkf2.add_task(amrfinder)
    #all of you tasks were added to you workflow, you can run it
    do_dag(wkf2)



"""
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 test01 /mnt/home/huanggy/raw_fq/2022-03/HZ200A014_20220316/21JS944006-1DL-DH-UDB-136_UDB-136.fq.gz
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 test02 /mnt/home/huanggy/raw_fq/2022-03/HZ200A014_20220316/21JS944009-1DL-DH-UDB-135_UDB-135.fq.gz

test
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 S200023771_L01_101 /mnt/data/NGSDATA/meta_data/S200023771_L01_101.fq.gz
python ~/idseq_wkf/runs/workflow_2.py /mnt/home/huanggy/project/20220417 S200023771_L01_102 /mnt/data/NGSDATA/meta_data/S200023771_L01_102.fq.gz

批次分析思路：sub_sh,每批次阴性控制先一步运行，后续临床样本并行运行
python ~/idseq_wkf/runs/workflow_2.py project_dir 阴性控制 fq 
python ~/idseq_wkf/runs/workflow_2.py project_dir clin_sp1 fq &
python ~/idseq_wkf/runs/workflow_2.py project_dir clin_sp2 fq &
...
"""


def main():
    if len(sys.argv) != 4:
        print("Usage:python sys.argv[0] <outdir> <sample_id> <fq>")
    else:
        if sys.argv[2].find("yinkong") != -1:
            print("this is 阴控样本%s"%(sys.argv[2]))
            run_negative_control()
        else:
            run_clin()


if __name__ == '__main__':
    main()
