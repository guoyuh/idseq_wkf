import glob

SAMPLE=["S200026380_L01_90","S200023932_L01_97"]
ind=["90","97"]
OUTDIR="/mnt/home/huanggy/project/20211117"
host_yanhaung="/mnt/data/yanhuang"
PriceSeqFilter="/usr/local/bin/PriceSeqFilter"
kneaddata="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata"
kneaddata_read_count_table="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata_read_count_table"
KRAKEN="/mnt/project/tools/kraken2/kraken2"
Bracken="/mnt/project/tools/Bracken/bracken"
rule all:
    ####想要哪一步的结果,都在这里买定义
    input:
        expand("{outdir}/result/{sample}/filter",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/filter/{sample}.f.fastq",outdir=OUTDIR,sample=SAMPLE)

rule kneaddata:
    input:
        #"/mnt/data/NGSDATA/meta_data/S200026380_L01_90.fq.gz",
        #"/mnt/data/NGSDATA/meta_data/S200023932_L01_97.fq.gz"
        #single = expand("/mnt/data/NGSDATA/meta_data/{sample}.fq.gz",sample=SAMPLE)
        single = "/mnt/data/NGSDATA/meta_data/{sample}.fq.gz"
    output:
        "{outdir}/result/{sample}/filter",
        "{outdir}/result/{sample}/filter/{sample}.fastq"

    params:
        thread_num = 16
    shell:
        "kneaddata \
        -i {input.single} \
        -o {output[0]}  \
        -v -t {params.thread_num} --remove-intermediate-output \
        --trimmomatic /mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10' \
        --bowtie2-options '--very-sensitive --dovetail' -db {host_yanhaung} \
        --output-prefix {wildcards.sample}"
##### ????????????????????shell 中用sample 这个变量 必须要用wildcards.sample ,单独用{sample} 报错如下：
# RuleException in line 27 of /mnt/home/huanggy/idseq_dag/report/mNGS.snakemake.py:
# NameError: The name 'sample' is unknown in this context. Did you mean 'wildcards.sample'?

rule kneaddata_stat:
    input:
        #expand("{outdir}/result/{sample}/filter",outdir=OUTDIR,sample=SAMPLE)
        "{outdir}/result/{sample}/filter"
    output:
        kneaddata_sum = "{outdir}/result/{sample}/filter/kneaddata.sum.txt"
    shell:
        "{kneaddata_read_count_table} --input {input} --output {output.kneaddata_sum}"


# rule PriceSeqFilter:
#     input:
#         "{outdir}/result/{sample}/filter/{sample}.fastq"
#     output:
#         "{outdir}/result/{sample}/filter/{sample}.f.fastq"
#     shell:
#         "{PriceSeqFilter} \
#         -a 12 -rnf 90 -log c -f {outdir}/result/{sample}/filter/{sample}.fastq \
#         -o {outdir}/result/{wildcards.sample}/filter/{wildcards.sample}.f.fastq -rqf 85 0.98 -lenf 36 && rm -f {outdir}/result/{wildcards.sample}/filter/{wildcards.sample}.fastq "



rule PriceSeqFilter:
    input:
        #"{outdir}/result/{sample}/filter/kneaddata.sum.txt",
        "{outdir}/result/{sample}/filter/{sample}.fastq"
    output:
        "{outdir}/result/{sample}/filter/{sample}.f.fastq"
    shell:
        "{PriceSeqFilter} -a 12 -rnf 90 -log c -f {input} -o {output} -rqf 85 0.98 -lenf 36 && rm -f {input}"
    
    