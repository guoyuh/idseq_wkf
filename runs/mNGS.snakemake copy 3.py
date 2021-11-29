SAMPLE=["S200026380_L01_90","S200023932_L01_97"]
ind=["90","97"]
OUTDIR="/mnt/home/huanggy/project/20211117"
host_yanhaung="/mnt/data/yanhuang"
PriceSeqFilter="/usr/local/bin/PriceSeqFilter"
kneaddata="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata"
kneaddata_read_count_table="/mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/kneaddata_read_count_table"
bowtie2="/mnt/project/tools/bowtie2-2.4.4"
KRAKEN="/mnt/project/tools/kraken2/kraken2"
Bracken="/mnt/project/tools/Bracken/bracken"

"""
$ ll /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/
总用量 89492
drwxrwxr-x 2 huanggy huanggy     4096 11月 17 20:16 ./
drwxrwxr-x 4 huanggy huanggy     4096 11月 17 20:15 ../
-rw-rw-r-- 1 huanggy huanggy  2124589 11月 17 20:16 S200023932_L01_97.fastq
-rw-rw-r-- 1 huanggy huanggy     8219 11月 17 20:16 S200023932_L01_97.log
-rw-rw-r-- 1 huanggy huanggy 30453478 11月 17 20:15 S200023932_L01_97.repeats.removed.fastq
-rw-rw-r-- 1 huanggy huanggy 30705093 11月 17 20:15 S200023932_L01_97.trimmed.fastq
-rw-rw-r-- 1 huanggy huanggy 28328889 11月 17 20:16 S200023932_L01_97_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq

Decompressing gzipped file ...

Decompressing gzipped file ...

Initial number of reads ( /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/decompressed_armkzm65_S200026380_L01_90.fq ): 250000.0
Running Trimmomatic ...

java -Xmx500m -jar /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/trimmomatic.jar SE -threads 16 -phred33 /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/decompressed_armkzm65_S200026380_L01_90.fq /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.trimmed.fastq SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10

Initial number of reads ( /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/decompressed_hpfz9fet_S200023932_L01_97.fq ): 250000.0
Running Trimmomatic ...

java -Xmx500m -jar /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/trimmomatic.jar SE -threads 16 -phred33 /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/decompressed_hpfz9fet_S200023932_L01_97.fq /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.trimmed.fastq SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10

Total reads after trimming ( /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.trimmed.fastq ): 209728.0
Total reads after trimming ( /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.trimmed.fastq ): 231551.0
Running trf ...

kneaddata_trf_parallel --input /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.trimmed.fasta --output /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.trimmed.fasta.trf.parameters.2.7.7.80.10.50.500.dat --trf-path /mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/trf --trf-options '2 7 7 80 10 50 500 -h -ngs' --nproc 16

Running trf ...

kneaddata_trf_parallel --input /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.trimmed.fasta --output /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.trimmed.fasta.trf.parameters.2.7.7.80.10.50.500.dat --trf-path /mnt/home/huanggy/miniconda3/envs/snakemake_py36/bin/trf --trf-options '2 7 7 80 10 50 500 -h -ngs' --nproc 16

Decontaminating ...
Running bowtie2 ...

/mnt/project/tools/bowtie2-2.4.4/bowtie2 --threads 16 --very-sensitive --dovetail --phred33 -x /mnt/data/yanhuang/yh_add_hg19_GRCh38_latest_rna -U /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.repeats.removed.fastq --un /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90_yh_add_hg19_GRCh38_latest_rna_bowtie2_clean.fastq --al /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq -S /dev/null

Decontaminating ...
Running bowtie2 ...

/mnt/project/tools/bowtie2-2.4.4/bowtie2 --threads 16 --very-sensitive --dovetail --phred33 -x /mnt/data/yanhuang/yh_add_hg19_GRCh38_latest_rna -U /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.repeats.removed.fastq --un /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97_yh_add_hg19_GRCh38_latest_rna_bowtie2_clean.fastq --al /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq -S /dev/null

Total contaminate sequences in file ( /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq ) : 196369.0
Total reads after removing those found in reference database ( /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90_yh_add_hg19_GRCh38_latest_rna_bowtie2_clean.fastq ): 11066.0
Total reads after merging results from multiple databases ( /mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.fastq ): 11066.0

Final output file created:
/mnt/home/huanggy/project/20211117/result/S200026380_L01_90/S200026380_L01_90.fastq

[Wed Nov 17 20:16:30 2021]
Finished job 1.
1 of 3 steps (33%) done
Total contaminate sequences in file ( /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq ) : 213547.0
Total reads after removing those found in reference database ( /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97_yh_add_hg19_GRCh38_latest_rna_bowtie2_clean.fastq ): 16103.0
Total reads after merging results from multiple databases ( /mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.fastq ): 16103.0

Final output file created:
/mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.fastq


"""
#snakemake_py36  环境bowtie2 报错https://www.jianshu.com/p/499ff4b90b07
##--remove-intermediate-output
rule all:
    ####想要哪一步的结果,都在这里买定义
    input:
        #expand("{outdir}/result/{sample}/",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/{sample}.fastq",outdir=OUTDIR,sample=SAMPLE)

rule kneaddata:
    input:
        #"/mnt/data/NGSDATA/meta_data/S200026380_L01_90.fq.gz",
        #"/mnt/data/NGSDATA/meta_data/S200023932_L01_97.fq.gz"
        #single = expand("/mnt/data/NGSDATA/meta_data/{sample}.fq.gz",sample=SAMPLE)
        single = "/mnt/data/NGSDATA/meta_data/test/{sample}.fq.gz"
    output:
        directory("{outdir}/result/{sample}/"),
        "{outdir}/result/{sample}/{sample}.fastq"

    params:
        thread_num = 16
    shell:
        "{kneaddata} \
        -i {input.single} \
        -o {output[0]}  \
        -v -t {params.thread_num} \
        --trimmomatic /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10' \
        --bowtie2 {bowtie2} --bowtie2-options '--very-sensitive --dovetail' -db {host_yanhaung} \
        --output-prefix {wildcards.sample}"
##### ????????????????????shell 中用sample 这个变量 必须要用wildcards.sample ,单独用{sample} 报错如下：
# RuleException in line 27 of /mnt/home/huanggy/idseq_dag/report/mNGS.snakemake.py:
# NameError: The name 'sample' is unknown in this context. Did you mean 'wildcards.sample'?

rule kneaddata_stat:
    input:
        #expand("{outdir}/result/{sample}/filter",outdir=OUTDIR,sample=SAMPLE)
        "{outdir}/result/{sample}/"
        #The flag 'directory' used in rule kneaddata_stat is only valid for outputs, not inputs.
    output:
        kneaddata_sum = "{outdir}/result/{sample}/kneaddata.sum.txt"
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
        "{outdir}/result/{sample}/{sample}.fastq"
    output:
        "{outdir}/result/{sample}/{sample}.f.fastq"
    shell:
        "{PriceSeqFilter} -a 12 -rnf 90 -log c -f {input} -o {output} -rqf 85 0.98 -lenf 36 && rm -f {input}"
    
