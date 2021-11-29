SAMPLE=["S200026380_L01_90","S200023932_L01_97"]
ind=["90","97"]
OUTDIR="/mnt/home/huanggy/project/20211120"

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
"""
#snakemake_py36  环境bowtie2 报错https://www.jianshu.com/p/499ff4b90b07
##--remove-intermediate-output
rule all:
    ####想要哪一步的结果,都在这里买定义
    input:
        #expand("{outdir}/result/{sample}/",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/filter/{sample}.fastq",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/filter/{sample}.f.fastq",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/kraken2/{sample}.report",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/kraken2/{sample}.s.bracken",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/kraken2/{sample}.s.bracken.csv",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/kraken2/{sample}.krona",outdir=OUTDIR,sample=SAMPLE),
        #expand("{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log",outdir=OUTDIR,sample=SAMPLE)
        expand("{outdir}/result/{sample}/align/{sample}.refseq.sam",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/align/{sample}.bam",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/align/{sample}.sorted.bam",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/align/{sample}.cove.tsv",outdir=OUTDIR,sample=SAMPLE),
        expand("{outdir}/result/{sample}/align/{sample}.refseq.cove.tsv",outdir=OUTDIR,sample=SAMPLE),



"""
https://stackoverflow.com/questions/58808354/childioexception-error-in-snake-make-after-running-flye
Creating specified working directory ag.
Building DAG of jobs...
ChildIOException:
File/directory is a child to another output:
('/mnt/home/huanggy/project/20211117/result/S200023932_L01_97', kneaddata)
('/mnt/home/huanggy/project/20211117/result/S200023932_L01_97/S200023932_L01_97.f.fastq', PriceSeqFilter)
"""


rule kneaddata:
    input:
        #"/mnt/data/NGSDATA/meta_data/S200026380_L01_90.fq.gz",
        #"/mnt/data/NGSDATA/meta_data/S200023932_L01_97.fq.gz"
        #single = expand("/mnt/data/NGSDATA/meta_data/{sample}.fq.gz",sample=SAMPLE)
        single = "/mnt/data/NGSDATA/meta_data/test/{sample}.fq.gz"
    output:
        #directory("{outdir}/result/{sample}/"),
        "{outdir}/result/{sample}/filter/{sample}.fastq"

    params:
        thread_num = 16,
        output_dir = "{outdir}/result/{sample}/filter/"
    shell:
        "{kneaddata} \
        -i {input.single} \
        -o {params.output_dir}  \
        -v -t {params.thread_num} \
        --trimmomatic /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10' \
        --bowtie2 {bowtie2_path} --bowtie2-options '--very-sensitive --dovetail' -db {host_yanhaung} \
        --output-prefix {wildcards.sample}"
##### ????????????????????shell 中用sample 这个变量 必须要用wildcards.sample ,单独用{sample} 报错如下：
# RuleException in line 27 of /mnt/home/huanggy/idseq_dag/report/mNGS.snakemake.py:
# NameError: The name 'sample' is unknown in this context. Did you mean 'wildcards.sample'?




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
        "{PriceSeqFilter} -a 12 -rnf 90 -log c -f {input} -o {output} -rqf 85 0.98 -lenf 36"
    

rule kneaddata_stat:
    input:
        #expand("{outdir}/result/{sample}/filter",outdir=OUTDIR,sample=SAMPLE)
        "{outdir}/result/{sample}/filter"
        #The flag 'directory' used in rule kneaddata_stat is only valid for outputs, not inputs.
    output:
        kneaddata_sum = "{outdir}/result/{sample}/filter/kneaddata.sum.txt"
    shell:
        "{kneaddata_read_count_table} --input {input} --output {output.kneaddata_sum}"



###snakemake -s ~/idseq_dag/report/mNGS.snakemake.py  --dag --cores 1 | dot -Tpdf > test.pdf


rule kraken2_classify:
    input:
        "{outdir}/result/{sample}/filter/{sample}.f.fastq"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.report",
        "{outdir}/result/{sample}/kraken2/{sample}.output"
    params:
        thread_num = 16
    shell:
        "{KRAKEN} --db {kraken2_db} {input} \
        --threads {params.thread_num} \
        --report {output[0]} \
        --output {output[1]}"

rule bracken_Abundance_estimation:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.report"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken",
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken.report",
        "{outdir}/result/{sample}/kraken2/{sample}.g.bracken",
        "{outdir}/result/{sample}/kraken2/{sample}.g.bracken.report"

    shell:
        "{Bracken} -d {kraken2_db} -i {input} \
        -o {output[0]} \
        -w {output[1]} -r 50 -l S \
        && {Bracken} -d {kraken2_db} -i {input} \
        -o {output[2]} \
        -w {output[3]} -r 50 -l G "


rule  bracken_format_rpm:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken",
        "{outdir}/result/{sample}/kraken2/{sample}.g.bracken"
    output:
        ###不是变量,仅作为rule all输出
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken.csv"
    shell:
        #NameError: The name 'outdir' is unknown in this context. Did you mean 'wildcards.outdir'?
        "python  ~/idseq_dag/kraken_reads_rpm.py -w {wildcards.outdir}  -p {wildcards.sample}"

rule  kreport2mpa:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.report"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt"
    shell:
        "python ~/idseq_dag/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r {input} -o {output}"

rule  kreport2krona:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.report"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.krona",
        "{outdir}/result/{sample}/kraken2/{sample}.krona.html"
    shell:
        "python  ~/idseq_dag/KrakenTools/kreport2krona.py --report-file  {input} -o {output[0]} \
        && perl ~/idseq_dag/Krona/bin/ktImportText {output[0]}  -o {output[1]}"


# rule clean_reads_bowtie2_bac_fungi_protozoa_viral:
#     input:
#         expand("{outdir}/result/{sample}/filter/{sample}.f.fastq",outdir=OUTDIR,sample=SAMPLE)
#     output:
#         expand("{outdir}/result/{sample}/align/{sample}.clean.refseq.fq.gz",outdir=OUTDIR,sample=SAMPLE),
#         expand("{outdir}/result/{sample}/align/{sample}.refseq.sam",outdir=OUTDIR,sample=SAMPLE)
#     params:
#         bowtie2_bin="{bowtie2}/bowtie2",
#         threads=64,
#         bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral"
#     log:
#         log_file"{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log",
#     shell:
#         "echo /mnt/project/tools/bowtie2-2.4.4/bowtie2 -q --quiet --sensitive --threads {params.threads} -x {paramas.bac_fungi_protozoa_viral_index} -U {input} --al-gz {output[0]} -S {output[1]} 1>{log.log_file} 2>&1"


# rule clean_reads_bowtie2_bac_fungi_protozoa_viral:
#     input:
#         "{outdir}/result/{sample}/filter/{sample}.f.fastq"
#     output:
#         "{outdir}/result/{sample}/align/{sample}.clean.refseq.fq.gz",
#         "{outdir}/result/{sample}/align/{sample}.refseq.sam"
#     params:
#         threads=64,
#         bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral",
#     # log:
#     #     log_file="{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log",
#     shell:
#         """
#         bowtie2 -q --quiet --sensitive --threads 16 -x {params.bac_fungi_protozoa_viral_index} -U {input} --al-gz {output[0]} -S {output[1]} 
#         """


rule clean_reads_bowtie2_bac_fungi_protozoa_viral:
    input:
        "{outdir}/result/{sample}/filter/{sample}.f.fastq"
    output:
        "{outdir}/result/{sample}/align/{sample}.clean.refseq.fq.gz",
        "{outdir}/result/{sample}/align/{sample}.refseq.sam"
    params:
        threads=64,
        bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral"
    log:
        log_file="{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log"
    shell:
        "{bowtie2_bin} -q --quiet --sensitive --threads {params.threads} -x {params.bac_fungi_protozoa_viral_index} -U {input} --al-gz {output[0]} -S {output[1]} 1>{log.log_file} 2>&1"

rule sam2bam:
    input:
        "{outdir}/result/{sample}/align/{sample}.refseq.sam"
    output:
        "{outdir}/result/{sample}/align/{sample}.bam"
    shell:
        "{samtools} view -bS {input} > {output}"


rule samtools_sort:
    input:
        "{outdir}/result/{sample}/align/{sample}.bam"
    output:
        "{outdir}/result/{sample}/align/{sample}.sorted",
        "{outdir}/result/{sample}/align/{sample}.sorted.bam"
    shell:
        "{samtools} sort -o {output[1]} {input}"

# rule bam_index:
#     input:
#         prefix = "{outdir}/result/{sample}/align/{sample}.sorted",
#         refseq_sam = "{outdir}/result/{sample}/align/{sample}.refseq.sam",
#         bam="{outdir}/result/{sample}/align/{sample}.bam"
#     output:
#         ###不做output ，
#         "{outdir}/result/{sample}/align/{sample}.sorted.bam"
#     shell:
#         "samtools index {input.prefix}.bam && rm -f {input.refseq_sam} {input.bam}"

rule Covem:
    input:
        "{outdir}/result/{sample}/align/{sample}.sorted.bam"
    output:
        "{outdir}/result/{sample}/align/{sample}.cove.tsv"
    shell:
        "{Coverm} genome -m covered_fraction variance length count rpkm --min-covered-fraction 0 \
            --bam-files {input} \
            --genome-fasta-directory {refsesq_dir} --genome-fasta-extension fna > {output} "

rule Covem_format:
    input:
        "{outdir}/result/{sample}/align/{sample}.cove.tsv"
    output:
        "{outdir}/result/{sample}/align/{sample}.refseq.cove.tsv"
    shell:
        "python ~/idseq_dag/coverm_format.py {input} /mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20211101.txt /mnt/data/NCBI_Refseq/viral/kraken2_viral_refseq_info "