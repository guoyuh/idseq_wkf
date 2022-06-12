config:"config.yaml"
# SAMPLE=["S200026380_L01_90","S200023932_L01_97"]
#OUTDIR="/mnt/home/huanggy/project/20211123"
#OUTDIR="/mnt/home/huanggy/project/20211123"

def get_input_fastqs(wildcards):
    return config["sample_se"][wildcards.sample]

OUTDIR=config["OUTDIR"]

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
        #expand("{outdir}/result/{sample}/filter/{sample}.log",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/filter/{sample}.f.fastq",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.report",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.s.bracken",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.s.bracken.csv",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.report.mpa.rpm.xls",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.krona",outdir=OUTDIR,sample=config["sample_se"]),
        #expand("{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log",outdir=OUTDIR,sample=config["sample_se"])
        expand("{outdir}/result/{sample}/align/{sample}.refseq.sam",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/align/{sample}.bam",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/align/{sample}.sorted.bam",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/align/{sample}.cove.tsv",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/align/{sample}.refseq.cove.xls",outdir=OUTDIR,sample=config["sample_se"]),
        #directory(expand("{outdir}/result/{sample}/filter/spades_out",outdir=OUTDIR,sample=config["sample_se"])),
        expand("{outdir}/result/{sample}/filter/spades_out/contigs.fasta",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/filter/spades_out/{sample}.amr",outdir=OUTDIR,sample=config["sample_se"]),
        #expand("{outdir}/result/stats/{sample}.log",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/filter/kneaddata.sum.txt",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq_output",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq80_output",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq_kreport",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt",outdir=OUTDIR,sample=config["sample_se"]),
        expand("{outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt",outdir=OUTDIR,sample=config["sample_se"])
        #expand("{outdir}/result/stat.xls",outdir=OUTDIR)

        
        

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
        single = get_input_fastqs
    output:
        #directory("{outdir}/result/{sample}/"),
        "{outdir}/result/{sample}/filter/{sample}.fastq",
        "{outdir}/result/{sample}/filter/{sample}.log"

    params:
        thread_num = 16,
        output_dir = "{outdir}/result/{sample}/filter"
    shell:
        "{kneaddata} \
        -i {input.single} \
        -o {params.output_dir}  \
        -v -t {params.thread_num} --remove-intermediate-output \
        --trimmomatic /mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/snakemake_py36/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10' \
        --bowtie2 {bowtie2_path} --bowtie2-options '--very-sensitive --dovetail' -db {host_hg37} \
        --output-prefix {wildcards.sample}"
##### ????????????????????shell 中用sample 这个变量 必须要用wildcards.sample ,单独用{sample} 报错如下：
# RuleException in line 27 of /mnt/home/huanggy/idseq_wkf/report/mNGS.snakemake.py:
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
        "{outdir}/result/{sample}/filter/{sample}.log"
        #The flag 'directory' used in rule kneaddata_stat is only valid for outputs, not inputs.
    output:
        kneaddata_sum = "{outdir}/result/{sample}/filter/kneaddata.sum.txt"
    params:
        indir = "{outdir}/result/{sample}/filter",
        outdir=config["OUTDIR"]
    shell:
        "{kneaddata_read_count_table} --input {params.indir} --output {output.kneaddata_sum}  && rm -f {params.outdir}/result/{wildcards.sample}/filter/{wildcards.sample}_hg37dec_v0.1_bowtie2_contam.fastq"


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
        -w {output[1]} -r 75 -l S \
        && {Bracken} -d {kraken2_db} -i {input} \
        -o {output[2]} \
        -w {output[3]} -r 75 -l G "


rule  bracken_format_rpm:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken",
        "{outdir}/result/{sample}/kraken2/{sample}.g.bracken"
    output:
        ###不是变量,仅作为rule all输出
        "{outdir}/result/{sample}/kraken2/{sample}.s.bracken.csv"
    shell:
        #NameError: The name 'outdir' is unknown in this context. Did you mean 'wildcards.outdir'?
        "python  ~/idseq_wkf/scripts/bracken_reads_rpm.py -w {wildcards.outdir}  -p {wildcards.sample} "

# rule  kreport2mpa:
#     input:
#         "{outdir}/result/{sample}/kraken2/{sample}.report"
#     output:
#         "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt"
#     shell:
#         "python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r {input} -o {output}"

rule kout2kreport:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.output",
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.kreport",
    shell:
        """
        python ~/idseq_wkf/KrakenTools/make_kreport.py -i {input[0]} -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt -o {output[0]}
        """

###阴性结果先跑？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？每一批次指定一个阴控 YXKZ20220416
rule kreport2mpa2rpm_1:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.kreport",
        #"{outdir}/result/{sample}/filter/{sample}.f.fastq"
        "{outdir}/result/{sample}/filter/kneaddata.sum.txt"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.kreport.mpa.txt",
        "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.rpm.xls"
    shell:
        """
        #python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r /mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.kreport -o ~/project/20220408/result/21JS944006/kraken2/21JS944006.kreport.mpa.txt

        python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r {input[0]} -o {output[0]}
        python  ~/idseq_wkf/scripts/kraken_reads_rpm_2.py {input[1]} {output[0]}

            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {input.uniq_kreport} -o {output.uniq_kreport_mpa} 
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {input.uniq80_kreport} -o {output.uniq80_kreport_mpa}
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {input.kneaddata_sum} {input.report_mpa} {output.uniq_kreport_mpa}
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {input.kneaddata_sum} {input.report_mpa} {output.uniq80_kreport_mpa}

            python  ~/idseq_wkf/scripts/kraken_reads_rpm_new.py \
                    /mnt/home/huanggy/project/20220408/result/21JS944006/filter/kneaddata.sum.txt \
                    /mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.kreport.mpa.txt \
                    /mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.kreport.mpa.txt \
                    /mnt/home/huanggy/project/20220408/result/21JS944007/kraken2/21JS944007.report.mpa.txt \
        """



rule  kreport2mpa:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.report",
        #"{outdir}/result/{sample}/filter/{sample}.f.fastq"
        "{outdir}/result/{sample}/filter/kneaddata.sum.txt"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt",
        "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.rpm.xls"
    shell:
        """
        python ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r {input[0]} -o {output[0]}
        python  ~/idseq_wkf/scripts/kraken_reads_rpm_2.py {input[1]} {output[0]}
        """

rule  kreport2krona:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.report"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.krona",
        "{outdir}/result/{sample}/kraken2/{sample}.krona.html"
    shell:
        "python  ~/idseq_wkf/KrakenTools/kreport2krona.py --report-file  {input} -o {output[0]} \
        && perl ~/idseq_wkf/Krona/bin/ktImportText {output[0]}  -o {output[1]}"


rule kraken2_filter:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.output"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.uniq_output",
        "{outdir}/result/{sample}/kraken2/{sample}.uniq80_output"
    shell:
        "python ~/idseq_wkf/scripts/kraken2_uniq_reads_filter.py {input}"


rule make_kreport:
    input:
        "{outdir}/result/{sample}/kraken2/{sample}.uniq_output",
        "{outdir}/result/{sample}/kraken2/{sample}.uniq80_output"
    output:
        "{outdir}/result/{sample}/kraken2/{sample}.uniq_kreport",
        "{outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport"
    shell:
        """
            python ~/idseq_wkf/KrakenTools/make_kreport.py -i {input[0]} -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt -o {output[0]}
            python ~/idseq_wkf/KrakenTools/make_kreport.py -i {input[1]} -t /mnt/data/kraken2_mask/taxonomy/mydb_taxonomy20220408merge.txt -o {output[1]}
        """


rule kreport2mpa2rpm:
    input:
        uniq_kreport = "{outdir}/result/{sample}/kraken2/{sample}.uniq_kreport",
        uniq80_kreport = "{outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport",
        kneaddata_sum = "{outdir}/result/{sample}/filter/kneaddata.sum.txt",
        report_mpa = "{outdir}/result/{sample}/kraken2/{sample}.report.mpa.txt",
    output:
        uniq_kreport_mpa = "{outdir}/result/{sample}/kraken2/{sample}.uniq_kreport.mpa.txt",
        uniq80_kreport_mpa = "{outdir}/result/{sample}/kraken2/{sample}.uniq80_kreport.mpa.txt"
    shell:
        """
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {input.uniq_kreport} -o {output.uniq_kreport_mpa} 
            python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count -r {input.uniq80_kreport} -o {output.uniq80_kreport_mpa}
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {input.kneaddata_sum} {input.report_mpa} {output.uniq_kreport_mpa}
            python  ~/idseq_wkf/scripts/kraken_reads_rpm_3.py {input.kneaddata_sum} {input.report_mpa} {output.uniq80_kreport_mpa}
        """
        



rule clean_reads_bowtie2_bac_fungi_protozoa_viral:
    input:
        "{outdir}/result/{sample}/filter/{sample}.f.fastq",
        "{outdir}/result/{sample}/kraken2/{sample}.krona.html"  ##with --latency-wait
    output:
        "{outdir}/result/{sample}/align/{sample}.clean.refseq.fq.gz",
        "{outdir}/result/{sample}/align/{sample}.refseq.sam"
    params:
        threads=64,
        bac_fungi_protozoa_viral_index="/mnt/data/bowtie2/bacteria_fungi_protozoa_viral"
    log:
        log_file="{outdir}/result/{sample}/align/{sample}.bowtie2.refseq.log"
    shell:
        "{bowtie2_bin} -q --quiet --sensitive --threads {params.threads} -x {params.bac_fungi_protozoa_viral_index} -U {input[0]} --al-gz {output[0]} -S {output[1]} 1>{log.log_file} 2>&1"

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
        "{outdir}/result/{sample}/align/{sample}.sorted.bam"
    shell:
        "{samtools} sort -o {output} {input}"

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
        "{Coverm} genome -m mean covered_bases covered_fraction variance length count rpkm --min-covered-fraction 0 --bam-files {input} --genome-fasta-directory {refsesq_dir} --genome-fasta-extension fna > {output} "

rule Covem_format:
    input:
        "{outdir}/result/{sample}/align/{sample}.cove.tsv"
    output:
        "{outdir}/result/{sample}/align/{sample}.refseq.cove.xls"
    shell:
        "python ~/idseq_wkf/scripts/coverm_format.py {input} /mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20220406.txt /mnt/data/NCBI_Refseq/viral/kraken2_viral_refseq_info "



rule SPAdes:
    input:
        "{outdir}/result/{sample}/filter/{sample}.f.fastq",
        "{outdir}/result/{sample}/align/{sample}.refseq.cove.xls" ###with --latency-wait
    output:
        #directory("{outdir}/result/{sample}/filter/spades_out"),
        # ChildIOException:
        # File/directory is a child to another output:
        # ('/mnt/home/huanggy/project/20211125_test/result/LX2110776/filter/spades_out', SPAdes)
        # ('/mnt/home/huanggy/project/20211125_test/result/LX2110776/filter/spades_out/LX2110776.amr', amrfinder)
        # 同一个rule 里面，output 出现  “{outdir}/result/{sample}/filter/spades_out”   与 "{outdir}/result/{sample}/filter/spades_out/contigs.fasta" 
        ## 对于snakemake 而言 当没有 "{outdir}/result/{sample}/filter/spades_out/contigs.fasta" 时候，自然而然的就创建了 {outdir}/result/{sample}/filter/spades_out/contigs.fasta 文件，当然也创建了上一级目录
        ###{outdir}/result/{sample}/filter/spades_out，此时在output 里面的同时出现，必然会使得Snakemake complains
        "{outdir}/result/{sample}/filter/spades_out/contigs.fasta"

    params:
        out_dir = directory("{outdir}/result/{sample}/filter/spades_out")
    shell:
        #SPAdes
        # if --meta  giveme error Sorry, current version of metaSPAdes can work either with single library (paired-end only) or in hybrid paired-end + (TSLR or PacBio or Nanopore) mode.
        "/mnt/project/tools/SPAdes-3.15.3/bin/spades.py -s {input[0]} -o {params.out_dir}"



rule amrfinder:
    #amrfinder 依赖blastx，blastx add $PATH
    input:
        "{outdir}/result/{sample}/filter/spades_out/contigs.fasta"
    output:
        "{outdir}/result/{sample}/filter/spades_out/{sample}.amr"
    shell:
        "amrfinder -n {input} -o {output}"


# rule kneaddata_log_link:
#     input:
#         "{outdir}/result/{sample}/filter/{sample}.log"

#     output:
#         "{outdir}/result/stats/{sample}.log"

#     shell:
#         "ln -s {input} {output}"

# rule stats:
#     input:
#         expand("{outdir}/result/stats/{sample}.log",outdir=config["OUTDIR"],sample=config["sample_se"])
#     output:
#         "{outdir}/result/stats/kneaddata.stats"

#     params:
#         input_dir = directory("{outdir}/result/stats")
#     shell:
#         "kneaddata_read_count_table --input {params.input_dir} --output {output} && rm -f {outdir}/result/{wildcards.sample}/filter/{wildcards.sample}_L01_yh_add_hg19_GRCh38_latest_rna_bowtie2_contam.fastq "



########################################################################################################
# rule merge_stats:
#python ~/idseq_wkf/scripts/kneaddata_merge.py /mnt/home/huanggy/project/20220325
#     input:
#         #"{outdir}/result/{sample}/filter/{sample}.log"
#         # $ sh ~/idseq_wkf/launch.sh ~/project/20211208_mNGS ~/idseq_wkf/runs/config_20211203.yaml 8
#         # Creating specified working directory /mnt/home/huanggy/project/20211208_mNGS.
#         # Building DAG of jobs...
#         # WildcardError in line 540 of /mnt/home/huanggy/idseq_wkf/runs/mNGS.snakemake.py:
#         # Wildcards in input files cannot be determined from output files:
#         # 'sample'
#         #########expand("{outdir}/result/{sample}/filter/{sample}.log",outdir=OUTDIR,sample=config["sample_se"])
#         expand("{outdir}/result/{sample}/filter/kneaddata.sum.txt",outdir=OUTDIR,sample=config["sample_se"])
#     output:
#         "{outdir}/result/stat.xls"
#     params:
#         outdir=config["OUTDIR"]
#     shell:
#         #"python ~/idseq_wkf/scripts/kneaddata_merge.py /mnt/home/huanggy/project/20211129_mNGS"
#         #RuleException in line 547 of /mnt/home/huanggy/idseq_wkf/runs/mNGS.snakemake.py:
#         #NameError: The name 'outdir' is unknown in this context. Did you mean 'wildcards.outdir'?
#         "python ~/idseq_wkf/scripts/kneaddata_merge.py {params.outdir} "
