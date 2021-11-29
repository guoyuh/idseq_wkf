#!/bin/bash
#conda activate metagenome
if [ $# -lt 3 ];then
    echo "Usage :sh $0 <outdir> <sample> <fq> "
    exit 1
fi

#which fastqc
#which multiqc
#which trimmomatic
#which bowtie2
#which seqtk
#which kraken2
#which bracken
#which samtools
#which coverm
gsnapl=/mnt/project/tools/gmap-2021-07-23/bin/gsnapl
hg37_genome_index=/mnt/db/kneaddata/human_genome/hg37dec_v0.1
bac_fungi_viral_index=/mnt/data/bowtie2/bac_fungi_viral


outdir=$1
sample_name=$2
fq=$3

if [ ! -d $outdir ];then
  mkdir -p $outdir
fi
#mkdir -p ${outdir}/result/{filter,qc,kraken2,align}
cd $outdir

#初始fastqc
#time fastqc seq/*.gz -t 6
##############trimatic##################
mkdir -p  ${outdir}/result/filter
java -Xmx500m -jar /mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/trimmomatic.jar \
SE -threads 16 ${fq} \
${outdir}/result/filter/${sample_name}.trimmed.fq.gz \
ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
echo "run status $?"
echo "trimmomatic finished"

################bowtie2 去宿主########
bowtie2 -q --quiet --sensitive --threads 64 \
-x ${hg37_genome_index} \
-U ${outdir}/result/filter/${sample_name}.trimmed.fq.gz \
--un-gz ${outdir}/result/filter/${sample_name}.clean.fq.gz \
1>${outdir}/result/filter/${sample_name}.bowtie2.log 2>&1
##-un: unpaired unaligned
##--un-conc: paired unaligned
echo "run status $?"
echo "bowtie2 remove hostgenome finished"

##test##
#bowtie2 -q --quiet --sensitive --threads 64 \
#-x /mnt/db/kneaddata/human_genome/hg37dec_v0.1 \
#-U /mnt/home/huanggy/project/20210926/result/filter/S200023771_L01_101.trimmed.fq.gz \
#--un-gz /mnt/home/huanggy/project/20210926/result/S200023771_L01_101.clean.fq.gz 1>/mnt/home/huanggy/project/20210926/result/bowtie2.log 2>&1


Initial number of reads ( /mnt/home/huanggy/project/20210926_test/result/temp/qc/decompressed_z5kjor0i_S200026521_L01_102.fq ): 41414585.0
Running Trimmomatic ...

java -Xmx500m -jar /mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/trimmomatic.jar SE -threads 3 -phred33 /mnt/home/huanggy/project/20210926_test/result/temp/qc/decompressed_z5kjor0i_S200026521_L01_102.fq /mnt/home/huanggy/project/20210926_test/result/temp/qc/S200026521_L01_102_kneaddata.trimmed.fastq SLIDINGWINDOW:4:15 MINLEN:36

Total reads after trimming ( /mnt/home/huanggy/project/20210926_test/result/temp/qc/S200026521_L01_102_kneaddata.trimmed.fastq ): 40970101.0
Decontaminating ...
Running bowtie2 ...

/mnt/home/huanggy/miniconda3/envs/metagenome/bin/bowtie2 --threads 3 --very-sensitive --dovetail --phred33 -x /mnt/db/kneaddata/human_genome/hg37dec_v0.1 -U /mnt/home/huanggy/project/20210926_test/result/temp/qc/S200026521_L01_102_kneaddata.trimmed.fastq --un /mnt/home/huanggy/project/20210926_test/result/temp/qc/S200026521_L01_102_kneaddata_hg37dec_v0.1_bowtie2_clean.fastq --al /mnt/home/huanggy/project/20210926_test/result/temp/qc/S200026521_L01_102_kneaddata_hg37dec_v0.1_bowtie2_contam.fastq -S /dev/null


kneaddata -i ${fq} \
-o temp/qc -v -t 32 --remove-intermediate-output \
--trimmomatic /mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1 \
--trimmomatic-options "SLIDINGWINDOW:4:15 MINLEN:36 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10" \
--bowtie2-options "--very-sensitive --dovetail" -db /mnt/db/kneaddata/human_genome





###kraken2 & bracken classificatioon and abuandance   #####
mkdir -p  ${outdir}/result/kraken2
kraken2 --db /mnt/data/kraken2_mask  ${outdir}/result/filter/${sample_name}.clean.fq.gz \
--threads 16 \
--report ${outdir}/result/kraken2/${sample_name}.report \
--output ${outdir}/result/kraken2/${sample_name}.output

###种水平丰度估计
/mnt/project/tools/Bracken/bracken -d /mnt/data/kraken2_mask -i ${outdir}/result/kraken2/${sample_name}.report \
-o ${outdir}/result/kraken2/${sample_name}.s.bracken \
-w ${outdir}/result/kraken2/${sample_name}.s.bracken.report -r 150 -l S 


###属水平丰度估计
/mnt/project/tools/Bracken/bracken -d /mnt/data/kraken2_mask -i ${outdir}/result/kraken2/${sample_name}.report \
-o ${outdir}/result/kraken2/${sample_name}.g.bracken \
-w ${outdir}/result/kraken2/${sample_name}.g.bracken.report -r 150 -l G  

echo "run status $?"
echo "kraken2 & bracken finished"

####质控
mkdir -p  ${outdir}/result/qc 
fastqc ${outdir}/result/filter/*.gz -t 6
#multiqc汇总#
multiqc -d ${outdir}/result/filter/ -o ${outdir}/result/qc




##fq2fa & blast
mkdir -p ${outdir}/result/align

/mnt/project/tools/seqtk/seqtk seq -a ${outdir}/result/filter/${sample_name}.clean.fq.gz > ${outdir}/result/filter/${sample_name}.clean.fa


#######基于clinrefseq mapping ########
$gsnapl \
-A m8 --batch=0 --use-shared-memory=0  --npaths=100 --ordered -t 48 --max-mismatches=40 \
-D /mnt/data/gsnap -d refseq_k16 ${outdir}/result/filter/${sample_name}.clean.fa > ${outdir}/result/align/${sample_name}.out
echo "run status $?"
echo "gsnapl mapping to blast finished"

python  ~/idseq_dag/m8_6.py ${outdir}/result/align/${sample_name}.out ${outdir}/result/align ${sample_name}

###extract
## gsnapl 比对报错？？？？
#$gsnapl \
#-A sam --batch=0 --use-shared-memory=0 --npaths=1 --ordered -t 48 --max-mismatches=40 \
#-D /mnt/data/gsnap -d refseq_k16 ${outdir}/result/filter/${sample_name}.clean.fa \
#> ${outdir}/result/align/${sample_name}.sam
#echo "run status $?"
#echo "gsnapl mapping to bam finished"


####bowtie2 构建索引(5h)  bowtie2-build  --large-index  -f /mnt/data/kraken2_db/bac_fungi_viral.fna /mnt/data/bowtie2/bac_fungi_viral --threads 64
##bowtie2 比对
bowtie2 -q --quiet --sensitive --threads 64 \
-x $bac_fungi_viral_index \
-U ${outdir}/result/filter/${sample_name}.clean.fq.gz \
--al-gz ${outdir}/result/align/${sample_name}.clean.refseq.fq.gz \
-S ${outdir}/result/align/${sample_name}.refseq.sam \
1>${outdir}/result/align/${sample_name}.bowtie2.refseq.log 2>&1
echo "run status $?"
echo "gsnapl mapping to bam finished"


samtools view -bS ${outdir}/result/align/${sample_name}.refseq.sam > ${outdir}/result/align/${sample_name}.bam  
samtools sort ${outdir}/result/align/${sample_name}.bam  ${outdir}/result/align/${sample_name}.sorted
samtools index ${outdir}/result/align/${sample_name}.sorted.bam && rm -f ${outdir}/result/align/${sample_name}.refseq.sam  ${outdir}/result/align/${sample_name}.bam 

coverm genome -m covered_fraction variance rpkm --bam-files ${outdir}/result/align/${sample_name}.sorted.bam --genome-fasta-directory /mnt/data/NCBI_Refseq/bac_fungi_viral_protozoa  > ${outdir}/result/align/${sample_name}.cove.tsv
awk -F"\t" '{if ($2!=0) print $0}' ${outdir}/result/align/${sample_name}.cove.tsv > ${outdir}/result/align/${sample_name}.f.cove.tsv  
#&& rm -f ${outdir}/result/align/${sample_name}.cove.tsv
echo "run status $?"
echo "coverm finished"


