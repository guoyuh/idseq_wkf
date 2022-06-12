#!/bin/bash
#conda activate metagenome
#####################修改########################
#1,删除GSNAP 比对及其解析
#2，covrm  修改参数min-covered-fraction == 0
##############################################
#which fastqc
#which multiqc
#which trimmomatic
#which bowtie2
#which PriceSeqFilter
#which seqtk
#which kraken2
#which bracken
#which samtools
#which coverm
#gsnapl=/mnt/project/tools/gmap-2021-07-23/bin/gsnapl
#hg37_genome_index=/mnt/db/kneaddata/human_genome/hg37dec_v0.1
bac_fungi_protozoa_viral_index=/mnt/data/bowtie2/bacteria_fungi_protozoa_viral
refsesq_dir=/mnt/data/NCBI_Refseq/bac_fungi_protozoa_viral
host_yanhaung=/mnt/data/yanhuang ### 由炎黄1号 + somatic-b37_Homo_sapiens_assembly19.fasta + GRCh38_latest_rna
host_hg37=/mnt/db/kneaddata/human_genome ### 由hg37 + human_contamination  组成from 官网
kraken2_db=/mnt/data/kraken2_mask
outdir=$1
# sample_name=$2
# fq=$3
config=$2

if [ $# -lt 2 ];then
    echo "Usage :sh $0 <outdir> config> "
    exit 1
fi

if [ ! -d $outdir ];then
  mkdir -p $outdir
fi
#mkdir -p ${outdir}/result/{filter,qc,kraken2,align}

##############trimmomatic(trim adapters sliding window ,MINLEN)  bowtie2 去宿主##################
mkdir -p ${outdir}/result/filter
cat $config | while read id ;
do
	arr=($id)
	sample_name=${arr[0]};fq=${arr[1]}
	echo "$sample_name $fq"
	kneaddata -i ${fq} \
	-o ${outdir}/result/filter  -v -t 8 --remove-intermediate-output \
	--trimmomatic /mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1 \
	--trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:40 ILLUMINACLIP:/mnt/home/huanggy/miniconda3/envs/metagenome/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10" \
	--bowtie2-options "--very-sensitive --dovetail" -db /mnt/db/kneaddata/human_genome \
	--output-prefix ${sample_name} 
done
echo "kneaddata:run status $?"
echo "kneaddata finished"

##################kneaddata 统计log#############################################################
kneaddata_read_count_table --input ${outdir}/result/filter --output ${outdir}/result/filter/kneaddata.sum.txt
rm -f ${outdir}/result/filter/${sample_name}_hg37dec*.fastq  && rm -f ${outdir}/result/filter/${sample_name}.trimmed*.fastq
echo "trimmomatic bowtie2 remove hostgenome:run status $?"
echo "trimmomatic bowtie2 remove hostgenome finished"

##############进一步去除低质量 ,短序列##########################################################
cat $config | while read id;
do
    arr=($id)
    sample_name=${arr[0]}
    PriceSeqFilter -a 12 -rnf 90 -log c -f ${outdir}/result/filter/${sample_name}.fastq  \
    -o ${outdir}/result/filter/${sample_name}.f.fastq -rqf 85 0.98 -lenf 36 
    #rm -f ${outdir}/result/filter/${sample_name}.fastq
done 

#############kraken2 & bracken classificatioon and abuandance   #########################
mkdir -p  ${outdir}/result/kraken2
ls ${outdir}/result/filter/*.f.fastq | while read id ;
do
    sp=$(basename -s .f.fastq $id)
    kraken2 --db /mnt/data/kraken2_mask  $id \
    --threads 4 \
    --report ${outdir}/result/kraken2/$sp.report \
    --output ${outdir}/result/kraken2/$sp.output 
done 

###########种水平丰度估计#########
ls ${outdir}/result/kraken2/*.report | while read id ;
do
    sp=$(basename -s .report $id)
    /mnt/project/tools/Bracken/bracken -d /mnt/data/kraken2_mask -i $id \
    -o ${outdir}/result/kraken2/$sp.s.bracken \
    -w ${outdir}/result/kraken2/$sp.s.bracken.report -r 75 -l S 
done


#########属水平丰度估计############
ls ${outdir}/result/kraken2/*.report | while read id ;
do
    sp=$(basename -s .report $id)
    /mnt/project/tools/Bracken/bracken -d /mnt/data/kraken2_mask -i $id \
    -o ${outdir}/result/kraken2/$sp.g.bracken \
    -w ${outdir}/result/kraken2/$sp.g.bracken.report -r 75 -l G
done
echo "kraken2 & bracken:run status $?"
echo "kraken2 & bracken finished"


##bracken format and kraken rpm
#python  ~/idseq_wkf/scripts/bracken_reads_rpm.py -w ${outdir}  -p ${sample_name}
#kraken2 report2mpa
#python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count  -r ${outdir}/result/${sample_name}/kraken2/${sample_name}.report -o ${outdir}/result/${sample_name}/kraken2/${sample_name}.report.mpa.txt
#python  ~/idseq_wkf/scripts/kraken_reads_rpm.py ${outdir}/result/${sample_name}/filter/${sample_name}.f.fastq  ${outdir}/result/${sample_name}/kraken2/${sample_name}.report.mpa.txt
#python  ~/idseq_wkf/scripts/kraken_reads_rpm_2.py ${outdir}/result/${sample_name}/filter/kneaddata.sum.txt  ${outdir}/result/${sample_name}/kraken2/${sample_name}.report.mpa.txt
######################kraken2 report2mpa   && report.mpa  标准化#######################
ls ${outdir}/result/kraken2/*.report | while read id ;
do
    sp=$(basename -s .report $id)
    python  ~/idseq_wkf/KrakenTools/kreport2mpa.py --no-intermediate-ranks --read_count \
    -r $id \
    -o ${outdir}/result/kraken2/$sp.report.mpa.txt
    #&& python  ~/idseq_wkf/scripts/kraken_reads_rpm_2.py ${outdir}/result/filter/kneaddata.sum.txt  ${outdir}/result/${sample_name}/kraken2/${sample_name}.report.mpa.txt
    #kraken2 report2krona 
    python  ~/idseq_wkf/KrakenTools/kreport2krona.py --report-file  $id  -o ${outdir}/result/kraken2/$sp.krona
    ##Krona 可视化
    perl ~/idseq_wkf/Krona/bin/ktImportText ${outdir}/result/kraken2/$sp.krona  -o ${outdir}/result/kraken2/$sp.krona.html
done
echo "kraken2 report2mpa status $?"
echo "kraken2 report2mpa finished"


##bowtie2 比对
mkdir -p ${outdir}/result/align
ls ${outdir}/result/filter/*.f.fastq | while read id ;
do
    sp=$(basename -s .f.fastq $id)
    bowtie2 -q --quiet --sensitive --threads 8 \
    -x ${bac_fungi_protozoa_viral_index} \
    -U $id \
    --al-gz ${outdir}/result/align/${sp}.clean.refseq.fq.gz \
    -S ${outdir}/result/align/${sp}.refseq.sam \
    1>${outdir}/result/align/${sp}.bowtie2.refseq.log 2>&1
done
echo "bowtie2 mapping clean fastq to refseq:run status $?"
echo "bowtie2 mapping to bam finished"

cat $config | while read id;
do
    arr=($id)
    sample_name=${arr[0]}
    samtools view -bS ${outdir}/result/align/${sample_name}.refseq.sam > ${outdir}/result/align/${sample_name}.bam  
    samtools sort ${outdir}/result/align/${sample_name}.bam  ${outdir}/result/align/${sample_name}.sorted
    samtools index ${outdir}/result/align/${sample_name}.sorted.bam && rm -f ${outdir}/result/align/${sample_name}.refseq.sam  ${outdir}/result/align/${sample_name}.bam
done


cat $config | while read id;
do
    arr=($id)
    sample_name=${arr[0]}
    coverm genome -m covered_fraction variance length count rpkm --min-covered-fraction 0 \
    --bam-files ${outdir}/result/align/${sample_name}.sorted.bam --genome-fasta-directory ${refsesq_dir} --genome-fasta-extension fna > ${outdir}/result/align/${sample_name}.cove.tsv

    ##asssemble2 species
    python ~/idseq_wkf/scripts/coverm_format.py \
    ${outdir}/result/align/${sample_name}.cove.tsv \
    /mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20211101.txt /mnt/data/NCBI_Refseq/viral/kraken2_viral_refseq_info
done
echo "coverm:run status $?"
echo "coverm finished"


#SPAdes assemble
