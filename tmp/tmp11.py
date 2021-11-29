

# with open("/mnt/home/huanggy/project/20210927_test/result/kraken2/LX2110919.s.bracken",'r') as fh:
#     next(fh)

#     for line in fh:
#         #print(line.split("\t"))
#         name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads = line.strip().split("\t")
#         print("======>", name,taxonomy_id, taxonomy_lvl, kraken_assigned_reads, new_est_reads, fraction_total_reads)



# from Bio import SeqIO

# file_path="/mnt/data/NCBI_Refseq/bacteria/my_bac_refseq_back/GCF_016727585.1_ASM1672758v1_genomic.fna"
# with open(file_path,'r') as handle:
#     record = SeqIO.parse(handle,"fasta")
#     for (ind,seq_record) in enumerate(record):
#         print(seq_record.id ,"------>",len(seq_record.seq))


import pandas as pd
import numpy as np

def assembleFile2strain(assemble_summary):
    assemble_genome2organism_name = {}
    with open(assemble_summary) as fh:
        for line in fh:
            if line.startswith("assembly_accession"):
                continue
            if not line:
                break
            else:
                organism_name = line.strip().split("\t")[7]
                ftp_path = line.strip().split('\t')[19].split('/')[-1] + "_genomic"
                #print(organism_name + "------>" + ftp_path)
                if ftp_path not in organism_name:
                    assemble_genome2organism_name[ftp_path] = organism_name
                else:
                    print("assemble_genome_name 有重复")

    return assemble_genome2organism_name


def run(coverm_file,assemble_summary):
    assemble_genome2organism_name = assembleFile2strain(assemble_summary)
    #print(assemble_genome2organism_name)
    tb = pd.read_table(coverm_file,sep='\t',header=0)
    tb.columns = ['Genome','Covered Fraction','Variance','Length','my_refseq_Read_Count','RPKM']
    # for i in range(len(tb)):
    #     print(tb.loc[i, 'Genome'])
    #     #if assemble_genome2organism_name.get(tb.loc[i, 'Genome']):
    #     tb.loc[i, 'Genome'] = assemble_genome2organism_name.get(tb.loc[i, 'Genome'])

    tb['Genome'] = tb['Genome'].map(assemble_genome2organism_name,na_action = "abc")
    #tb = tb.reindex(index=tb['Genome'])
    tmp = tb[['Covered Fraction','Variance','Length','my_refseq_Read_Count','RPKM']]]
    print(tb)
    



run(coverm_file = '/mnt/home/huanggy/project/20211101/result/LX2110582/align/LX2110582.cove.tsv',assemble_summary = '/mnt/data/NCBI_Refseq/my_assembly_summary_refseq.info20211101.txt')