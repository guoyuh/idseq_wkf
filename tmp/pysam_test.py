import pysam


samfile = pysam.AlignmentFile("/mnt/home/huanggy/project/20211101/result/LX2110967/align/LX2110967.sorted.bam", "rb")

samfile2 = pysam.AlignmentFile("/mnt/home/huanggy/project/20211101/result/LX2110582/align/LX2110582.sorted.bam", "rb")

for read in samfile2.fetch('NC_007795.1', 1266392, 1729170):
    print(read)

# for pileupcolumn in samfile2.pileup('NC_007795.1', 1266392, 1729170):
#     print ("\ncoverage at base %s = %s" %
#            (pileupcolumn.pos, pileupcolumn.n))
#     for pileupread in pileupcolumn.pileups:
#         if not pileupread.is_del and not pileupread.is_refskip:
#             # query position is None if is_del or is_refskip is set.
#             print ('\tbase in read %s = %s' %
#                   (pileupread.alignment.query_name,
#                    pileupread.alignment.query_sequence[pileupread.query_position]))



print(samfile2.count_coverage('NC_007795.1',start=0,stop=1729170))