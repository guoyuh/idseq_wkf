set -vex
hostname
date
cd /mnt/home/huanggy/idseq_wkf/runs
echo task start
makeblastdb -in /mnt/home/huanggy/idseq_wkf/runs/db.fasta -dbtype nucl
touch /mnt/home/huanggy/idseq_wkf/runs/make_db_done
echo task done
date
