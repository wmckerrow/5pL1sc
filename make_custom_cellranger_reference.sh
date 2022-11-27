hg38_fasta=$1
L1HS_and_L1PA_bed=$2
L1HS_and_L1PA_consensus_fa=$3
hg38_genes_gtf=$4

bedtools maskfasta -fi $hg38_fasta -fo temp.fa -bed $L1HS_and_L1PA_bed
cat temp.fa $L1HS_and_L1PA_consensus_fa > L1HS_L1PA_seperated_hg38.fa
cellranger mkref --genome=L1HS_L1PA_seperated_hg38 --fasta=L1HS_L1PA_seperated_hg38.fa --genes=$hg38_genes_gtf --nthreads=16
