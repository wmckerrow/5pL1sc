genome_fasta=$1
reps_to_mask_bed=$2
rep_consensus_fa=$3
genes_gtf=$4
output_name=$5

bedtools maskfasta -fi $genome_fasta -fo temp.fa -bed $reps_to_mask_bed
cat temp.fa $rep_consensus_fa > $output_name.fa
cellranger mkref --genome=$output_name --fasta=$output_name.fa --genes=$genes_gtf --nthreads=16
rm temp.fa
