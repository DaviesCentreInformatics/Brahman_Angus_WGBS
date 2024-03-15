#!/bin/bash -l
fasta=$1
divergence=$2
num_seq_to_make=$3
name_prefix=$4

chromosome=$(basename $fasta .fa)

chr_len=$(tail -n +2 $fasta | tr -d '\n' | wc -m)
num_snps_to_mutate=$(echo "$chr_len * $divergence" | bc | cut -d '.' -f 1)

snpmutator -r 12 \
	-n 1 \
	-s $num_snps_to_mutate \
	-o metadata/${name_prefix}_${chromosome}.summary \
	-v metadata/${name_prefix}_${chromosome}.vcf \
	-M metadata/${name_prefix}_${chromosome}.metrics \
	$fasta
