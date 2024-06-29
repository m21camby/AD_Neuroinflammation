#!/bin/bash

# written by Seung Kim and Marcel Schilling
# March 2 2020

###########
# version
###########
# sambamba 0.6.8
# bedtools v2.27.1
# bedGraphToBigWig v4

#####################
# method description
#####################
# Sambamba (0.6.8) was used to sort BAM file produced from 10X cellranger count
# We extracted only primary aligment reads from Sorted BAM file
# We created bedgraph file using bedtools (2.27.1) and splited file by strand specific
# We created BigWig file using bedGraphToBigWig from Kent utils

############
# input
############
# FILE is BAM file from 10X cellranger count
# PREFIX is sample ID to add as prefix

FILE=$1
PREFIX=$2

#################################
# 1. Convert to BAM file and sort
#################################
sambamba view --format=bam --compression-level=9 --nthreads=12 '/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/'$FILE | \
sambamba sort --compression-level=9 --nthreads=12 --memory-limit=10G --tmpdir='/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/BigWig_Files/' --out=/dev/stdout /dev/stdin >  $PREFIX'.sorted.bam'


####################################
# 2. Extract only primary alignment
####################################
sambamba view --format=bam --compression-level=9 --nthreads=12 --filter='not secondary_alignment' /dev/stdin < $PREFIX'.sorted.bam' > $PREFIX'_unimappers.sorted.bam'


#########################
# 3. convert to bedgraph 
#########################

# split whole BAM file
bedtools genomecov -ibam - -bg -split -strand '+' -scale $(sambamba flagstat /dev/stdin < $PREFIX'.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'.coverage.plus.rpm.bedgraph'


bedtools genomecov -ibam - -bg -split -strand '-' -scale $(sambamba flagstat /dev/stdin < $PREFIX'.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'.coverage.minus.rpm.bedgraph'


# split unimapper BAM file
bedtools genomecov -ibam - -bg -split -strand '+' -scale $(sambamba flagstat /dev/stdin < $PREFIX'_unimappers.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'_unimappers.coverage.plus.rpm.bedgraph'


bedtools genomecov -ibam - -bg -split -strand '-' -scale $(sambamba flagstat /dev/stdin < $PREFIX'_unimappers.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'_unimappers.coverage.minus.rpm.bedgraph'



####################
# 4. create BigWig
####################
./bedGraphToBigWig $PREFIX'_unimappers.coverage.plus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size $PREFIX'_unimappers.coverage.plus.rpm_test.bw'

./bedGraphToBigWig $PREFIX'_unimappers.coverage.minus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size $PREFIX'_unimappers.coverage.minus.rpm_test.bw'


