#!/bin/bash

# written by Marcel and modified by Seung
# Oct 28 2019
#######################################################################################
# Error from original file (BigWig_Normalized_Marcel.sh)
# This file is for testing below error from original file
# ./BigWig_Normalized_Marcel.sh: line 35: syntax error near unexpected token `done'
# ./BigWig_Normalized_Marcel.sh: line 35: `done'
# for loop of STAR mapping using bash create above errors 
#######################################################################################

#####################################################
# second error message during sambamba sort:
# not enough data in stream
# Failed to open BAM file -
# needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
###################################################### 

###########
# version
###########
# sambamba 0.6.8
# bedtools v2.27.1

FILE=$1
PREFIX=$2

#################################
# 1. Convert to BAM file and sort
#################################

#sambamba sort --compression-level=9 --nthreads=12 $FILE > $PREFIX'.sorted.bam'

####################################
# old way
###################################
sambamba view --format=bam --compression-level=9 --nthreads=12 '/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/'$FILE | \
sambamba sort --compression-level=9 --nthreads=12 --memory-limit=10G --tmpdir='/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/BigWig_Files/' --out=/dev/stdout /dev/stdin >  $PREFIX'.sorted.bam'


####################################
# 2. Extract only primary alignment
####################################
sambamba view --format=bam --compression-level=9 --nthreads=12 --filter='not secondary_alignment' /dev/stdin < $PREFIX'.sorted.bam' > $PREFIX'_unimappers.sorted.bam'



###################################
# 3. bedgraph 
###################################

# whole BAM file
bedtools genomecov -ibam - -bg -split -strand '+' -scale $(sambamba flagstat /dev/stdin < $PREFIX'.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'.coverage.plus.rpm.bedgraph'


bedtools genomecov -ibam - -bg -split -strand '-' -scale $(sambamba flagstat /dev/stdin < $PREFIX'.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'.coverage.minus.rpm.bedgraph'


# unimapper BAM file
bedtools genomecov -ibam - -bg -split -strand '+' -scale $(sambamba flagstat /dev/stdin < $PREFIX'_unimappers.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'_unimappers.coverage.plus.rpm.bedgraph'


bedtools genomecov -ibam - -bg -split -strand '-' -scale $(sambamba flagstat /dev/stdin < $PREFIX'_unimappers.sorted.bam' | gawk '$4=="mapped"&&$0=10^6/$1') < $PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > $PREFIX'_unimappers.coverage.minus.rpm.bedgraph'



#################
# 4. BigWig
#################
./bedGraphToBigWig $PREFIX'_unimappers.coverage.plus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size $PREFIX'_unimappers.coverage.plus.rpm_test.bw'

./bedGraphToBigWig $PREFIX'_unimappers.coverage.minus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size $PREFIX'_unimappers.coverage.minus.rpm_test.bw'


