#!/bin/bash

# written by Marcel and modified by Seung
# This file is for Raw read genome browser 

###########
# version
###########
# sambamba 0.6.8
# bedtools v2.27.1

#FILE=$1
PREFIX=$1

# Run e.g.
#./BigWig_Create_Raw_Count.sh SP064_022


# unimapper BAM file
bedtools genomecov -ibam - -bg -split -strand '+' < BAM_files/$PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > BAM_files/$PREFIX'_unimappers.not_normalized.coverage.plus.rpm.bedgraph'


bedtools genomecov -ibam - -bg -split -strand '-' < BAM_files/$PREFIX'_unimappers.sorted.bam' | sort --parallel=12 --key=1,1 --key=2,3n > BAM_files/$PREFIX'_unimappers.not_normalized.coverage.minus.rpm.bedgraph'


#################
# 4. BigWig
#################

./bedGraphToBigWig BAM_files/$PREFIX'_unimappers.not_normalized.coverage.plus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size BAM_files/$PREFIX'_unimappers.not_normalized.coverage.plus.rpm_test.bw'

./bedGraphToBigWig BAM_files/$PREFIX'_unimappers.not_normalized.coverage.minus.rpm.bedgraph' /data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/GRCm38-3.1.0_mrna/fasta/genome.size BAM_files/$PREFIX'_unimappers.not_normalized.coverage.minus.rpm_test.bw'
