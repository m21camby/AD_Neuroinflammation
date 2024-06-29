#!/bin/bash

#############################################
# Written by SJK For Maryna
# example of bash file for Drop-seq pipeline
############################################

###############################
# 3 arguments used in the file
# File1 is first fastq file
# File2 is second fastq file
# ARG1 is prefix of the file
###############################

File1=$1
File2=$2
ARG1=$3

############################################################################################
# Example of how to use the file
# example data:  NR_SS_001 
# Should change permission of this file before use (e.g. chmod +x DropSeq_Pipe.sh)
# Command: ./DropSeq_Pipe.sh ./Schneeberger/NR_SS_001/ss01_1_S1_R1_001.fastq.gz ./Schneeberger/NR_SS_001/ss01_1_S1_R2_001.fastq.gz ss_001
# Should aware that this file is for mouse. So you should change some parameters of each step
# e.g. STAR mapping --genomeDir is my own index so you should change to your own
############################################################################################ 

###################
# 0. fastq to sam
###################
java -jar /data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar FastqToSam \
F1=$File1 \
F2=$File2 \
SM=$ARG1 \
O=$ARG1'.bam' 2> S0_fastqToSam.STDERR$ARG1'.txt'

######################
# 1. Cell barcode tag
######################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
INPUT=$ARG1'.bam' \
OUTPUT=unaligned_tagged_Cell_$ARG1'.bam' \
SUMMARY=unaligned_tagged_Cellular.bam_summary_$ARG1'.txt' \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=false \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 \
COMPRESSION_LEVEL=0 2> S1_CellBarcodeTag.STDERR$ARG1'.txt'

###########################
# 2. Molecular barcode tag
###########################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
SUMMARY=unaligned_tagged_Molecular.bam_summary$ARG1'.txt' \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=true \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=unaligned_tagged_Cell_$ARG1'.bam' \
OUTPUT=unaligned_tagged_CellMole_$ARG1'.bam' \
COMPRESSION_LEVEL=0 2> S2_MoleBarcodeTag.STDERR$ARG1'txt'

#######################################################
# 3. XQ tag removed (quality score below will remove)
#######################################################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/FilterBam \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMole_$ARG1'.bam' \
OUTPUT=unaligned_tagged_filtered_$ARG1'.bam' \
COMPRESSION_LEVEL=0 2> S3_XQTag.STDERR$ARG1'.txt'

#############################
# 4. SMART tag removed (5')
#############################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/TrimStartingSequence \
OUTPUT_SUMMARY=adapter_trimming_report_$ARG1'.txt' \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5 \
INPUT=unaligned_tagged_filtered_$ARG1'.bam' \
OUTPUT=unaligned_tagged_trimmed_smart_$ARG1'.bam' \
COMPRESSION_LEVEL=0 2> S4_SMART_Tag.STDERR$ARG1'.txt'

#############################
# 5. PolyA remove (3')
#############################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/PolyATrimmer \
OUTPUT_SUMMARY=polyA_trimming_report_$ARG1'.txt' \
MISMATCHES=0 \
NUM_BASES=6 \
INPUT=unaligned_tagged_trimmed_smart_$ARG1'.bam' \
OUTPUT=unaligned_mc_tagged_polyA_filtered_$ARG1'.bam' \
COMPRESSION_LEVEL=0 2> S5_PolyA.STDERR$ARG1'.txt'

###################
# 6. Back to fastq
###################
java -jar /data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered_$ARG1'.bam' \
FASTQ=$ARG1'.fastq' 2> S6_SamToFastq.STDERR$ARG1'.txt'

###################################
# 7. STAR mapping
# I create STAR Index in advance
###################################
STAR --runThreadN 24 --genomeDir /data/rajewsky/indices/mm10_GRCm38.p6_star_2.7.0a/ \
--readFilesIn ./$ARG1'.fastq' \
--outFileNamePrefix $ARG1 2> S7_STARMapping.STDERR$ARG1'.txt'

#################
# 8. SortSam 
#################
java -jar /data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar SortSam \
I=$ARG1'Aligned.out.sam' \
O=$ARG1'Aligned.out.sorted.bam' \
SO=queryname 2> S8_SortSam_STDERR$ARG1'.txt'

####################################
# 9. MergaBAMAlgiment 
####################################

java -jar /data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar MergeBamAlignment \
REFERENCE_SEQUENCE=../mm10/mm10.fa \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_$ARG1'.bam' \
ALIGNED_BAM=$ARG1'Aligned.out.sorted.bam' \
OUTPUT=merged_$ARG1'.bam' \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false 2> S9_MergeBam.STDERR$ARG1'.txt'

#################################
# 10. TagReadWithGeneExon
#################################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
I=merged_$ARG1'.bam' \
O=star_gene_exon_tagged_$ARG1'.bam' \
ANNOTATIONS_FILE=GRCm38.M21.refFlat 2> S10_TagReadWith.STDERR$ARG1'.txt';

###################################
# 11. DetectBeadSubstitutionErrors
###################################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors \
I=star_gene_exon_tagged_$ARG1'.bam' \
O=my_clean_subtitution_$ARG1'.bam' \
OUTPUT_REPORT=my_clean_substitution_report_$ARG1'.txt' 2> S11_DetectBead_Sub_$ARG1'.txt'


################################
# 12. DetectBeadSynthesisErrors
################################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors \
I=my_clean_subtitution_$ARG1'.bam' \
O=my_clean_$ARG1'.bam' \
REPORT=my_clean.indel_report_$ARG1'.txt' \
OUTPUT_STATS=my.synthesis_stats_$ARG1'.txt' \
SUMMARY=my.synthesis_stats.summary_$ARG1'.txt' \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC 2> S12_DetectBead_Errors.STDERR_$ARG1'.txt'

###################################
# 13. BAMTAgHistogram
###################################
/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0/BamTagHistogram \
I=my_clean_$ARG1'.bam' \
O=out_cell_readcounts_$ARG1'.txt.gz' \
TAG=XC 2> S13_BAMTagHisto.STDERR_$ARG1'.txt'


################################################
# 14. Digital Gene expression
# Do separately as it requires NUM_CORE_BARCODES
################################################



