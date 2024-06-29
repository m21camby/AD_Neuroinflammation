#!/bin/bash

File1=$1
File2=$2
ARG=$3
# from NR_140 ~ NR_143 
# e.g. ./DropSeq_pipe.sh ./Schneeberger/NR_140/ds076_S1_R1_001.fastq.gz ./Schneeberger/NR_140/ds076_S1_R2_001.fastq.gz NR_140

# 0. fastq to sam
java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar FastqToSam F1=$File1 F2=$File2 SM=$ARG  O=$ARG'.bam' 2> fastqToSam.STDERR$ARG'.txt' 

# 1. Cell barcode tag
/data/murphy/shared_bins/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended \ 
SUMMARY=Reports/unaligned_tagged_Cellular.bam_summary_$ARG'.txt' \ 
BASE_RANGE=1-12 \ 
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=false \ 
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 \ 
INPUT=$ARG'.bam' \
OUTPUT=unaligned_tagged_Cell_$ARG'.bam' \ 
COMPRESSION_LEVEL=0 2> CellBarcodeTag.STDERR$ARG'.txt' 

# 2. Molecular barcode tag
/data/murphy/shared_bins/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended \
SUMMARY=Reports/unaligned_tagged_Molecular.bam_summary$ARG'.txt' \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=true \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=unaligned_tagged_Cell_$ARG'.bam' \
OUTPUT=unaligned_tagged_CellMole_$ARG'.bam'\
COMPRESSION_LEVEL=0 2> MoleBarcodeTag.STDERR$ARG'txt'

# 3. XQ tag removed (quality score below will remove)
/data/murphy/shared_bins/Drop-seq_tools-1.12/FilterBAM \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMole_$ARG'.bam' \
OUTPUT=unaligned_tagged_filtered_$ARG'.bam' \
COMPRESSION_LEVEL=0 2> XQTag.STDERR$ARG'.txt'

# 4. SMART tag removed (5')
/data/murphy/shared_bins/Drop-seq_tools-1.12/TrimStartingSequence \
OUTPUT_SUMMARY=Reports/adapter_trimming_report_$ARG'.txt' \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5 \
INPUT=unaligned_tagged_filtered_$ARG'.bam' \
OUTPUT=unaligned_tagged_trimmed_smart_$ARG'.bam' \
COMPRESSION_LEVEL=0 2> SMART_Tag.STDERR$ARG'.txt'

# 5. PolyA remove (3')
/data/murphy/shared_bins/Drop-seq_tools-1.12/PolyATrimmer \
OUTPUT_SUMMARY=Reports/polyA_trimming_report_$ARG'.txt' \
MISMATCHES=0 \
NUM_BASES=6 \
INPUT=unaligned_tagged_trimmed_smart_$ARG'.bam' \
OUTPUT=unaligned_mc_tagged_polyA_filtered_$ARG'.bam' \
COMPRESSION_LEVEL=0 2> PolyA.STDERR$ARG'.txt'

# 6. Back to fastq
java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered_$ARG'.bam' \
FASTQ=$ARG'.fastq' 2> SamToFastq.STDERR$ARG'.txt'

# 7. STAR mapping
# I create STAR Index in advance
STAR \
--runThreadN 24 \
--genomeDir /data/rajewsky/indices/mm10_GRCm38.p5_STAR_2.6.0a/ \
--readFilesIn /data/murphy/home/skim/Microglia_Heppner/$ARG'.fastq' \
--outFileNamePrefix $ARG 2> STARMapping.STDERR$ARG'.txt'

# 8. SortSam
# This will continue on DropSeq_Pipe 2nd bash file






