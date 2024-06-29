#!/bin/bash

ARG1=$1

sambamba view -S --format=bam --compression-level=0 $ARG1'Aligned.out.sam' > $ARG1'Aligned.out.bam' &

sambamba sort -m 8G -t 8 -o $ARG1'Aligned.out.sorted.bam' --tmpdir=tmp $ARG1'Aligned.out.bam' 2> sambambaSort$ARG1'.txt' &

java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=mm10/mm10.fa \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_$ARG1'.bam' \
ALIGNED_BAM=$ARG1'Aligned.out.sorted.bam' \
OUTPUT=merged_$ARG1'.bam' \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false 2> MergeBam.STDERR$ARG1'.txt' &

/data/murphy/shared_bins/Drop-seq_tools-1.12/TagReadWithGeneExon \
I=merged_$ARG1'.bam' \
O=star_gene_exon_tagged_$ARG1'.bam' \
ANNOTATIONS_FILE=gencode.vM15.annotation.refFlat \
TAG=GE 2> TagReadWith.STDERR$ARG'.txt' &

# e.g. NR_140: 3440403 x2 ~ 7000000

samtools view unaligned_tagged_CellMole_$ARG1'.bam' |awk '{print $12}' |sort | uniq|wc -l > $ARG1'cell_number.txt' &

#ARG2= $(awk NF $ARG1'cell_number.txt')
#ARG3= $(($ARG2*2))

#/data/murphy/shared_bins/Drop-seq_tools-1.12/DetectBeadSynthesisErrors \
#I=star_gene_exon_tagged_$ARG1'.bam' \
#O=my_clean_$ARG1'.bam' \
#OUTPUT_STATS=my.synthesis_stats_$ARG1'.txt' \
#SUMMARY=my.synthesis_stats.summary_$ARG1'.txt' \
#NUM_BARCODES=$ARG3 \
#PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC 2> DetectBead.STDERR_$ARG1'.txt' &

#/data/murphy/shared_bins/Drop-seq_tools-1.12/BAMTagHistogram \
#I=my_clean_$ARG1'.bam' \
#O=out_cell_readcounts_$ARG1'.txt.gz' \
#TAG=XC 2> BAMTagHisto.STDERR_$ARG1'.txt' &


####################################################################
# java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar SortSam I=NR_140Aligned.out.sam O=NR_140Aligned.out.sorted.bam SO=queryname
java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar SortSam I=NR_141Aligned.out.sam O=NR_141Aligned.out.sorted.bam SO=queryname  (wd: /data/murphy/home/skim/Microglia_Heppner)
java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar SortSam I=NR_142Aligned.out.sam O=NR_142Aligned.out.sorted.bam SO=queryname  (wd: /data/murphy/home/skim/Microglia_Heppner)
java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar SortSam I=NR_143Aligned.out.sam O=NR_143Aligned.out.sorted.bam SO=queryname  (wd: /data/murphy/home/skim/Microglia_Heppner)

java -jar /data/murphy/shared_bins/picard-tools-2.9.0/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=mm10/mm10.fa UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_NR_141.bam ALIGNED_BAM=NR_141Aligned.out.sorted.bam OUTPUT=merged_NR_141.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false 2> MergeBam.STDERRNR_141.txt &

/data/murphy/shared_bins/Drop-seq_tools-1.12/TagReadWithGeneExon I=merged_NR_141.bam O=star_gene_exon_tagged_NR_141.bam ANNOTATIONS_FILE=gencode.vM15.annotation.refFlat TAG=GE 2> TagReadWith.STDERRNR_141.txt &

samtools view unaligned_tagged_CellMole_NR_141.bam |awk '{print $12}' |sort | uniq|wc -l > NR_141cell_number.txt &

NR_141: 2370877
NR_142: 3296806
NR_143: 3289375

/data/murphy/shared_bins/Drop-seq_tools-1.12/DetectBeadSynthesisErrors I=star_gene_exon_tagged_NR_141.bam O=my_clean_NR_141.bam OUTPUT_STATS=my.synthesis_stats_NR_141.txt SUMMARY=my.synthesis_stats.summary_NR_141.txt NUM_BARCODES=4800000 PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC 2> DetectBead.STDERR_NR_141.txt &

/data/murphy/shared_bins/Drop-seq_tools-1.12/BAMTagHistogram I=my_clean_NR_141.bam O=out_cell_readcounts_NR_141.txt.gz TAG=XC 2> BAMTagHisto.STDERR_NR_141.txt &
        
/data/murphy/shared_bins/Drop-seq_tools-1.12/DigitalExpression I=my_clean_NR_140.bam O=out_gene_exon_tagged.dge_NR_140.txt.gz SUMMARY=out_gene_exon_tagged.dge_NR_140.summary.txt NUM_CORE_BARCODES=99208

NR_141 = 138112
NR_142 = 100493
NR_143 = 136540
        












