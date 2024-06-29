#!/bin/bash


Samples="SP064_022 SP064_023 SP064_024 SP064_025 SP064_026 SP064_027 SP064_028 SP064_029 SP064_030"

for i in "${Samples[@]}"
do

	echo $i 

	/data/rajewsky/shared_bins/cellranger-3.1.0/cellranger count --id=$i \
	--fastqs=/data/rajewsky/sequencing/mouse/191218_K00302_0148_BHFVK3BBXY/HFVK3BBXY/$i"/" \
	--sample=$i \
	--transcriptome=GRCm38-3.1.0_premrna \
	--localcores=24 \
	--localmem=200 

done
