
############################
# create for exon reference
############################

#/data/rajewsky/shared_bins/cellranger-3.1.0/cellranger mkref --genome=GRCm38-3.1.0_mrna \
#                   --fasta=/data/rajewsky/genomes/GRCm38_98/Mus_musculus.GRCm38.dna.primary_assembly.fa \
#                   --genes=/data/rajewsky/annotation/GRCm38/Mus_musculus.GRCm38.98.gtf

############################
# create premrna gtf file
############################

#awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' \
#              GRCm38-3.1.0_mrna/genes/genes.gtf > Mus_musculus.GRCm38.98.premrna.gtf

###########################
# create for intron reference
###########################

/data/rajewsky/shared_bins/cellranger-3.1.0/cellranger mkref --genome=GRCm38-3.1.0_premrna \
                   --fasta=/data/rajewsky/genomes/GRCm38_98/Mus_musculus.GRCm38.dna.primary_assembly.fa \
                   --genes=Mus_musculus.GRCm38.98.premrna.gtf

