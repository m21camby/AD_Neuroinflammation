#!/bin/bash

#####################################
# output filtered file: other biotypes such as gene_biotype:pseudogene are excluded from the GTF annotation.
# ref: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_3.1.0
#####################################

/data/rajewsky/shared_bins/cellranger-3.1.0/cellranger mkgtf /data/rajewsky/annotation/GRCm38/Mus_musculus.GRCm38.98.gtf Mus_musculus.GRCm38.98.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene
