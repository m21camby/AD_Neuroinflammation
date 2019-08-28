
Ref: https://github.com/mdozmorov/scRNA-seq_notes/blob/master/README.md

## Preprocessing

`SEQC` - Single-Cell Sequencing **Quality Control** and **Processing Software**, a general purpose method to build a count matrix from single cell sequencing reads, able to process data from inDrop, drop-seq, 10X, and Mars-Seq2 technologies.

## Quality Control

`scrublet` - Detect doublets in single-cell RNA-seq data
 
## batch effect

`conos` - joint analysis of scRNA-seq datasets through inter-sample mapping **(mutual nearest-neighbor mapping)** and constructing a joint graph

'MNN` - **mutual nearest neighbors method** for single-cell batch correction. Assumptions: MNN exist between batches, batch is orthogonal to the biology. Cosine normalization, Euclidean distance, a pair-specific barch-correction vector as a vector difference between the expression profiles of the paired cells using selected genes of interest and hypervariable genes.


## Clustering

`Biscuit` - a Bayesian clustering and normalization method


## Trajectory

`Monocle2` - Temporal ordering of single cell gene expression profiles. Reversed graph embedding (DDRTree), finding low-dimensional mapping of differential genes while learning the graph in this reduced space. Allows for the selection of root.

`destiny` - diffusion maps-based visualization of single-cell data

## Differential Expression

`MAST` - scRNA-seq DEG analysis. CDR - the fraction of genes that are detectably expressed in each cell - added to the hurdle model that explicitly parameterizes distributions of expressed and non-expressed genes. Generalized linear model, log2(TPM+1), Gaussian. Regression coeffs are estimated using Bayesian approach.

 




