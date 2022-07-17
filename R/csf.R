#' Hashtag oligonucleotide (HTO) counts from 2,590 droplets
#'
#' Cerebral spinal fluid (CSF) cells and peripheral blood mononuclear cells
#' (PBMCs) were pooled and prepared for single-cell sequencing using the
#' 10x Chromium System. Due to the low numbers of cells obtained from CSF, only
#' the PBMCs but not the CSF cels were stained using oligonucleotide-labeled
#' antibodies (BioLegend TotalSeq-A0257). CSF cells and PBMCs in this dataset
#' were obtained from two genetically diverse individuals so that genetic
#' demultiplexing could be used to validate the HTO-based demultiplexing.
#' Genetic demultiplexing was performed with freemuxlet, which is part of the
#' popscle software package.
#'
#' @format A data frame with 2,590 rows and 4 variables:
#' \describe{
#'   \item{HTO}{Number of HTO counts observed}
#'   \item{NumGenes}{Number of genes detected in the cell}
#'   \item{freemuxlet}{Genetic demultiplexing result}
#'   \item{freemuxlet.prob}{Posterior probability from genetic demultiplexing
#'     in logarithmic scale}
#' }
#' 
#' Raw sequencing data was aligned and processed using Cell Ranger 6.0.1. All
#' droplets that passed Cell Ranger's default filtering step were read in. Genes
#' with at least one read were considered as detected. Since Cell Ranger's
#' threshold to identify non-empty droplets is relatively lenient, some droplets
#' have as few as 30-50 genes detected. For most analyses, it is recommended
#' to remove droplets with less than about 200 detected genes
#' before demultiplexing.
#'
#' @source Center for Translational and Computational Neuroimmunology,
#'   Department of Neurology, Columbia University Irving Medical Center,
#'   contact: Hans-Ulrich Klein (hk2948@cumc.columbia.edu)
#' @docType data
#' @name csf
#' @usage data(csf)
#' @examples
#' data(csf)
#' csf <- csf[csf$NumGenes >= 200, ]
#' hto <- t(matrix(csf$HTO, dimnames=list(rownames(csf), "HTO")))
#' dmm <- demuxmix(hto, model="naive")
#' summary(dmm)
#' certain <- exp(csf$freemuxlet.prob) >= 0.999
#' table(dmmClassify(dmm)$HTO[certain], csf$freemuxlet[certain])
"csf"

# Code used to read in the dataset:
# library(SingleCellExperiment)
# library(DropletUtils)
# fm <- read.table("freemuxlet/HK03/HK03.clust1.samples.gz", header=TRUE, stringsAsFactors=FALSE)
# rownames(fm) <- fm$BARCODE
# 
# sce <- read10xCounts("cellRangerCount/HK03/outs/filtered_feature_bc_matrix", col.names=TRUE)
# ncol(sce) # 2590
# fm <- fm[colnames(sce), ]
# 
# csf <- data.frame(HTO=assay(sce)["HTO_7",],
#                   NumGenes=apply(assay(sce)[!grepl("HTO_7", rownames(sce)),] > 0, 2, sum),
#                   freemuxlet=fm$BEST.GUESS,
#                   freemuxlet.prob=fm$BEST.POSTERIOR)
# save(csf, file="csf.RData")