#' Demultiplexing using mixture models
#'
#' This methods uses mixture models as a probabilistic framework to assign
#' hashtags to cells and to identify multiplets based on counts obtained from
#' sequencing a hastag oligo (HTO) library. If the numbers of detected genes
#' from the corresponding RNA library are passed in the second argument,
#' regression mixture models are used, which usually improves the
#' classification accuracy by leveraging the relationship between HTO and RNA
#' counts.
#' 
#' @param object A matrix of HTO counts where each row corresponds to a hashtag
#'   and each column to a cell barcode.
#' @param rna An optional numeric vector with the number of genes detected in
#'   the RNA library for each cell. Same length as columns in \code{object}.
#'   If missing, a naive instead of a regression mixture model is used.
#' @param p.acpt Acceptance probability that must be reached in order to
#'   assign a cell to a hashtag. Cells with lower probabilities are classified
#'   as "uncertain".
#' @param alpha Threshold for defining the left tail of the mixture
#'   distribution where cells should not be classified as "positive". Either a
#'   single value between 0 and 1 or a vector with one entry for each hashtag
#'   in the dataset. See details.
#' @param beta Threshold for defining the right tail of the mixture
#'   distribution where cells should not be classified as "negative". Either a
#'   single value between 0 and 1 or a vector with one entry for each hashtag
#'   in the dataset. See details.
#' @param classifyTails If \code{TRUE}, cells meeting the threshold defined by
#'   \code{alpha} (\code{beta}) are classified as "negative" ("positive") even
#'   if the mixture model suggests a different classification. See details.
#' @param tol Convergence criterion for the EM algorithm fitting the mixture
#'   model(s). The algorithm stops when the relative increase of the log
#'   likelihood is less than or equal to \code{tol}. Either a single value
#'   or a vector with one entry for each hashtag in the dataset.
#' @param maxIter Maximum number of iterations for the EM algorithm and
#'   for the alternating iteration process fitting the nb regression models
#'   within each EM iteration. Either a single value or a vector with one
#'   entry for each hashtag in the dataset.
#' @param k.hto Factor to define outliers in the HTO counts. Among cells
#'   positive for the hashtag based on initial clustering, HTO counts
#'   larger than the 0.75 quantile + \code{k.hto} * IQR are considered
#'   outliers. Either a single value or a vector with one entry for each
#'   hashtag in the dataset. See details.
#' @param k.rna Factor to define outliers in the numbers of detected genes.
#'   Numbers of detected genes larger than the 0.75 quantile +
#'   \code{k.rna} * IQR are considered outliers. Either a single value or
#'   a vector with one entry for each hashtag in the dataset. See details.
#' 
#' @details The single cell dataset should undergo basic filtering to
#'   remove low quality or empty droplets before calling this function,
#'   but the HTO counts should not be transformed or pre-processed otherwise.
#'   The number of detected genes passed via the \code{rna} argument is
#'   typically defined as the number of genes in the RNA library with at least
#'   one read.
#'   
#'   The method itself consists of three steps:
#'   
#' @return demux
#' @examples atze <- 3
#' @aliases demuxmix,matrix,missing-method demuxmix,matrix,numeric-method
#'   demuxmix,Matrix,missing-method demuxmix,Matrix,numeric-method
#' @export
setGeneric("demuxmix",
           function(object, rna, p.acpt=0.9, alpha=0.9, beta=0.9, classifyTails=TRUE, tol=10^-5, maxIter=100, k.hto=1.5, k.rna=1.5)
             standardGeneric("demuxmix"),
           signature=c("object", "rna"))

#' @export
setGeneric("dmmOverlap",
           function(model, tol=0.001)
            standardGeneric("dmmOverlap"),
           signature="model")

#' @export
setGeneric("dmmSummary",
           function(dmmResults, p.acpt)
             standardGeneric("dmmSummary"),
           signature="dmmResults")

#' @export
setGeneric("plotDmmHistogram",
           function(model,  quantile=0.95, binwidth=50)
             standardGeneric("plotDmmHistogram"),
           signature="model")

#' @export
setGeneric("plotDmmScatter",
           function(model, log=TRUE, pointsize=1.2)
             standardGeneric("plotDmmScatter"),
           signature="model")

#' @export
setGeneric("plotDmmPosteriorP",
           function(model, log=FALSE, bins=50)
             standardGeneric("plotDmmPosteriorP"),
           signature="model")