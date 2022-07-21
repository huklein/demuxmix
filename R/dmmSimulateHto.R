#' @importFrom stats rnbinom
.dmmSimulateHto <- function(class, mu=150, theta=20, muAmbient=30, thetaAmbient=10, muRna=2500, thetaRna=40) {
  
  stopifnot(typeof(class) == "logical")
  stopifnot(!any(is.na(class)))
  nHtos <- nrow(class)
  if (is.null(rownames(class))) {
    rownames(class) <- paste("HTO", seq_len(nHtos), sep="_")
  } else {
    stopifnot(!any(duplicated(rownames(class))))
  }
  
  mu <- rep_len(mu, length.out=nHtos)
  theta <- rep_len(theta, length.out=nHtos)
  muAmbient <- rep_len(muAmbient, length.out=nHtos)
  thetaAmbient <- rep_len(thetaAmbient, length.out=nHtos)
  stopifnot(length(muRna) == 1)
  stopifnot(length(thetaRna) == 1)
  
  rna <- rnbinom(n=ncol(class), mu=muRna, size=thetaRna)
  counts <- matrix(NA, nrow=nrow(class), ncol=ncol(class))
  rownames(counts) <- rownames(class)
  for (i in seq_len(nrow(counts))) {
    muH <- round(rna * mu[i] / muRna)
    muAmbH <- round(rna * muAmbient[i] / muRna)
    counts[i, class[i, ]] <- rnbinom(n=sum(class[i, ]), mu=muH[class[i, ]], size=theta[i])
    counts[i, !class[i, ]] <- rnbinom(n=sum(!class[i, ]), mu=muAmbH[!class[i, ]], size=thetaAmbient[i])
  }
  
  groundTruth <- apply(class, 2, function (ind) {
    if (sum(ind) == 0) {
      return("negative")
    } else {
      return(paste(rownames(counts)[ind], collapse=","))
    }
  })
  
  return(list(hto=counts, rna=rna, groundTruth=groundTruth))
}


setMethod("dmmSimulateHto",
          signature=c(class="matrix"),
          .dmmSimulateHto)