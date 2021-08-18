#' @importFrom stats rnbinom
.dmmSimulateHto <- function(class, mu=150, theta=20, muAmbient=30, thetaAmbient=10, muRna=2500, thetaRna=40) {
  
  stopifnot(typeof(class) == "logical")
  stopifnot(!any(is.na(class)))
  nHtos <- nrow(class)
  if (is.null(rownames(class))) {
    rownames(class) <- paste("HTO", 1:nHtos, sep="_")
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
  for (i in 1:nrow(counts)) {
    muH <- round(rna * mu[i] / muRna)
    muAmbH <- round(rna * muAmbient[i] / muRna)
    counts[i, class[i, ]] <- rnbinom(n=sum(class[i, ]), mu=muH, size=theta[i])
    counts[i, !class[i, ]] <- rnbinom(n=sum(!class[i, ]), mu=muAmbH, size=thetaAmbient[i])
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
#
# set.seed(2642)
# simdata <- .dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#                                       c(rep(FALSE, 200), rep(TRUE, 220))))
# 
# dmm <- demuxmix(simdata$hto, p.acpt=0.9)
# table(dmm$results$class, simdata$groundTruth)
# 
# dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
# table(dmmreg$results$class, simdata$groundTruth)


# 
# class <- matrix(FALSE, nrow=4, ncol=2050)
# class[1, 1:500] <- TRUE
# class[2, 501:1000] <- TRUE
# class[3, 1001:1500] <- TRUE
# class[4, 1501:2000] <- TRUE
# class[c(1,2), 2001:2050] <- TRUE
# 
# sim <- .dmmSimulateHto(class)
# 
# dmm <- demuxmix(sim$hto, sim$rna)
# table(dmm$results$class, sim$groundTruth)
# 
# dmm <- demuxmix(sim$hto)
# table(dmm$results$class, sim$groundTruth)
# 
# plot(sim$rna, apply(sim$hto, 2, sum))
# cor(sim$rna, apply(sim$hto, 2, sum))
# 
# table(dmm$results$class, sim$groundTruth)
# dmmSummary(dmm$results, p.acpt = 0.9)
#  
# plotDmmHistogram(dmm$model, binwidth = 5)
# 
# class <- factor(rep(c("A", "B", "C", "D"), times=c(500, 500, 500, 500)), levels=c("A", "B", "C", "D"))
# mu <- 300
# theta <- 10
# ambientMu <- 30
# ambientTheta <- 10
# dispersionExplRna <- 0.1
