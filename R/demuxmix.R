# Internal method dmmPreprocess
#' @importFrom stats kmeans median quantile
dmmPreprocess <- function(hto, rna, k.hto=1.5, k.rna=1.5) {
  
  # Initial k-means
  # Cluster ID=1 identifies cluster with lower hto counts
  km <- kmeans(log(hto + 10), centers=2)
  ind <- km$cluster == 1
  if (median(hto[ind]) <= median(hto[!ind])) {
    clusterInit <- km$cluster
  } else {
    clusterInit <- c(2, 1)[km$cluster]
  }
  
  # HTO outliers
  qrtls <- quantile(hto[clusterInit == 2], probs=c(0.25, 0.75))
  th <- qrtls[2] + (k.hto * (qrtls[2] - qrtls[1]))
  outlierHto <- hto > th
  
  # RNA outliers
  if (missing(rna)) {
    if (sum(outlierHto)/length(outlierHto) >= 0.03) {
      warning("More than 3% of the cells have been identified as outliers and excluded from model fitting. Consider increasing k.hto and check the initial clustering by demuxmix:::dmmPreprocess.")
    }
    return(data.frame(hto=hto, cluster=clusterInit, outlier=outlierHto))
  } else {
    qrtls <- quantile(rna, probs=c(0.25, 0.75))
    th <- qrtls[2] + (k.hto * (qrtls[2] - qrtls[1]))
    outlierRna <- rna > th
    
    if (sum(outlierHto | outlierRna)/length(outlierHto) >= 0.05) {
      warning("More than 5% of the cells have been identified as outliers and excluded from model fitting. Consider increasing k.hto and/or k.rna and check the initial clustering by demuxmix:::dmmPreprocess.")
    }
    return(data.frame(hto=hto, rna=rna, cluster=clusterInit, outlier=outlierHto | outlierRna))
  }
}



# Internal method dmmFitNaive
#' @importFrom stats var dnbinom predict glm.control coefficients
dmmFitNaive <- function(hto, clusterInit, tol=10^-5, maxIter=100, htoId="HTO") {
  
  stopifnot(all(clusterInit %in% c(1,2)))
  stopifnot(all(!is.na(hto)))
  
  # Initialize variables
  ind <- clusterInit == 1
  pi <- c(sum(ind), sum(!ind)) / length(ind)
  mu1 <- mean(hto[ind])
  mu2 <- mean(hto[!ind])
  var1 <- var(hto[ind])
  var2 <- var(hto[!ind])
  theta1 <- mu1^2 / (var1 - mu1)
  theta2 <- mu2^2 / (var2 - mu2)
  if (mu1 > mu2) {
    stop("Mean count of cluster 1 must be smaller than mean of cluster 2. Please swap cluster IDs in the clusterInit parameter.")
  }
  
  # Run EM algorithm
  iter <- 1
  logLik <- 1
  deltaLogLik <- Inf
  logData <- data.frame()
  while (iter < maxIter & abs(deltaLogLik/logLik) > tol) {
    
    # E-step
    f <- matrix(nrow=length(hto), ncol=2)
    f[, 1] <- dnbinom(hto, mu=mu1, size=theta1)
    f[, 2] <- dnbinom(hto, mu=mu2, size=theta2)
    pi.f <- t(pi * t(f))
    posteriorProb <- pi.f / apply(pi.f, 1, sum)
    pi <- apply(posteriorProb, 2, sum) / nrow(posteriorProb)
    
    # M-step
    w <- apply(posteriorProb, 2, function(x) {x / sum(x)})
    mu1 <- sum(w[, 1] * hto)
    var1 <- sum(w[, 1] * (hto - mu1)^2)
    theta1 <- mu1^2/ (var1 - mu1)
    mu2 <- sum(w[, 2] * hto)
    var2 <- sum(w[, 2] * (hto - mu2)^2)
    theta2 <- mu2^2/ (var2 - mu2)
    
    # Update counter and keep some logs
    iter <- iter + 1
    oldLogLik <- logLik
    logLik <- sum(log(apply(pi.f, 1, sum)))
    deltaLogLik <- logLik - oldLogLik
    logData <- rbind(logData,
                     data.frame(Iteration=iter-1,
                                LogLikelihood=logLik,
                                DeltaLogLikelihood=deltaLogLik,
                                mu1=mu1,
                                mu2=mu2,
                                theta1=theta1,
                                theta2=theta2))
  }
  
  if (!iter < maxIter) {
    warning(paste("Maximum number of iterations reached before convergence (", htoId, ").", sep=""))
  }
  rownames(logData) <- NULL
  
  mixModel <- NaiveMixModel(mu1=mu1, mu2=mu2, theta1=theta1, theta2=theta2, pi=pi,
                            hto=hto, htoId=htoId,
                            parameters=list(tol=tol, maxIter=maxIter),
                            log=logData, converged=iter<maxIter)
  return(mixModel)
}



# Internal method dmmFitReg
#' @importFrom stats var dnbinom predict glm.control coefficients
#' @importFrom MASS glm.nb
dmmFitReg <- function(hto, rna, clusterInit, regRnaNegComp=TRUE, tol=10^-5, maxIter=100, htoId="HTO") {
  
  stopifnot(all(clusterInit %in% c(1,2)))
  stopifnot(all(!is.na(hto)))
  stopifnot(all(!is.na(rna)))
  stopifnot(all(rna > 0))
  rna <- log(rna)
  
  # Initialize variables
  ind <- clusterInit == 1
  pi <- c(sum(ind), sum(!ind)) / length(ind)
  mu1 <- mean(hto[ind])
  mu2 <- mean(hto[!ind])
  var1 <- var(hto[ind])
  var2 <- var(hto[!ind])
  theta1 <- mu1^2 / (var1 - mu1)
  theta2 <- mu2^2 / (var2 - mu2)
  if (mu1 > mu2) {
    stop("Mean count of cluster 1 must be smaller than mean of cluster 2. Please swap cluster IDs in the clusterInit parameter.")
  }
  
  # Run EM algorithm
  iter <- 1
  logLik <- 1
  deltaLogLik <- Inf
  logData <- data.frame()
  while (iter < maxIter & abs(deltaLogLik/logLik) > tol) {
    
    # E-step
    f <- matrix(nrow=length(hto), ncol=2)
    f[, 1] <- dnbinom(hto, mu=mu1, size=theta1)
    f[, 2] <- dnbinom(hto, mu=mu2, size=theta2)
    pi.f <- t(pi * t(f))
    posteriorProb <- pi.f / apply(pi.f, 1, sum)
    pi <- apply(posteriorProb, 2, sum) / nrow(posteriorProb)
    
    # M-step
    if (regRnaNegComp) {
      fit1 <- glm.nb(hto ~ rna, weights=posteriorProb[, 1], control=glm.control(maxit=maxIter))
    } else {
      fit1 <- glm.nb(hto ~ 1, weights=posteriorProb[, 1], control=glm.control(maxit=maxIter))
    }
    mu1 <- predict(fit1, type="response")
    theta1 <- fit1$theta
    fit2 <- glm.nb(hto ~ rna, weights=posteriorProb[, 2], control=glm.control(maxit=maxIter))
    mu2 <- predict(fit2, type="response")
    theta2 <- fit2$theta
    
    # Update counter and keep some logs
    iter <- iter + 1
    oldLogLik <- logLik
    logLik <- sum(log(apply(pi.f, 1, sum)))
    deltaLogLik <- logLik - oldLogLik
    logData <- rbind(logData,
                     data.frame(Iteration=iter-1,
                                LogLikelihood=logLik,
                                DeltaLogLikelihood=deltaLogLik,
                                mu1=predict(fit1, newdata=data.frame(rna=sum(rna * posteriorProb[, 1]) / sum(posteriorProb[, 1])), type="response"),
                                mu2=predict(fit2, newdata=data.frame(rna=sum(rna * posteriorProb[, 2]) / sum(posteriorProb[, 2])), type="response"),
                                theta1=theta1,
                                theta2=theta2,
                                intercept1=coefficients(fit1)[1],
                                intercept2=coefficients(fit2)[1],
                                coef1=coefficients(fit1)[2],
                                coef2=coefficients(fit2)[2]))
  }
  
  if (!iter < maxIter) {
    warning(paste("Maximum number of iterations reached before convergence (", htoId, ").", sep=""))
  }
  rownames(logData) <- NULL

  mixModel <- RegMixModel(fit1=fit1, fit2=fit2, pi=pi, htoId=htoId,
                          parameters=list(regRnaNegComp=regRnaNegComp, tol=tol, maxIter=maxIter),
                          log=logData, converged=iter<maxIter)
  return(mixModel)
}



# Internal method dmmApplyModel
#' @importFrom stats dnbinom pnbinom predict
dmmApplyModel <- function(model, hto, rna, alpha=0.9, beta=0.9, correctTails=TRUE) {
  
  if (is(model, "NaiveMixModel")) {
    stopifnot(all(!is.na(hto)))
    mu1 <- model@mu1
    mu2 <- model@mu2
    theta1 <- model@theta1
    theta2 <- model@theta2
    pi <- model@pi
    
  } else if (is(model, "RegMixModel")) {
    stopifnot(all(!is.na(hto)))
    stopifnot(all(!is.na(rna)))
    stopifnot(all(rna > 0))
    newData <- data.frame(hto=hto, rna=log(rna))
    mu1 <- predict(model@fit1, newdata=newData, type="response")
    mu2 <- predict(model@fit2, newdata=newData, type="response")
    theta1 <- model@fit1$theta
    theta2 <- model@fit2$theta
    pi <- model@pi
  }
  
  f <- matrix(nrow=length(hto), ncol=2)
  f[, 1] <- dnbinom(hto, mu=mu1, size=theta1)
  f[, 2] <- dnbinom(hto, mu=mu2, size=theta2)
  pi.f <- t(pi * t(f))
  posteriorP <- (pi.f / apply(pi.f, 1, sum))[, 2]
  Pneg.grX <- pnbinom(hto - 1, mu=mu1, size=theta1, lower.tail=FALSE)  # P(X >= x)
  Ppos.smX <- pnbinom(hto, mu=mu2, size=theta2, lower.tail=TRUE)       # P(X <= x)
  posteriorP[is.na(posteriorP) & Ppos.smX == 1] <- 1  # Extreme large HTO counts where both components have P(X=x) = 0
  stopifnot(!any(is.na(posteriorP)))
  
  tailException <- rep(FALSE, length(hto))
  if (correctTails) {  # adjust posterior probability
    ind <- posteriorP >= 0.5 & Pneg.grX > alpha # positive component with heavy left tail
    tailException[ind] <- TRUE
    posteriorP[ind] <- 0
    ind <- posteriorP < 0.5 & Ppos.smX > beta  # negative component with heavy right tail (never observed in real data)
    tailException[ind] <- TRUE
    posteriorP[ind] <- 1
  } else {             # just note the problematic classifications but don't correct anything
    tailException[posteriorP >= 0.5 & Pneg.grX > alpha] <- TRUE
    tailException[posteriorP < 0.5 & Ppos.smX > beta] <- TRUE
  }

  return(data.frame(posteriorP=posteriorP, Pneg.grX=Pneg.grX, Ppos.smX=Ppos.smX, tailException=tailException))
}


# Internal method .demuxmix
.demuxmix <- function(hto, rna, p.acpt=0.9^nrow(hto), model="auto", alpha=0.9, beta=0.9, correctTails=TRUE, tol=10^-5, maxIter=100, k.hto=1.5, k.rna=1.5) {
  
  # Check parameters
  n <- nrow(hto)
  if (is.null(rownames(hto)) || any(duplicated(rownames(hto)))) {
    stop("Matrix hto must have unique row names identifying the HTOs.")
  }
  if (any(hto < 0) | any(is.na(hto))) {
    stop("Matrix hto must not contain negative or missing values.")
  }
  if (!all(is.element(model, c("auto", "naive", "reg", "regpos")))) {
    stop("Parameter model must be either \"auto\", \"naive\", \"reg\" or \"regpos\".")
  }
  if (any(is.element(model, c("auto", "reg", "regpos"))) & missing(rna)) {
    stop("Try to set model=\"naive\" to run demumix without RNA data. Parameter \"rna\" must be given if model is \"auto\" (default), \"reg\" or \"regpos\".")
  }
  if (any(is.element(model, c("auto", "reg", "regpos"))) && ncol(hto) != length(rna)) {
    stop("Length of rna must equal the number of columns in hto.")
  }
  if (any(is.element(model, c("auto", "reg", "regpos"))) && any(rna < 1)) {
    stop("Values in rna must be larger or equal to 1.")
  }
  if (length(p.acpt) != 1) {
    stop("Only one value for p.acpt must be given.")
  }
  if (p.acpt < 0 | p.acpt > 1) {
    stop("p.acpt must be between 0 and 1.")
  }
  model <- rep_len(model, length.out=n)
  alpha <- rep_len(alpha, length.out=n)
  beta <- rep_len(beta, length.out=n)
  correctTails <- rep_len(correctTails, length.out=n)
  tol <- rep_len(tol, length.out=n)
  maxIter <- rep_len(maxIter, length.out=n)
  k.hto <- rep_len(k.hto, length.out=n)
  k.rna <- rep_len(k.rna, length.out=n)
  
  # Run preprocessing, model fitting, and classification for each hashtag
  outliers <- matrix(NA, nrow=nrow(hto), ncol=ncol(hto))
  rownames(outliers) <- rownames(hto)
  colnames(outliers) <- colnames(hto)
  clusterInit <- matrix(NA, nrow=nrow(hto), ncol=ncol(hto))
  rownames(clusterInit) <- rownames(hto)
  colnames(clusterInit) <- colnames(hto)
  posteriorProb <- matrix(NA, nrow=nrow(hto), ncol=ncol(hto))
  rownames(posteriorProb) <- rownames(hto)
  colnames(posteriorProb) <- colnames(hto)
  tailException <- matrix(NA, nrow=nrow(hto), ncol=ncol(hto))
  rownames(tailException) <- rownames(hto)
  colnames(tailException) <- colnames(hto)
  mixModels <- list()
  modelSelection <- data.frame()
  for (i in 1:n) {
    if (model[i] == "naive") {
      prepData <- dmmPreprocess(hto[i, ], k.hto=k.hto[i])
      mixModels[[i]] <- dmmFitNaive(hto=prepData$hto[!prepData$outlier],
                                    clusterInit=prepData$cluster[!prepData$outlier],
                                    tol=tol[i], maxIter=maxIter[i],
                                    htoId=rownames(hto)[i])
      outliers[i, ] <- prepData$outlier
      clusterInit[i, ] <- prepData$cluster
      res <- dmmApplyModel(mixModels[[i]], hto=hto[i, ], alpha=alpha[i], beta=beta[i], correctTails=correctTails[i])
      
    } else if (model[i] == "reg" | model[i] == "regpos") {
      prepData <- dmmPreprocess(hto[i, ], rna, k.hto=k.hto[i], k.rna=k.rna[i])
      mixModels[[i]] <- dmmFitReg(hto=prepData$hto[!prepData$outlier], rna=prepData$rna[!prepData$outlier],
                                  clusterInit=prepData$cluster[!prepData$outlier],
                                  regRnaNegComp=model[i] == "reg",
                                  tol=tol[i], maxIter=maxIter[i],
                                  htoId=rownames(hto)[i])
      outliers[i, ] <- prepData$outlier
      clusterInit[i, ] <- prepData$cluster
      res <- dmmApplyModel(mixModels[[i]], hto=hto[i, ], rna=rna, alpha=alpha[i], beta=beta[i], correctTails=correctTails[i])
      
    } else if (model[i] == "auto") {
      prepDataN <- dmmPreprocess(hto[i, ], k.hto=k.hto[i])
      prepDataR <- dmmPreprocess(hto[i, ], rna, k.hto=k.hto[i], k.rna=k.rna[i])
      mmNaive <- dmmFitNaive(hto=prepDataN$hto[!prepDataN$outlier],
                             clusterInit=prepDataN$cluster[!prepDataN$outlier],
                             tol=tol[i], maxIter=maxIter[i],
                             htoId=rownames(hto)[i])
      mmRegpos <- dmmFitReg(hto=prepDataR$hto[!prepDataR$outlier], rna=prepDataR$rna[!prepDataR$outlier],
                            clusterInit=prepDataR$cluster[!prepDataR$outlier],
                            regRnaNegComp=FALSE,
                            tol=tol[i], maxIter=maxIter[i],
                            htoId=rownames(hto)[i])
      mmReg <- dmmFitReg(hto=prepDataR$hto[!prepDataR$outlier], rna=prepDataR$rna[!prepDataR$outlier],
                         clusterInit=prepDataR$cluster[!prepDataR$outlier],
                         regRnaNegComp=TRUE,
                         tol=tol[i], maxIter=maxIter[i],
                         htoId=rownames(hto)[i])
      ov.mmNaive <- .dmmOverlap(mmNaive)
      ov.mmRegpos <- .dmmOverlap(mmRegpos)
      ov.mmReg <- .dmmOverlap(mmReg)
      
      if (ov.mmNaive <= ov.mmRegpos & ov.mmNaive <= ov.mmReg) {
        bestModel <- "naive"
        mixModels[[i]] <- mmNaive
        outliers[i, ] <- prepDataN$outlier
        clusterInit[i, ] <- prepDataN$cluster
        res <- dmmApplyModel(mmNaive, hto=hto[i, ], alpha=alpha[i], beta=beta[i], correctTails=correctTails[i])
      } else if (ov.mmRegpos <= ov.mmReg) {
        bestModel <- "regpos"
        mixModels[[i]] <- mmRegpos
        outliers[i, ] <- prepDataR$outlier
        clusterInit[i, ] <- prepDataR$cluster
        res <- dmmApplyModel(mmRegpos, hto=hto[i, ], rna=rna, alpha=alpha[i], beta=beta[i], correctTails=correctTails[i])
      } else {
        bestModel <- "reg"
        mixModels[[i]] <- mmReg
        outliers[i, ] <- prepDataR$outlier
        clusterInit[i, ] <- prepDataR$cluster
        res <- dmmApplyModel(mmReg, hto=hto[i, ], rna=rna, alpha=alpha[i], beta=beta[i], correctTails=correctTails[i])
      }
      modelSelection <- rbind(modelSelection,
                              data.frame(hto=rownames(hto)[i], ov.naive=ov.mmNaive, ov.regpos=ov.mmRegpos, ov.reg=ov.mmReg, best=bestModel))
    
    } else {
      stop("Unknown model specified.") # Should never happen.
    }
  
    tailException[i, ] <- res$tailException
    posteriorProb[i, ] <- res$posteriorP
  }
      
  names(mixModels) <- rownames(hto)
  rownames(modelSelection) <- NULL

  parameters <- list(p.acpt=p.acpt, model=model, alpha=alpha, beta=beta, correctTails=correctTails, tol=tol, maxIter=maxIter, k.hto=k.hto, k.rna=k.rna)
  
  dmm <- Demuxmix(models=mixModels,
                  outliers=outliers,
                  clusterInit=clusterInit,
                  posteriorProb=posteriorProb,
                  tailException=tailException,
                  modelSelection=modelSelection,
                  parameters=parameters)
  return(dmm)
}



#' @importFrom methods setMethod
#' @importClassesFrom Matrix Matrix
setMethod("demuxmix", signature=c(hto="Matrix", rna="missing"), .demuxmix)
setMethod("demuxmix", signature=c(hto="matrix", rna="missing"), .demuxmix)
setMethod("demuxmix", signature=c(hto="Matrix", rna="numeric"), .demuxmix)
setMethod("demuxmix", signature=c(hto="matrix", rna="numeric"), .demuxmix)