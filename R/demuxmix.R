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
    if (sum(outlierHto)/length(outlierHto) >= 0.01) {
      warning("More than 1% of the cells have been identified as outliers and excluded from model fitting. Consider increasing k.hto and check the initial clustering by demuxmix:::dmmPreprocess.")
    }
    return(data.frame(hto=hto, cluster=clusterInit, outlier=outlierHto))
  } else {
    qrtls <- quantile(rna, probs=c(0.25, 0.75))
    th <- qrtls[2] + (k.hto * (qrtls[2] - qrtls[1]))
    outlierRna <- rna > th
    
    if (sum(outlierHto | outlierRna)/length(outlierHto) >= 0.03) {
      warning("More than 3% of the cells have been identified as outliers and excluded from model fitting. Consider increasing k.hto and/or k.rna and check the initial clustering by demuxmix:::dmmPreprocess.")
    }
    return(data.frame(hto=hto, rna=rna, cluster=clusterInit, outlier=outlierHto | outlierRna))
  }
}


# Internal method dmmFit
#' @importFrom stats var dnbinom predict glm.control coefficients
#' @importFrom MASS glm.nb
dmmFit <- function(hto, rna, clusterInit, tol=10^-5, maxIter=100) {
  
  stopifnot(all(clusterInit %in% c(1,2)))
  stopifnot(all(!is.na(hto)))
  if (missing(rna)) {
    regmm <- FALSE
  } else {
    regmm <- TRUE
    stopifnot(all(!is.na(rna)))
    rna <- log(rna)
  }
  
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
    if (regmm) { # regression mixture model
      fit1 <- glm.nb(hto ~ rna, weights=posteriorProb[, 1], control=glm.control(maxit=maxIter))
      mu1 <- predict(fit1, type="response")
      theta1 <- fit1$theta
      fit2 <- glm.nb(hto ~ rna, weights=posteriorProb[, 2], control=glm.control(maxit=maxIter))
      mu2 <- predict(fit2, type="response")
      theta2 <- fit2$theta
    } else {     # naive mixture model
      w <- apply(posteriorProb, 2, function(x) {x / sum(x)})
      mu1 <- sum(w[, 1] * hto)
      var1 <- sum(w[, 1] * (hto - mu1)^2)
      theta1 <- mu1^2/ (var1 - mu1)
      mu2 <- sum(w[, 2] * hto)
      var2 <- sum(w[, 2] * (hto - mu2)^2)
      theta2 <- mu2^2/ (var2 - mu2)
    }
    
    # Update counter and keep some logs
    iter <- iter + 1
    oldLogLik <- logLik
    logLik <- sum(log(apply(pi.f, 1, sum)))
    deltaLogLik <- logLik - oldLogLik
    if (regmm) {
      logData <- rbind(logData,
                       data.frame(Iteration=iter - 1,
                                  LogLikelihood=logLik,
                                  DeltaLogLikelihood=deltaLogLik,
                                  mu1=mean(mu1),
                                  mu2=mean(mu2),
                                  theta1=theta1,
                                  theta2=theta2,
                                  intercept1=coefficients(fit1)[1],
                                  intercept2=coefficients(fit2)[1],
                                  coef1=coefficients(fit1)[2],
                                  coef2=coefficients(fit2)[2]))
    } else {
      logData <- rbind(logData,
                       data.frame(Iteration=iter - 1,
                                  LogLikelihood=logLik,
                                  DeltaLogLikelihood=deltaLogLik,
                                  mu1=mu1,
                                  mu2=mu2,
                                  theta1=theta1,
                                  theta2=theta2))
    }
  }
  
  if (!iter < maxIter) {
    warning("Maximum number of iterations reached before convergence.")
  }
  
  rownames(logData) <- NULL
  if (regmm) {
    mixModel <- list(fit1=fit1, fit2=fit2, pi=pi,
                     logData=logData, converged=iter < maxIter, numCells=length(hto))
  } else {
    mixModel <- list(mu1=mu1, mu2=mu2, theta1=theta1, theta2=theta2, pi=pi, y=hto,
                     logData=logData, converged=iter < maxIter, numCells=length(hto))
  }
  
  return(mixModel)
}


# Internal method dmmClassify
#' @importFrom stats dnbinom pnbinom predict
dmmClassify <- function(hto, rna, mModel, p.acpt=0.9, alpha=0.9, beta=0.9, classifyTails=TRUE) {
  
  stopifnot(all(!is.na(hto)))
  if (missing(rna)) {
    regmm <- FALSE
    stopifnot(all(c("mu1", "mu2", "theta1", "theta2", "pi") %in% names(mModel)))
    mu1 <- mModel$mu1
    mu2 <- mModel$mu2
    theta1 <- mModel$theta1
    theta2 <- mModel$theta2
    pi <- mModel$pi
  } else {
    regmm <- TRUE
    stopifnot(all(!is.na(rna)))
    stopifnot(all(c("fit1", "fit2") %in% names(mModel)))
    predData <- data.frame(hto=hto, rna=log(rna))
    mu1 <- predict(mModel$fit1, newdata=predData, type="response")
    mu2 <- predict(mModel$fit2, newdata=predData, type="response")
    theta1 <- mModel$fit1$theta
    theta2 <- mModel$fit2$theta
    pi <- mModel$pi
  }
  
  f <- matrix(nrow=length(hto), ncol=2)
  f[, 1] <- dnbinom(hto, mu=mu1, size=theta1)
  f[, 2] <- dnbinom(hto, mu=mu2, size=theta2)
  pi.f <- t(pi * t(f))
  posteriorP <- (pi.f / apply(pi.f, 1, sum))[, 2]
  Pneg.grX <- pnbinom(hto, mu=mu1, size=theta1, lower.tail=FALSE)
  Ppos.smX <- pnbinom(hto, mu=mu2, size=theta2, lower.tail=TRUE)
  
  classification <- rep("uncertain", length(hto))
  classification[posteriorP >= p.acpt] <- "positive"
  classification[posteriorP < 1 - p.acpt] <- "negative"
  if (classifyTails) {
    ind <- classification %in% c("positive", "uncertain") & Pneg.grX > alpha # positive component with heavy left tail
    classification[ind] <- "negative"
    posteriorP[ind] <- 0
    ind <- classification %in% c("negative", "uncertain") & Ppos.smX > beta  # negative component with heavy right tail (never observed in real data)
    classification[ind] <- "positive"
    posteriorP[ind] <- 1
  } else {
    classification[classification == "positive" & Pneg.grX > alpha] <- "uncertain"
    classification[classification == "negative" & Ppos.smX > beta] <- "uncertain"
  }

  classRes <- data.frame(class=classification, posteriorP=posteriorP,
                         Pneg.grX=Pneg.grX, Ppos.smX=Ppos.smX)
  return(classRes)
}


# Internal method .demuxmix
.demuxmix <- function(object, rna, p.acpt=0.9, alpha=0.9, beta=0.9, classifyTails=TRUE, tol=10^-5, maxIter=100, k.hto=1.5, k.rna=1.5) {
  
  # Check parameters
  n <- nrow(object)
  stopifnot(!is.null(rownames(object)) & !any(duplicated(rownames(object))))
  if (missing(rna)) {
    regmm <- FALSE
  } else {
    regmm <- TRUE
    stopifnot(length(rna) == ncol(object))
  }
  stopifnot(length(p.acpt) == 1)
  stopifnot(length(alpha) == 1 | length(alpha) == n)
  if (length(alpha) == 1) {
    alpha <- rep(alpha, n)
  }
  stopifnot(length(beta) == 1 | length(beta) == n)
  if (length(beta) == 1) {
    beta <- rep(beta, n)
  }
  stopifnot(length(classifyTails) == 1 | length(classifyTails) == n)
  if (length(classifyTails) == 1) {
    classifyTails <- rep(classifyTails, n)
  }
  stopifnot(length(tol) == 1 | length(tol) == n)
  if (length(tol) == 1) {
    tol <- rep(tol, n)
  }
  stopifnot(length(maxIter) == 1 | length(maxIter) == n)
  if (length(maxIter) == 1) {
    maxIter <- rep(maxIter, n)
  }
  stopifnot(length(k.hto) == 1 | length(k.hto) == n)
  if (length(k.hto) == 1) {
    k.hto <- rep(k.hto, n)
  }
  stopifnot(length(k.rna) == 1 | length(k.rna) == n)
  if (length(k.rna) == 1) {
    k.rna <- rep(k.rna, n)
  }
  
  # Run preprocessing, model fitting, and classification for each hashtag
  mModel <- list()
  classRes <- list()
  for (i in 1:n) {
    if (regmm) {
      prepData <- dmmPreprocess(object[i, ], rna, k.hto=k.hto[i], k.rna=k.rna[i])
      mModel[[i]] <- dmmFit(hto=prepData$hto[!prepData$outlier], rna=prepData$rna[!prepData$outlier],
                            clusterInit=prepData$cluster[!prepData$outlier], tol=tol[i], maxIter=maxIter[i])
      mModel[[i]]$htoId <- rownames(object)[i]
      classRes[[i]] <- dmmClassify(hto=prepData$hto, rna=prepData$rna, mModel=mModel[[i]],
                                   p.acpt=p.acpt, alpha=alpha[i], beta=beta[i], classifyTails=classifyTails[i])
      colnames(classRes[[i]]) <- paste(rownames(object)[i], colnames(classRes[[i]]), sep=".")
    } else {
      prepData <- dmmPreprocess(object[i, ], k.hto=k.hto[i], k.rna=k.rna[i])
      mModel[[i]] <- dmmFit(hto=prepData$hto[!prepData$outlier],
                            clusterInit=prepData$cluster[!prepData$outlier], tol=tol[i], maxIter=maxIter[i])
      mModel[[i]]$htoId <- rownames(object)[i]
      classRes[[i]] <- dmmClassify(hto=prepData$hto, mModel=mModel[[i]],
                                   p.acpt=p.acpt[i], alpha=alpha[i], beta=beta[i], classifyTails=classifyTails[i])
      colnames(classRes[[i]]) <- paste(rownames(object)[i], colnames(classRes[[i]]), sep=".")
    }
  }
  names(mModel) <- rownames(object)
  
  # Combine results from all hashtags
  results <- do.call(cbind, classRes)
  posteriorP <- results[, paste(rownames(object), "posteriorP", sep="."), drop=FALSE]
  mlClassMat <- posteriorP >= 0.5
  posteriorP[!mlClassMat] <- 1 - posteriorP[!mlClassMat]
  mlClass <- apply(mlClassMat, 1, function(c) {return(paste(rownames(object)[c], collapse=","))})
  mlClass[mlClass == ""] <- "negative"
  mlClassP <- apply(posteriorP, 1, prod)
  mlClassType <- c("negative", "singleton", "multiplet")[pmin(apply(mlClassMat, 1, sum), 2) + 1]
  class <- mlClass
  class[mlClassP < p.acpt] <- "uncertain"
  classification <- data.frame(class=class, mlClass=mlClass, mlClassP=mlClassP, mlClassType=mlClassType, stringsAsFactors=FALSE)
  classification <- cbind(classification, results)
  
  parameters=list(p.acpt=p.acpt, alpha=alpha, beta=beta, classifyTails=classifyTails, tol=tol, maxIter=maxIter, k.hto=k.hto, k.rna=k.rna)
  return(list(results=classification, model=mModel, parameters=parameters))
}


#' @importClassesFrom Matrix Matrix
setMethod("demuxmix",
          signature=c(object="Matrix", rna="missing"),
          .demuxmix)

setMethod("demuxmix",
          signature=c(object="matrix", rna="missing"),
          .demuxmix)

#' @importClassesFrom Matrix Matrix
setMethod("demuxmix",
          signature=c(object="Matrix", rna="numeric"),
          .demuxmix)

setMethod("demuxmix",
          signature=c(object="matrix", rna="numeric"),
          .demuxmix)