#' @importFrom stats dnbinom predict
#' @importFrom ggplot2 ggplot geom_point geom_line scale_colour_gradient2 xlab ylab aes
.plotDmmScatter <- function(model, log=TRUE, pointsize=1.2, plotDecBoundary=TRUE, tol=0.01) {

  hto <- getHto(model, standardize=FALSE)
  rna <- exp(model@fit2$model$rna)
  posteriorProb <- getPosteriorProbability(model)
  
  if (log) {
    df <- data.frame(hto=log10(hto+1), rna=log10(rna+1), posteriorProb=posteriorProb[, 2])
    xlab <- substitute(paste(log[10], nn, sep=""), list(nn=paste("(", model@htoId, " HTO counts)", sep="")))
    ylab <- expression('log'[10]*'(RNA features detected)')
  
  } else {
    df <- data.frame(hto=hto, rna=rna, posteriorProb=posteriorProb[, 2])
    xlab <- paste(model@htoId, " HTO counts")
    ylab <- "RNA features detected"
  }

  p <- ggplot(df, aes(x=hto, y=rna, color=posteriorProb)) +
    geom_point(size=pointsize) +
    scale_colour_gradient2(name="Posterior\nprobability", low="dodgerblue", mid="black", high="firebrick2", midpoint=0.5) +
    xlab(xlab) +
    ylab(ylab)

  if (plotDecBoundary) {
    rnaRange <- c(min(rna), max(rna))
    htoRange <- c(min(hto), max(hto))
    db <- .getDecisionBoundary(model, rnaRange=rnaRange, htoRange=htoRange, tol=tol)
    if (log) {
      db$rna <- log10(db$rna + 1)
      db$hto <- log10(db$hto + 1)
    }
    p <- p + geom_line(aes(x=hto, y=rna, color=0.5), data=db, linetype="dashed") 
  }
  
  return(p)
}


# This methods numerically calculates the hto values with a posterior
# probability of 0.5 for 'npoints' within a given range 'rnaRange' of RNA
# values (number detected genes). hto values outside of the given 'htoRange'
# are removed.
.getDecisionBoundary <- function(model, rnaRange, htoRange, tol=0.01, npoints=10, maxIter=100) {
  
  rna <- seq(rnaRange[1], rnaRange[2], length.out=npoints)
  mu1 <- predict(model@fit1, newdata=data.frame(rna=log(rna)), type="response")
  mu2 <- predict(model@fit2, newdata=data.frame(rna=log(rna)), type="response")
  theta1 <- getTheta1(model)
  theta2 <- getTheta2(model)
  pi <- model@pi
  
  # htoRange is the (plotting) range of HTO values within we search for the
  # decision boundary. If the decision boundary is not withing the range for a
  # given RNA value, remove the RNA value.
  htoRange <- matrix(c(rep(htoRange[1], npoints), rep(htoRange[2], npoints)), ncol=2)
  lowerRange <- pi[2] * dnbinom(htoRange[, 1], mu=mu2, size=theta2) - pi[1] * dnbinom(htoRange[, 1], mu=mu1, size=theta1)
  upperRange <- pi[2] * dnbinom(htoRange[, 2], mu=mu2, size=theta2) - pi[1] * dnbinom(htoRange[, 2], mu=mu1, size=theta1)
  indOutside <- lowerRange >= 0 | upperRange <= 0
  htoRange <- htoRange[!indOutside, ]
  mu1 <- mu1[!indOutside]
  mu2 <- mu2[!indOutside]
  npoints <- sum(!indOutside)
  if (npoints == 0) {
    stop("Decision boundary outside plotting range.") # Should never happen.
  }
  
  granularity=10 # number of points within HTO range per iteration
  converged <- FALSE
  iter <- 1
  htoBoundary <- rep(NA, npoints)
  
  while (!converged) {
    
    for (j in 1:npoints) {
      # Calculate posterior prob - 0.5 for each point in the grid
      htos <- round(seq(htoRange[j, 1], htoRange[j, 2], length.out=granularity))
      comp1 <- pi[1] * dnbinom(htos, mu=mu1[j], size=theta1)
      comp2 <- pi[2] * dnbinom(htos, mu=mu2[j], size=theta2)
      ppValue <- comp2 / (comp1 + comp2) - 0.5
      
      # Find sign switch and define new hto boundaries
      indFirstPos <- min(which(ppValue > 0))
      htoRange[j, ] <- c(htos[indFirstPos - 1], htos[indFirstPos])
      
      # Interpolate boundary between hto limits
      relPos0 <- abs(ppValue[indFirstPos - 1]) / (abs(ppValue[indFirstPos - 1]) + ppValue[indFirstPos])
      htoBoundary[j] <- htos[indFirstPos - 1] + relPos0 * (htos[indFirstPos] - htos[indFirstPos - 1])
    }
    
    # Check convergence criteria
    # 1. All HTO ranges are reduced to just 1 HTO read
    if (all(htoRange[, 2] - htoRange[, 1] <= 1)) {
      converged <- TRUE
    }
    # 2. Error is smaller than given tol
    comp1 <- pi[1] * dnbinom(round(htoBoundary), mu=mu1, size=theta1)
    comp2 <- pi[2] * dnbinom(round(htoBoundary), mu=mu2, size=theta2)
    error <- comp2 / (comp1 + comp2) - 0.5
    if (all(abs(error) <= tol)) {
      converged <- TRUE
    }
    # 3. Max iterations reached
    if (iter == maxIter) {
      converged <- TRUE
      warning("Estimation of decision boundary did not converge.")
    }
    
    iter <- iter + 1
    
  }
  
  decisionBoundary <- data.frame(rna=rna, hto=htoBoundary, error=error)
  return(decisionBoundary)
}


#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
#' @importFrom methods is
setMethod("plotDmmScatter", signature=c(object="Demuxmix", hto="missing"),
  function (object, hto, log=TRUE, pointsize=1.2, plotDecBoundary=TRUE, tol=0.01) {
    ind <- sapply(object@models, is, "RegMixModel")
    if (sum(ind) == 0) {
      stop("The given object does not contain regression mixture models.")
    } else if (any(!ind)) {
      warning(paste("Plots are only generated for regression mixture models (",
              paste(names(ind)[ind], collapse=", "), ").", sep=""))
    }
    plots <- lapply(object@models[ind], .plotDmmScatter, log=log, pointsize=pointsize, plotDecBoundary=plotDecBoundary, tol=tol)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)

setMethod("plotDmmScatter", signature=c(object="Demuxmix", hto="ANY"),
  function (object, hto, log=TRUE, pointsize=1.2, plotDecBoundary=TRUE, tol=0.01) {
    if (is.numeric(hto) & any(hto > length(object@models))) {
      stop("Invalid HTO identifier.")
    }
    if (is.character(hto) & any(!is.element(hto, names(object@models)))) {
      stop("Invalid HTO identifier.")
    }
    models <- object@models[hto]
    ind <- sapply(models, is, "RegMixModel")
    if (sum(ind) == 0) {
      stop("The specified models are not regression mixture models.")
    } else if (any(!ind)) {
      warning(paste("Plots are only generated for regression mixture models (",
                    paste(names(ind)[ind], collapse=", "), ").", sep=""))
    }
    plots <- lapply(models[ind], .plotDmmScatter, log=log, pointsize=pointsize, plotDecBoundary=plotDecBoundary, tol=tol)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)