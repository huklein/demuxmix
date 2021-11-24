#' @importFrom stats dnbinom predict
#' @importFrom ggplot2 ggplot geom_point scale_colour_gradient2 xlab ylab aes
.plotDmmScatter <- function(model, log=TRUE, pointsize=1.2) {

  hto <- getHto(model, standardize=FALSE)
  rna <- model@fit2$model$rna
  posteriorProb <- getPosteriorProbability(model)
  
  if (log) {
    df <- data.frame(hto=log10(hto+1), rna=log10(rna+1), posteriorProb=posteriorProb[, 2])
    xlab <- substitute(paste(log[10], nn, sep=""), list(nn=paste("(", model@htoId, " counts)", sep="")))
    ylab <- expression('log'[10]*'(RNA features detected)')
  
  } else {
    df <- data.frame(hto=hto, rna=rna, posteriorProb=posteriorProb[, 2])
    xlab <- paste(model@htoId, "count")
    ylab <- "RNA features detected"
  }
  
  p <- ggplot(df, aes(x=hto, y=rna, color=posteriorProb)) +
    geom_point(size=pointsize) +
    scale_colour_gradient2(name="Posterior\nprobability", low="dodgerblue", mid="black", high="firebrick2", midpoint=0.5) +
    xlab(xlab) +
    ylab(ylab)
  
  return(p)
}


#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
setMethod("plotDmmScatter", signature=c(object="Demuxmix", hto="missing"),
  function (object, hto, log=TRUE, pointsize=1.2) {
    ind <- sapply(object@models, isa, "RegMixModel")
    if (sum(ind) == 0) {
      stop("The given object does not contain regression mixture models.")
    } else if (any(!ind)) {
      warning(paste("Plots are only generated for regression mixture models (",
              paste(names(ind)[ind], collapse=", "), ").", sep=""))
    }
    plots <- lapply(object@models[ind], .plotDmmScatter, log=log, pointsize=pointsize)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)

setMethod("plotDmmScatter", signature=c(object="Demuxmix", hto="ANY"),
  function (object, hto, log=TRUE, pointsize=1.2) {
    if (is.numeric(hto) & any(hto > length(object@models))) {
      stop("Invalid HTO identifier.")
    }
    if (is.character(hto) & any(!is.element(hto, names(object@models)))) {
      stop("Invalid HTO identifier.")
    }
    models <- object@models[hto]
    ind <- sapply(models, isa, "RegMixModel")
    if (sum(ind) == 0) {
      stop("The specified models are not regression mixture models.")
    } else if (any(!ind)) {
      warning(paste("Plots are only generated for regression mixture models (",
                    paste(names(ind)[ind], collapse=", "), ").", sep=""))
    }
    plots <- lapply(models[ind], .plotDmmScatter, log=log, pointsize=pointsize)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)