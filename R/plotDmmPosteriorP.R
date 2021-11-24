#' @importFrom ggplot2 ggplot geom_histogram xlab ylab aes
.plotDmmPosteriorP <- function(model, bins=50) {
  
  posteriorProb <- getPosteriorProbability(model)
  df <- data.frame(posteriorProb=posteriorProb[, 2])
  
  p <- ggplot(df, aes(x=posteriorProb)) +
    geom_histogram(bins=bins) +
    xlab("Posterior probability") +
    ylab("Number of cells")

  return(p)
}


#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
setMethod("plotDmmPosteriorP", signature=c(object="Demuxmix", hto="missing"),
  function (object, hto, bins=50) {
    plots <- lapply(object@models, .plotDmmPosteriorP, bins=bins)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)

setMethod("plotDmmPosteriorP", signature=c(object="Demuxmix", hto="ANY"),
  function (object, hto, bins=50) {
    if (is.numeric(hto) & any(hto > length(object@models))) {
      stop("Invalid HTO identifier.")
    }
    if (is.character(hto) & any(!is.element(hto, names(object@models)))) {
      stop("Invalid HTO identifier.")
    }
    plots <- lapply(object@models[hto], .plotDmmPosteriorP, bins=bins)
    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(do.call(grid.arrange, plots))
    }
  }
)