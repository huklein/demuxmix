#' @importFrom stats dnbinom pnbinom qnbinom
.dmmOverlap <- function(model, tol=0.001) {

  mu1 <- getMu1(model, standardize=TRUE)
  mu2 <- getMu2(model, standardize=TRUE)
  theta1 <- getTheta1(model)
  theta2 <- getTheta2(model)
  
  q1 <- qnbinom(1-tol, mu=mu1, size=theta1)
  q2 <- qnbinom(1-tol, mu=mu2, size=theta2)
  maxq <- max(q1, q2)
  area <- sum(pmin(dnbinom(0:maxq, mu=mu1, size=theta1), dnbinom(0:maxq, mu=mu2, size=theta2)))
  
  if (q1 <= q2) {  # density one has smaller tail at max(q1, q2)
    tail <- pnbinom(q=maxq, mu=mu1, size=theta1, lower.tail=FALSE)
  } else {         # density two has smaller tail at max(q1, q2)
    tail <- pnbinom(q=maxq, mu=mu2, size=theta2, lower.tail=FALSE)
  }
    
  return(area + tail)
}


#' @importFrom methods setMethod
setMethod("dmmOverlap", signature=c(object="Demuxmix", hto="missing"),
  function (object, hto, tol=0.001) {
    return(vapply(object@models, .dmmOverlap, numeric(1), tol=tol))
  }
)

setMethod("dmmOverlap", signature=c(object="Demuxmix", hto="ANY"),
  function (object, hto, tol=0.001) {
    if (is.numeric(hto) & any(hto > length(object@models))) {
      stop("Invalid HTO identifier.")
    }
    if (is.character(hto) & any(!is.element(hto, names(object@models)))) {
      stop("Invalid HTO identifier.")
    }
    return(vapply(object@models[hto], .dmmOverlap, numeric(1), tol=tol))
  }
)