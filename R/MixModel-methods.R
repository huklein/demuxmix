# Internal methods for the internal classes NaiveMixModel and RegMixModel.
# These are mainly simple getter-methods so that higher higher-level methods
# do not have to distinguish between the two different types of mixture models.

# Methods for class NaiveMixModel
#' @importFrom methods setMethod
setMethod("getHto", signature=c(model="NaiveMixModel"),
  function(model, standardize=FALSE) {
    return(model@hto) # can't standardize for RNA in naive model
  }
)

setMethod("getMu1", signature=c(model="NaiveMixModel"),
  function(model, standardize=FALSE) {
    if (standardize) {
      return(model@mu1)
    } else {
      return(rep(model@mu1, times=length(getHto(model))))
    }
  }
)

setMethod("getMu2", signature=c(model="NaiveMixModel"),
  function(model, standardize=FALSE) {
    if (standardize) {
      return(model@mu2)
    } else {
      return(rep(model@mu2, times=length(getHto(model))))
    }
  }
)

setMethod("getTheta1", signature=c(model="NaiveMixModel"),
  function(model) {
    return(model@theta1)
  }
)

setMethod("getTheta2", signature=c(model="NaiveMixModel"),
  function(model) {
    return(model@theta2)
  }
)

setMethod("getPi", signature=c(model="NaiveMixModel"),
  function(model) {
    return(model@pi)
  }
)



# Methods for class RegMixModel
#' @importFrom stats predict residuals
#' @importFrom methods setMethod
setMethod("getHto", signature=c(model="RegMixModel"),
  function(model, standardize=FALSE) {
    if (standardize) { # standardize: return hto counts adjusted for different RNA profiles
      posteriorProb <- getPosteriorProbability(model)
      mu1 <- getMu1(model, standardize=TRUE)
      mu2 <- getMu2(model, standardize=TRUE)
      hto <- ((mu1 + residuals(model@fit1, type="response")) * posteriorProb[, 1]) +
             ((mu2 + residuals(model@fit2, type="response")) * posteriorProb[, 2])
      return(hto)
    } else {
      return(model@fit2$y)
    }
  }
)

setMethod("getMu1", signature=c(model="RegMixModel"),
  function(model, standardize=FALSE) {
    if (standardize) { # standardize: return mu for a cell with "average" RNA profile
      if (model@parameters$regRnaNegComp == FALSE) {
        mu <- predict(model@fit1, type="response")[1]
        names(mu) <- NULL
        return(mu)
      } else {
        posteriorProb <- getPosteriorProbability(model)
        rna <- sum(model@fit1$model$rna * posteriorProb[, 1]) / sum(posteriorProb[, 1])
        predData <- data.frame(hto=getHto(model), rna=rna)
        mu <- predict(model@fit1, newdata=predData, type="response")[1]
        names(mu) <- NULL
        return(mu)
      }
    } else {
      mu <- predict(model@fit1, type="response")
      names(mu) <- NULL
      return(mu)
    }
  }
)

setMethod("getMu2", signature=c(model="RegMixModel"),
  function(model, standardize=FALSE) {
    if (standardize) { # standardize: return mu for a cell with "average" RNA profile
      posteriorProb <- getPosteriorProbability(model)
      rna <- sum(model@fit2$model$rna * posteriorProb[, 2]) / sum(posteriorProb[, 2])
      predData <- data.frame(hto=getHto(model), rna=rna)
      mu <- predict(model@fit2, newdata=predData, type="response")[1]
      names(mu) <- NULL
      return(mu)
    } else {
      mu <- predict(model@fit2, type="response")
      names(mu) <- NULL
      return(mu)
    }
  }
)

setMethod("getTheta1", signature=c(model="RegMixModel"),
  function(model) {
    return(model@fit1$theta)
  }
)

setMethod("getTheta2", signature=c(model="RegMixModel"),
  function(model) {
    return(model@fit2$theta)
  }
)

setMethod("getPi", signature=c(model="RegMixModel"),
  function(model) {
    return(model@pi)
  }
)



# Methods for both classes Naive- and RegMixModel
#' @importFrom stats dnbinom
#' @importFrom methods setMethod
.getPosteriorProbability <- function(model) {
  hto <- getHto(model)
  f <- matrix(nrow=length(hto), ncol=2)
  f[, 1] <- dnbinom(hto, mu=getMu1(model), size=getTheta1(model))
  f[, 2] <- dnbinom(hto, mu=getMu2(model), size=getTheta2(model))
  pi.f <- t(getPi(model) * t(f))
  posteriorProb <- pi.f / apply(pi.f, 1, sum)
  return(posteriorProb)
}
setMethod("getPosteriorProbability", signature=c(model="NaiveMixModel"), .getPosteriorProbability)
setMethod("getPosteriorProbability", signature=c(model="RegMixModel"), .getPosteriorProbability)