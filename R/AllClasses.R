#
# ADD DOCUMENTATION
#
Demuxmix <- setClass("Demuxmix",
  slots=c(models="list",
          outliers="matrix",
          clusterInit="matrix",
          posteriorProb="matrix",
          tailException="matrix",
          modelSelection="data.frame",
          parameters="list"),
  validity=function(object) {
    if(!all(sapply(object@models, is, "NaiveMixModel") |
            sapply(object@models, is, "RegMixModel"))) {
      return("models must contain objects of class NaivMixModel or RegMixModel.")
    }
    if (any(is.na(object@posteriorProb))) {
      return("posteriorProb must not contain missing values.")
    }
    if (any(object@posteriorProb < 0) | any(object@posteriorProb > 1)) {
      return("All values in posteriorProb must be in [0, 1].")
    }
    if (!is.element("p.acpt", names(object@parameters))) {
      return("p.acpt must be element in parameters.")
    }
    if (object@parameters$p.acpt < 0 | object@parameters$p.acpt > 1) {
      return("p.acpt given in parameters must be in [0, 1].")
    }
    return(TRUE)
  }
)


# Internal class representing a naive mixture model
# fitted to a single hastag oligo.
NaiveMixModel <- setClass("NaiveMixModel",
  slots=c(mu1="numeric",
          mu2="numeric",
          theta1="numeric",
          theta2="numeric",
          pi="numeric",
          hto="numeric",
          htoId="character",
          parameters="list",
          log="data.frame",
          converged="logical"),
  validity=function(object) {
    if (length(object@mu1) > 1 |  length(object@mu2) > 1 |
        length(object@theta1) > 1 | length(object@theta1) > 1) {
      return("mu and theta must have length 1.")
    }
    if (object@mu1 > object@mu2) {
      return("mu1 must be smaller than mu2.")
    }
    if (object@theta1 <= 0 | object@theta2 <= 0) {
      return("theta1 and theta2 must be larger than 0.")
    }
    if (length(object@pi) != 2) {
      return("pi must be of length 2.")
    }
    if (any(object@pi < 0) | any(object@pi > 1)) {
      return("pi must be 0 <= pi <= 1.")
    }
    if (any(object@hto < 0)) {
      return("hto contain positive numbers.")
    }
    return(TRUE)
  }
)


# Internal class representing a regression mixture model
# fitted to a single hastag oligo.
RegMixModel <- setClass("RegMixModel",
  slots=c(fit1="ANY",
          fit2="ANY",
          pi="numeric",
          htoId="character",
          parameters="list",
          log="data.frame",
          converged="logical"),
          validity=function(object) {
            if (length(object@fit1$y) != length(object@fit2$y)) {
              return("Both models must be fitted to the same data.")
            }
            if (! "rna" %in% colnames(object@fit2$model)) {
              return("rna must be a covariate of the second model.")
            }
            if (length(object@pi) !=2) {
              return("pi must be of length 2.")
            }
            if (any(object@pi < 0) | any(object@pi > 1)) {
              return("pi must be 0 <= pi <= 1.")
            }
            return(TRUE)
          }
)