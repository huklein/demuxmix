#' @describeIn Demuxmix-class Displays the object on the command line.
setMethod("show", signature=c(object="Demuxmix"),
  function(object) {
    numCells <- ncol(object@posteriorProb)
    htos <- rownames(object@posteriorProb)
    if (length(htos) == 1) {
      cat(paste("Demuxmix object with", length(htos), "HTO and", numCells, "cells.\n"))
    } else {
      cat(paste("Demuxmix object with", length(htos), "HTOs and", numCells, "cells.\n"))
    }
    for (i in 1:length(htos)) {
      type <- class(object@models[[i]])
      if (type == "RegMixModel" && !object@models[[i]]@parameters$regRnaNegComp) {
          type <- "RegMixModel (no regression for negative comp.)"
      }
      if (object@models[[i]]@converged) {
        converged <- "converged"
      } else {
        converged <- "not converged!"
      }
      cat(paste("  ", htos[i], ": ", type, "; ", converged, "\n", sep=""))
      cat(paste(paste(rep(" ", times=nchar(htos[i])+4), collapse=""),
                "mu=(", format(getMu1(object@models[[i]], standardize=TRUE), digits=5), ", ",
                format(getMu2(object@models[[i]], standardize=TRUE), digits=5), ")\n", sep=""))
      if (type == "RegMixModel") {
        cat(paste(paste(rep(" ", times=nchar(htos[i])+4), collapse=""),
                  "RNA coef. negative comp: ",
                  format(coefficients(summary(object@models[[i]]@fit1))["rna", "Estimate"], digits=4),
                  ", p=", format(coefficients(summary(object@models[[i]]@fit1))["rna", "Pr(>|z|)"], digits=3), "\n", sep=""))
        cat(paste(paste(rep(" ", times=nchar(htos[i])+4), collapse=""),
                  "          positive comp: ",
                  format(coefficients(summary(object@models[[i]]@fit2))["rna", "Estimate"], digits=4),
                  ", p=", format(coefficients(summary(object@models[[i]]@fit2))["rna", "Pr(>|z|)"], digits=3), "\n", sep=""))
      }
      if (type == "RegMixModel (no regression for negative comp.)") {
        cat(paste(paste(rep(" ", times=nchar(htos[i])+4), collapse=""),
                  "RNA coef. positive comp: ",
                  format(coefficients(summary(object@models[[i]]@fit2))["rna", "Estimate"], digits=4),
                  ", p=", format(coefficients(summary(object@models[[i]]@fit2))["rna", "Pr(>|z|)"], digits=3), "\n", sep=""))
      }
      
    }
  }
)



# setMethod("getPosteriorProbability", signature=c(model="Demuxmix"),
#   function(model) {
#     return(model@posteriorProb)
#   }
# )


#' @describeIn Demuxmix-class Returns the acceptance probability \code{p.acpt}.
setMethod("p.acpt", signature=c(model="Demuxmix"),
  function(model) {
    return(model@parameters$p.acpt)
  }
)


#' @describeIn Demuxmix-class Sets a new acceptance probability \code{p.acpt}.
setMethod("p.acpt<-", signature=c(model="Demuxmix", value="numeric"),
  function(model, value) {
    if (length(value) > 1) {
      value <- value[1]
      warning("Only first element will be used.")
    }
    if (value < 0 | value > 1) {
      stop("New acceptance probability (p.acpt) must be in [0, 1].")
    }
    model@parameters$p.acpt <- value
    return(model)
  }
)


# setGeneric("summary") has been set in AllGenerics
summary.Demuxmix <- function(object, ...) {
  dmmResults <- dmmClassify(object)
  
  # Label all multiplets as "multiplet"
  dmmResults$HTO[dmmResults$Type == "multiplet" & dmmResults$HTO != "uncertain"] <- "multiplet"
  
  # Get all classes and sort them
  otherClasses <- c("singleton", "multiplet", "negative", "uncertain")
  htoClasses <- sort(unique(dmmResults$HTO))
  htoClasses <- htoClasses[! htoClasses %in% otherClasses]
  
  # Remove uncertain droplets
  nUncertain <- sum(dmmResults$HTO == "uncertain")
  dmmResults <- dmmResults[dmmResults$HTO != "uncertain", ]
  
  # Summarize classification results
  classSummary <- data.frame(Class=c(htoClasses, otherClasses), NumObs=NA, RelFreq=NA, MedProb=NA, ExpFPs=NA, FDR=NA)
  for (class in c(htoClasses, c("multiplet", "negative"))) {
    ind <- dmmResults$HTO == class
    classSummary$NumObs[classSummary$Class == class] <- sum(ind)
    classSummary$RelFreq[classSummary$Class == class] <- sum(ind) / length(ind)
    classSummary$MedProb[classSummary$Class == class] <- median(dmmResults$Prob[ind])
    classSummary$ExpFPs[classSummary$Class == class] <- sum(1 - dmmResults$Prob[ind])
    classSummary$FDR[classSummary$Class == class] <- classSummary$ExpFPs[classSummary$Class == class] / classSummary$NumObs[classSummary$Class == class]
  }
  # all singletons combined
  ind <- dmmResults$HTO %in% htoClasses
  classSummary$NumObs[classSummary$Class == "singleton"] <- sum(ind)
  classSummary$RelFreq[classSummary$Class == "singleton"] <- sum(ind) / length(ind)
  classSummary$MedProb[classSummary$Class == "singleton"] <- median(dmmResults$Prob[ind])
  classSummary$ExpFPs[classSummary$Class == "singleton"] <- sum(1 - dmmResults$Prob[ind])
  classSummary$FDR[classSummary$Class == "singleton"] <- classSummary$ExpFPs[classSummary$Class == "singleton"] / classSummary$NumObs[classSummary$Class == "singleton"]
  # uncertain
  classSummary$NumObs[classSummary$Class == "uncertain"] <- nUncertain
  return(classSummary)
}

#' @describeIn Demuxmix-class Summarizes the classification results and
#'   estimates error rates.
setMethod("summary", signature=c(object="Demuxmix"), summary.Demuxmix)