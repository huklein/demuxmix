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
          type <- "RegMixModel (not negative comp.)"
      }
      if (object@models[[i]]@converged) {
        converged <- "converged"
      } else {
        converged <- "not converged!"
      }
      cat(paste("  ", htos[i], ": ", type, "; ", converged, "\n", sep=""))
    }
  }
)


setMethod("getPosteriorProbability", signature=c(model="Demuxmix"),
  function(model) {
    return(model@posteriorProb)
  }
)

setMethod("p.acpt", signature=c(model="Demuxmix"),
  function(model) {
    return(model@parameters$p.acpt)
  }
)

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

setMethod("dmmClassify", signature=c(model="Demuxmix"),
  function(model) {
    p <- p.acpt(model)
    posteriorProb <- getPosteriorProbability(model)
    posHashtags <- posteriorProb >= 0.5
    posteriorProb[!posHashtags] <- 1 - posteriorProb[!posHashtags]
    hto <- apply(posHashtags, 2, function(h) {return(paste(rownames(posHashtags)[h], collapse=","))})
    hto[hto == ""] <- "negative"
    prob <- apply(posteriorProb, 2, prod)
    type <- c("negative", "singleton", "multiplet")[pmin(apply(posHashtags, 2, sum), 2) + 1]
    hto[prob < p] <- "uncertain"
    return(data.frame(HTO=hto, Prob=prob, Type=type))
  }
)


# setGeneric("summary") has been set in AllGenerics
summary.Demuxmix <- function(object, ...) {
  dmmResults <- dmmClassification(object)
  
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

setMethod("summary", signature=c(object="Demuxmix"), summary.Demuxmix)
