.dmmSummary <- function(dmmResults, p.acpt) {
   stopifnot(all(c("class", "mlClassP", "mlClassType") %in% colnames(dmmResults)))
   
   # Re-classify uncertain droplets with new p.acpt if given
   if (!missing(p.acpt)) {
      dmmResults$class <- dmmResults$mlClass
      dmmResults$class[dmmResults$mlClassP < p.acpt] <- "uncertain"
   }
   
   # Label all multiplets as "multiplet"
   dmmResults$class[dmmResults$mlClassType == "multiplet" & dmmResults$class != "uncertain"] <- "multiplet"
   
   # Get all classes and sort them
   otherClasses <- c("singleton", "multiplet", "negative", "uncertain")
   htoClasses <- sort(unique(dmmResults$class))
   htoClasses <- htoClasses[! htoClasses %in% otherClasses]
   
   # Remove uncertain from results
   nUncertain <- sum(dmmResults$class == "uncertain")
   dmmResults <- dmmResults[dmmResults$class != "uncertain", ]
   
   # Summarize classification results
   classSummary <- data.frame(Class=c(htoClasses, otherClasses), NumObs=NA, RelFreq=NA, MedProb=NA, ExpFPs=NA, FDR=NA,
                              stringsAsFactors=FALSE)
   for (class in c(htoClasses, c("multiplet", "negative"))) {
      ind <- dmmResults$class == class
      classSummary$NumObs[classSummary$Class == class] <- sum(ind)
      classSummary$RelFreq[classSummary$Class == class] <- sum(ind) / length(ind)
      classSummary$MedProb[classSummary$Class == class] <- median(dmmResults$mlClassP[ind])
      classSummary$ExpFPs[classSummary$Class == class] <- sum(1 - dmmResults$mlClassP[ind])
      classSummary$FDR[classSummary$Class == class] <- classSummary$ExpFPs[classSummary$Class == class] / classSummary$NumObs[classSummary$Class == class]
   }
   # all singletons
   ind <- dmmResults$class %in% htoClasses
   classSummary$NumObs[classSummary$Class == "singleton"] <- sum(ind)
   classSummary$RelFreq[classSummary$Class == "singleton"] <- sum(ind) / length(ind)
   classSummary$MedProb[classSummary$Class == "singleton"] <- median(dmmResults$mlClassP[ind])
   classSummary$ExpFPs[classSummary$Class == "singleton"] <- sum(1 - dmmResults$mlClassP[ind])
   classSummary$FDR[classSummary$Class == "singleton"] <- classSummary$ExpFPs[classSummary$Class == "singleton"] / classSummary$NumObs[classSummary$Class == "singleton"]
   # uncertain
   classSummary$NumObs[classSummary$Class == "uncertain"] <- nUncertain
      
   return(classSummary)
}


setMethod("dmmSummary",
  signature=c(dmmResults="data.frame"),
  .dmmSummary)
