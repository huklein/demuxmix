.dmmClassify <- function(object) {
  p <- p.acpt(object)
  posteriorProb <- object@posteriorProb
  posHashtags <- posteriorProb >= 0.5
  posteriorProb[!posHashtags] <- 1 - posteriorProb[!posHashtags]
  hto <- apply(posHashtags, 2, function(h) {return(paste(rownames(posHashtags)[h], collapse=","))})
  hto[hto == ""] <- "negative"
  prob <- apply(posteriorProb, 2, prod)
  type <- c("negative", "singlet", "multiplet")[pmin(apply(posHashtags, 2, sum), 2) + 1]
  hto[prob < p] <- "uncertain"
  type[prob < p] <- "uncertain"
  return(data.frame(HTO=hto, Prob=prob, Type=type))
}


#' @importFrom methods setMethod
setMethod("dmmClassify", signature=c(object="Demuxmix"), .dmmClassify)