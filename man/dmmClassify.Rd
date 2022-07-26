% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{dmmClassify}
\alias{dmmClassify}
\alias{dmmClassify,Demuxmix-method}
\title{Return classification results from a Demuxmix object}
\usage{
dmmClassify(object)
}
\arguments{
\item{object}{An object of class \code{\link{Demuxmix}}.}
}
\value{
A \code{data.frame} with 3 columns and one row for each droplet
  in the dataset. The first column gives the class (HTO) the droplet has been
  assigned to. The second column contains the posterior probability. And the
  third column specifies the type of the assigned class, i.e., "singlet",
  "multiplet", "negative" or "uncertain".
}
\description{
This method uses the posterior probabilities from the given demuxmix model
to assign each droplet to the most likely class, either a single HTO,
a combination of HTOs (multiplet) or the negative class (non-labeled
cells, empty droplets, cell debris). If the assignment cannot be made with
certainty above a defined threshold, the droplet is labeled as "uncertain".
}
\details{
A droplet is labeled as "uncertain" if the posterior probability of
  the most likely class is smaller than the threshold \code{pAcpt}, which
  is stored in the given \code{\link{Demuxmix}} object. The acceptance
  probability \code{pAcpt} can be inspected and set to a different value
  by applying the getter/setter method \code{\link{pAcpt}} to the
  \code{\link{Demuxmix}} object before calling this method. The method
  \code{\link{summary}} is useful to inspect classification results and to
  estimate error rates for different values of \code{pAcpt}.
}
\examples{
set.seed(2642)
simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
                                        c(rep(FALSE, 200), rep(TRUE, 220))))

dmm <- demuxmix(simdata$hto, rna = simdata$rna)
head(dmmClassify(dmm))
table(dmmClassify(dmm)$HTO, simdata$groundTruth)

pAcpt(dmm) <- 0.5
sum(dmmClassify(dmm)$HTO == "uncertain")
pAcpt(dmm) <- 0.9999
sum(dmmClassify(dmm)$HTO == "uncertain")

}
\seealso{
\code{\link{demuxmix}}
}
