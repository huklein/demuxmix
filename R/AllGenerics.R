#' Demultiplexing using mixture models
#'
#' This methods uses mixture models as a probabilistic framework to assign
#' hashtags to cells and to identify multiplets based on counts obtained from
#' sequencing a hastag oligo (HTO) library. If the numbers of detected genes
#' from the corresponding RNA library are passed in the second argument,
#' regression mixture models are used, which usually improves the
#' classification accuracy by leveraging the relationship between HTO and RNA
#' counts.
#' 
#' @param object A matrix of HTO counts where each row corresponds to a hashtag
#'   and each column to a cell barcode.
#' @param rna An optional numeric vector with the number of genes detected in
#'   the RNA library for each cell. Same length as columns in \code{object}.
#'   If missing, a naive instead of a regression mixture model is used.
#' @param p.acpt Acceptance probability that must be reached in order to
#'   assign a cell to a hashtag. Cells with lower probabilities are classified
#'   as "uncertain".
#' @param alpha Threshold for defining the left tail of the mixture
#'   distribution where cells should not be classified as "positive". Either a
#'   single value between 0 and 1 or a vector with one entry for each hashtag
#'   in the dataset. See details.
#' @param beta Threshold for defining the right tail of the mixture
#'   distribution where cells should not be classified as "negative". Either a
#'   single value between 0 and 1 or a vector with one entry for each hashtag
#'   in the dataset. See details.
#' @param correctTails If \code{TRUE}, cells meeting the threshold defined by
#'   \code{alpha} (\code{beta}) are classified as "negative" ("positive") even
#'   if the mixture model suggests a different classification. See details.
#' @param tol Convergence criterion for the EM algorithm fitting the mixture
#'   model(s). The algorithm stops when the relative increase of the log
#'   likelihood is less than or equal to \code{tol}. Either a single value
#'   or a vector with one entry for each hashtag in the dataset.
#' @param maxIter Maximum number of iterations for the EM algorithm and
#'   for the alternating iteration process fitting the nb regression models
#'   within each EM iteration. Either a single value or a vector with one
#'   entry for each hashtag in the dataset.
#' @param k.hto Factor to define outliers in the HTO counts. Among cells
#'   positive for the hashtag based on initial clustering, HTO counts
#'   larger than the 0.75 quantile + \code{k.hto} * IQR are considered
#'   outliers. Either a single value or a vector with one entry for each
#'   hashtag in the dataset. See details.
#' @param k.rna Factor to define outliers in the numbers of detected genes.
#'   Numbers of detected genes larger than the 0.75 quantile +
#'   \code{k.rna} * IQR are considered outliers. Either a single value or
#'   a vector with one entry for each hashtag in the dataset. See details.
#' 
#' @details The single cell dataset should undergo basic filtering to
#'   remove low quality or empty droplets before calling this function,
#'   but the HTO counts should not be transformed or pre-processed otherwise.
#'   The number of detected genes passed via the optional argument \code{rna}
#'   is typically defined as the number of genes in the RNA library with at
#'   least one read.
#'   
#'   The method itself consists of three steps, which can be modified by the
#'   remaining parameters. The default settings work well for a wide
#'   range of datasets, and in most cases only the acceptance probability
#'   \code{p.acpt} might be of interest (step 3). Steps 1 and 2 are
#'   executed for each hashtag separately using only the hto counts for that
#'   hashtag, step 3 classifies the cells based on results from all hashtags:
#'   
#'   \enumerate{
#'     \item Preprocessing (\code{k.hto}, \code{k.rna}). Cells are clustered
#'     into a negative and a positive group based on the hto counts of the
#'     hashtag using k-means. Cells in the positive group with hto counts
#'     larger than \code{k.hto} times the IQR of the hto counts in the positive
#'     group are marked as outliers. Outliers are still classified but will
#'     not be used to fit the mixture model for this hashtag in step 2. If
#'     the parameter \code{rna} is given, all cells (both groups) with number
#'     of detected genes larger than k.rna times the IQR are marked as
#'     outliers, too. Outliers with very large read counts are often
#'     observed in sequencing data and can negatively affect model fitting in
#'     step 2. If more than 3\% of the cells are marked as outliers, a warning
#'     message is printed and larger values for \code{k.hto} and \code{k.rna}
#'     might be preferable. If the model fit seems to be affected by a few
#'     large values (very high variance of the positive component), smaller
#'     values should be chosen.
#'     
#'     \item Model fitting (\code{alpha}, \code{beta}, \code{correctTails},
#'     \code{tol}, \code{maxIter}). An EM algorithm is used to fit a mixture
#'     model to the cells not marked as outliers in step 1. \code{maxIter}
#'     defines the maximum number of iterations of the EM algorithm, and,
#'     if parameter \code{rna} is given, it also defines the maximum number
#'     of iterations to fit the negative binomial regression models within
#'     each EM iteration. \code{tol} defines the convergence criterion for
#'     the EM algorithm. The algorithm stops if \eqn{\DeltaLL/LL \le}
#'     \code{tol}. After the mixture model is fit, the posterior probability
#'     that cell X is positive for the hashtag \eqn{P(X = pos)} is
#'     calculated. Depending on the given data, these probabilities can be
#'     inaccurate at the far tails of the mixture distribution. Specifically,
#'     a positive component with large variance can be larger than the
#'     negative component close to zero, if the negative component is
#'     narrow and shifted to the right to model ambient hashtags.
#'     If \code{correctTails} is \code{TRUE}, the following two rules are
#'     applied to avoid false classifications at the at the far tails. First,
#'     if a cell X is classified as positive based on the posterior
#'     probability, but the probability to detected more than the observed
#'     \eqn{h_x} hashtags in a negative cell is \eqn{P(H \ge h_x | neg)} >
#'     \code{alpha}, then \eqn{P(X = pos)} is set to 0 (left tail). Second,
#'     if a cell X is classified as negative, but \eqn{P(H \le h_x | pos)} >
#'     \code{beta}, \eqn{P(X = pos)} is set to 1 (right tail). For most
#'     datasets, these rules will not apply and it's recommend not to change
#'     these values. If \code{correctTails} is \code{FALSE}, posterior
#'     probabilities will not be altered, but potential problems at the tails
#'     will still be logged. See Value.
#'     
#'     \item Classification (\code{p.acpt}). The posterior probabilities
#'     obtained from the models fitted to each hashtag separately are 
#'     used to calculate the most likely class for each cell. The following
#'     classes are considered: one class for each hashtag (singletons), one
#'     class for each possible multiplet, and a negative class representing
#'     droplets negative for all hashtags (i.e. empty droplets or droplets
#'     containing only parts of cells). Each cell is assigned to the most
#'     likely class unless the probability is smaller than /code{p.acpt},
#'     in which case the cell is assigned to the class "uncertain". The
#'     calculation of the posterior probabilities assumes that
#'     the HTO counts of the different hashtags are distributed
#'     independently. This assumption is unlikely for real HTO data
#'     where underlying factors such as cell size usually cause a positive
#'     correlation. Consequently, a mild overestimation of multiplets and
#'     negative droplets can be expected, whereas the estimated error rates for
#'     singletons, which are selected for downstream analyses, are conservative
#'     and are likely an upper boundary.
#'   }
#'   
#' @return \code{demuxmix} returns a list with four elements. The first
#'   element contains the classification results in a \code{data.frame}. The
#'   other elements are helpful for assessing the model fit, but not intended
#'   to be directly used in most cases. In detail, the following elements
#'   are returned:
#'   
#'   \describe{
#'     \item{results}{A \code{data.frame} with four columns and one row
#'       per cell. The columns "mlClass", "mlClassP", and "mlClassType"
#'       contain the most likely class, the probability of the most likely
#'       class, and its type (singleton, multiplet, or negative). The column
#'       "class" contains the final classification and is either equal to
#'       "mlClass" or set to "uncertain", if "mlClassP" is smaller than the
#'       given threshold \code{p.acpt}. Classification results for a different
#'       threshold \code{p.acpt} can be quickly obtained using "mlClass" and
#'       "mlClassP" instead of re-running the method.}
#'   
#'     \item{htoResults}{A \code{data.frame} with four columns for each
#'       hashtag in the dataset and one row per cell storing information about
#'       the separate mixture models for each hashtag. The four columns for
#'       each hashtag contain the posterior probability that the cell is
#'       positive for the hashtag, the probability \eqn{P(H \ge h_x | neg)},
#'       the probability \eqn{P(H \le h_x | pos)}, and a \code{logical}
#'       value indicating whether one of the two rules indicating potential
#'       false classification at the tails was met. The logical value is also
#'       set to \code{TRUE} if \code{correctTails} is \code{FALSE}. See
#'       Details for more information about the last three columns.}
#'       
#'     \item{model}{A \code{list} with one element for each hashtag. Each list
#'       element stores the fitted mixture model together with the data used
#'       to fit the model and convergence information. In case of a regression
#'       mixture model, the list contains objects of class \code{negbin}. The
#'       list "model" or any of its elements can be passed to one of the
#'       plotting function, e.g. \code{\link{plotDmmHistogram}}, to assess
#'       the fit of the model.}
#'       
#'     \item{parameters}{A \code{list} storing all parameters passed to 
#'      \code{demuxmix}.}
#'  }
#'  
#' @seealso \code{\link{dmmSummary}} to summarize the classification results.
#'   \code{\link{plotDmmHistogram}}, \code{\link{plotDmmScatter}},
#'   \code{\link{plotDmmPosteriorP}}, and \code{\link{dmmOverlap}} to assess
#'   the model fit.
#'   
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, p.acpt=0.9)
#' table(dmm$results$class, simdata$groundTruth)
#' 
#' dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' table(dmmreg$results$class, simdata$groundTruth)
#' dmmSummary(dmmreg$results)
#' dmmOverlap(dmmreg$model)
#' \donttest{
#' plotDmmHistogram(dmmreg$model)
#' plotDmmScatter(dmmreg$model[[1]])}
#' 
#' @aliases demuxmix,matrix,missing-method demuxmix,matrix,numeric-method
#'   demuxmix,Matrix,missing-method demuxmix,Matrix,numeric-method
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("demuxmix",
           function(hto, rna, p.acpt=0.9^nrow(hto), model="auto", alpha=0.9, beta=0.9, correctTails=TRUE, tol=10^-5, maxIter=100, k.hto=1.5, k.rna=1.5)
             standardGeneric("demuxmix"),
           signature=c("hto", "rna"))



#' Calculate area intersected by two components of a mixture model
#' 
#' \code{dmmOverlap} integrates over the area intersected by the two
#' components of the given mixture model. The integral should be
#' close to 0 if the HTO labeling experiment was successful.
#'
#' @param object A mixture model or a list of mixture models as
#'   returned by \code{\link{demuxmix}}.
#' @param tol The maximum acceptable error when calculating the
#'   area.
#'   
#' @details The area under both the negative and positive component is an
#'   informative quality metric for the hashtag labeling efficiency. Values
#'   under 0.03 can be considered as good, values larger than 0.1 are
#'   problematic.
#'   
#'   The definition of the area is not obvious for a regression mixture
#'   model since the distributions' means depend on the covariate, i.e.,
#'   the number of detected genes in the RNA library. This method 
#'   calculates the weighted mean number of detected genes in the cells
#'   for each component, which are then used to calculate the expected
#'   number of HTO counts for the negative and positive component.
#' 
#' @return A numeric vector with the intersection area for each given
#'   mixture model.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples 
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, p.acpt=0.9)
#' dmmOverlap(dmm$model)
#' 
#' dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' dmmOverlap(dmmreg$model)
#' dmmOverlap(dmmreg$model[["HTO_1"]])
#' 
#' @aliases dmmOverlap,list-method
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("dmmOverlap",
           function(object, hto, tol=0.001)
             standardGeneric("dmmOverlap"),
           signature=c("object", "hto"))



#' Simulate HTO sequencing data
#' 
#' This method simulates HTO count data and corresponding numbers of detected
#' RNA features using the negative binomial distribution. The purpose of this
#' method is to provide simple example datasets for testing and documentation.
#'
#' @param class A \code{matrix} of type logical defining the number of
#'   hashtags, the number of cells, and the cells' class memberships, i.e.,
#'   which cells have been tagged with which hashtag. Each row corresponds
#'   to one hashtag and each column to a cell. Negative cells (all entries
#'   in the column are \code{FALSE}) and multiplets (more than one entry are
#'   \code{TRUE}) are allowed. If the matrix has row names, the names must be
#'   unique and are used as hashtag names.
#' @param mu Vector of expectation values of the HTO counts if a
#'   cell is positive for the hashtag. Values are recycled if \code{mu} is
#'   shorter than number of hashtags defined by \code{class}.
#' @param theta Vector of dispersion parameters of the HTO counts if a
#'   cell is positive for the hashtag. Values are recycled if \code{theta} is
#'   shorter than number of hashtags defined by \code{class}.
#' @param muAmbient Vector of expectation values of the HTO counts if a
#'   cell is negative for the hashtag. Values are recycled if \code{mu} is
#'   shorter than number of hashtags defined by \code{class}.
#' @param thetaAmbient Vector of dispersion parameters of the HTO counts if a
#'   cell is negative for the hashtag. Values are recycled if \code{theta} is
#'   shorter than number of hashtags defined by \code{class}.
#' @param muRna Single expectation value for the number of detected RNA
#'   features.
#' @param thetaRna Single dispersion parameter for the number of detected RNA
#'   features.
#'
#' @details A vector \eqn{r} of detected RNA features (same length as columns
#'   in \code{class}) is simulated using \code{\link[stats]{rnbinom}} with
#'   \code{muRna} and \code{thetaRna} as parameters. HTO counts of positive
#'   cells are then simulated using \code{\link[stats]{rnbinom}} with 
#'   \eqn{r} \code{mu}/\code{muRna} as expectation value and \code{theta}
#'   as dispersion. If a cell is negative for the hastag,
#'   \eqn{r} \code{muAmbient}/\code{muRna} and \code{thetaAmbient}
#'   are used respectively.
#' 
#' @return A list with three elements: "hto" is a matrix of same dimension as
#'   the given \code{class} matrix and contains the simulated HTO counts.
#'   "rna" is a vector of simulated detected number of genes (same length as
#'   "hto" has columns). "groundTruth" is a character vector encoding the
#'   class labels given by \code{class} as character strings for convenience.
#'
#' @seealso \code{\link{demuxmix}} 
#' 
#' @examples
#' set.seed(2642)
#' class <- rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                c(rep(FALSE, 200), rep(TRUE, 220)))
#' simdata <- dmmSimulateHto(class=class, mu=c(150, 300), theta=c(15, 20),
#'                           muAmbient=c(30, 30), thetaAmbient=c(10, 10),
#'                           muRna=3000, thetaRna=30)
#' dim(simdata$hto)
#' table(simdata$groundTruth)
#' 
#' mean(simdata$rna) # muRna
#' var(simdata$rna)  # muRna + muRna^2/thetaRna
#' 
#' mean(simdata$hto[1, class[1, ]])  # mu[1]
#' mean(simdata$hto[1, !class[1, ]]) # muAmbient[1]
#' var(simdata$hto[1, class[1, ]])   # > mu[1] + mu[1]^2/theta[1]
#' 
#' 
#' cor(simdata$rna[class[1, ]], simdata$hto[1, class[1, ]])
#' 
#' @aliases dmmSimulateHto,matrix-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("dmmSimulateHto",
           function(class, mu=180, theta=15, muAmbient=30, thetaAmbient=10, muRna=3000, thetaRna=30)
             standardGeneric("dmmSimulateHto"),
           signature="class")



#' Summarize classification results from demuxmix
#' 
#' This method takes the demultiplexing results from an HTO experiment
#' returned by \code{\link{demuxmix}} and returns a \code{data.frame}
#' summarizing the classification results and classification certainty.
#' 
#' @param dmmResults A "results" \code{data.frame} generated by
#'   \code{\link{demuxmix}}.
#' @param p.acpt An optional value between 0 and 1 defining the acceptance
#'   probability that must be reached in order to assign a cell to a hashtag.
#'   If missing, the original acceptance probability passed to
#'   \code{\link{demuxmix}} when generating \code{dmmResults} is used.
#'   
#' @details Results are summarized for the individual hashtags, for all
#'   singletons combined, for all multiplets, and for the negative class.
#'   Relative frequencies are calculated after excluding the "uncertain"
#'   class. The estimated number of false positive cells and the
#'   estimated FDR are based on several assumptions, one of which is the
#'   independence of the HTO counts from the different hashtags. This
#'   assumption unlikely for real data where HTO counts of different hashtags
#'   are obtained from the same cells. Error rates for singletons can be 
#'   expected to be conservative whereas error rates for multiplets and
#'   negative cells are likely underestimated.
#' 
#' @return A \code{data.frame} with one row per class showing the number
#'   of cells in the class (NumObs), the relative frequency of the class
#'   (RelFreq), the median probability that a cell assigned to the class
#'   is actually a member of the class (MedProb), the estimated number of
#'   cells falsely assigned to the class (ExpFPs), and the corresponding
#'   false discovery rate (FDR).
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples 
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' dmmSummary(dmm$results)
#' dmmSummary(dmm$results, p.acpt=0.5)
#' 
#' @aliases dmmSummary,data.frame-method
#' 
#' @importFrom methods setGeneric
###' @export
setGeneric("summary")
# setGeneric("dmmSummary",
#            function(dmmResults, p.acpt)
#              standardGeneric("dmmSummary"),
#            signature="dmmResults")



#' Plotting a histogram with mixture densities
#' 
#' This methods plots the mixture density and the components' densities
#' on top of a histogram of the HTO counts used to fit the mixture model.
#' The mixture model must be generated by \code{\link{demuxmix}}.
#' 
#' @param model A mixture model or a list of mixture models returned by
#'   \code{\link{demuxmix}}.
#' @param quantile Quantile of the mixture distribution which is used as
#'   right limit of the x axis of the plot.
#' @param binwidth Width of the bins of the histogram.
#' 
#' @details A histogram overlaid with the densities from the mixture model
#'   is often helpful to assess the mixture model fit. How such a plot
#'   should generated for regression mixture models is not obvious. This
#'   method calculates the weighted mean number of detected genes in the
#'   cells for each component, which are then used to calculate the
#'   expected number of HTO counts for the negative and positive component.
#'   The HTO counts shown in the histogram are adjusted to account for
#'   different numbers of detected genes by replacing the original HTO
#'   counts with the expected counts given the mean number of detected
#'   genes plus the residual counts from the regression model.
#'   
#'   Sometimes it is useful to zoom into the plot to obtain a better view
#'   of the fit. To restrict the plot to a certain range on the x or y axis,
#'   the method \code{\link[ggplot2]{coord_cartesian}} from the \code{ggplot2}
#'   package should be used (see examples).
#'   
#' @return An object of class \code{ggplot} is returned. Or, if multiple
#'   mixture models are passed to this function, an object of class
#'   \code{gtable} is returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, p.acpt=0.9)
#' \donttest{plotDmmHistogram(dmm$model)}
#' 
#' dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' \donttest{plotDmmHistogram(dmmreg$model)}
#' 
#' p <- plotDmmHistogram(dmmreg$model[1])
#' \donttest{p + ggplot2::coord_cartesian(xlim=c(25, 100), ylim=c(0, 0.01))}
#'
#' @aliases plotDmmHistogram,list-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmHistogram",
           function(object, hto, quantile=0.95, binwidth=5)
             standardGeneric("plotDmmHistogram"),
           signature=c("object", "hto"))



#' Plotting RNA features versus HTO counts
#' 
#' This methods plots the numbers of detected RNA features versus the numbers
#' of HTO counts and indicates the posterior probability obtained from the
#' mixture model by different colors. The mixture model passed to this
#' function must be a regression mixture model generated by
#' \code{\link{demuxmix}}.
#' 
#' @param model A regression mixture model or a list of regression mixture
#'   models returned by \code{\link{demuxmix}}.
#' @param log Logical value indicating whether both HTO counts and number
#'   of RNA features should be log transformed.
#' @param pointsize Numeric value specifying the size of the points.
#' 
#' @details The scatterplot produced by this method is helpful to assess
#'   the relation between the number of RNA features detected and the number
#'   of HTO counts obtained for a cell. A positive association is usually
#'   visible for the positive cells with a large posterior probability
#'   indicated by red color.
#'   
#' @return An object of class \code{ggplot} is returned. Or, if multiple
#'   mixture models are passed to this function, an object of class
#'   \code{gtable} is returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' \donttest{plotDmmScatter(dmmreg$model)}
#' \donttest{plotDmmScatter(dmmreg$model[1])}
#'
#' @aliases plotDmmScatter,list-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmScatter",
           function(object, hto, log=TRUE, pointsize=1.2)
             standardGeneric("plotDmmScatter"),
           signature=c("object", "hto"))



#' Plotting a histogram of posterior probabilities
#' 
#' This methods plots a histogram of posterior probabilities obtained
#' from the given mixture model. The probabilities indicate whether the
#' cell originates from the positive component of the mixture model (i.e.
#' tagged with the hashtag) or from the negative component. The mixture
#' model passed to this function must be generated by
#' \code{\link{demuxmix}}.
#' 
#' @param model A mixture model or a list of mixture models returned
#'   \code{\link{demuxmix}}.
#' @param log Logical value indicating whether a logarithmic scale
#'   should be used for the frequencies (y axis).
#' @param bins The number of bins of the histogram.
#' 
#' @details The histogram visualizes how well the positive cells can be
#'   separated from the cells negative for the hashtag. Ideally, the
#'   histogram shows many cells with a posterior probability very close
#'   to 0 and many cells close to 1, but no or very few cells with
#'   probabilities somewhere in between. The histogram can be useful to
#'   guide the choice of the acceptance probability \code{p.acpt} passed
#'   to \code{\link{demuxmix}}.
#'   
#' @return An object of class \code{ggplot} is returned. Or, if multiple
#'   mixture models are passed to this function, an object of class
#'   \code{gtable} is returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class=rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                       c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, p.acpt=0.9)
#' \donttest{plotDmmPosteriorP(dmm$model)}
#' 
#' dmmreg <- demuxmix(simdata$hto, rna=simdata$rna, p.acpt=0.9)
#' \donttest{plotDmmPosteriorP(dmmreg$model)}
#' \donttest{plotDmmPosteriorP(dmmreg$model[1])}
#'
#' @aliases plotDmmPosteriorP,list-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmPosteriorP",
           function(object, hto, bins=50)
             standardGeneric("plotDmmPosteriorP"),
           signature=c("object", "hto"))

#' @export
setGeneric("p.acpt",
           function(model)
             standardGeneric("p.acpt"),
           signature="model")

#' @export
setGeneric("p.acpt<-",
           function(model, value)
             standardGeneric("p.acpt<-"),
           signature=c("model", "value"))

#' @export
setGeneric("dmmClassify",
           function(model)
             standardGeneric("dmmClassify"),
           signature="model")

#' @export
setGeneric("getPosteriorProbability",
           function(model)
             standardGeneric("getPosteriorProbability"),
           signature="model")



# Internal methods for NaiveMixModel and RegMixModel
setGeneric("getHto",
           function(model, standardize=FALSE)
             standardGeneric("getHto"),
           signature="model")

setGeneric("getMu1",
           function(model, standardize=FALSE)
             standardGeneric("getMu1"),
           signature="model")

setGeneric("getMu2",
           function(model, standardize=FALSE)
             standardGeneric("getMu2"),
           signature="model")

setGeneric("getTheta1",
           function(model)
             standardGeneric("getTheta1"),
           signature="model")

setGeneric("getTheta2",
           function(model)
             standardGeneric("getTheta2"),
           signature="model")

setGeneric("getPi",
           function(model)
             standardGeneric("getPi"),
           signature="model")