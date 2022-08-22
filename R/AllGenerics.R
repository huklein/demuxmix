#' Demultiplexing using mixture models
#'
#' This method uses mixture models as probabilistic framework to assign
#' droplets to hashtags and to identify multiplets based on counts obtained from
#' a hashtag oligonucleotide (HTO) library. If the numbers of detected genes
#' from the corresponding RNA library are passed as second argument,
#' regression mixture models may be used, which often improves the
#' classification accuracy by leveraging the relationship between HTO and RNA
#' read counts.
#' 
#' @param hto A matrix of HTO counts where each row corresponds to a hashtag
#'   and each column to a droplet. The matrix must have unique row names.
#' @param rna An optional numeric vector with the number of genes detected in
#'   the RNA library for each droplet. Same length as columns in \code{hto}.
#'   If missing, parameter \code{model} must be set to "naive".
#' @param pAcpt Acceptance probability that must be reached in order to
#'   assign a droplet to a hashtag. Droplets with lower probabilities are
#'   classified as "uncertain". This parameter can be changed after running
#'   demuxmix by applying \code{\link{pAcpt<-}} to the returned object.
#' @param model A character specifying the type of mixture model to be used.
#'   Either "naive", "regpos", "reg" or "auto". The last three options
#'   require parameter \code{rna} to be specified. "auto" selects the best
#'   model based on the classification error probability summed over all
#'   droplets.
#' @param alpha Threshold defining the left tail of the mixture
#'   distribution where droplets should not be classified as "positive".
#'   Threshold must be between 0 and 1. See details.
#' @param beta Threshold for defining the right tail of the mixture
#'   distribution where droplets should not be classified as "negative".
#'   Threshold must be between 0 and 1. See details.
#' @param correctTails If \code{TRUE}, droplets meeting the threshold defined by
#'   \code{alpha} (\code{beta}) are classified as "negative" ("positive") even
#'   if the mixture model suggests a different classification. See details.
#' @param tol Convergence criterion for the EM algorithm used to fit the mixture
#'   models. The algorithm stops when the relative increase of the log
#'   likelihood is less than or equal to \code{tol}.
#' @param maxIter Maximum number of iterations for the EM algorithm and
#'   for the alternating iteration process fitting the NB regression models
#'   within each EM iteration.
#' @param k.hto Factor to define outliers in the HTO counts. Among droplets
#'   positive for the hashtag based on initial clustering, HTO counts
#'   larger than the 0.75 quantile + \code{k.hto} * IQR are considered
#'   outliers. See details.
#' @param k.rna Factor to define outliers in the numbers of detected genes.
#'   Numbers of detected genes larger than the 0.75 quantile +
#'   \code{k.rna} * IQR are considered outliers. See details.
#' 
#' @details The single cell dataset should undergo basic filtering to
#'   remove low quality or empty droplets before calling this function,
#'   but the HTO counts should not be transformed or pre-processed
#'   otherwise. The number of detected genes passed via the optional argument
#'   \code{rna} is typically defined as the number of genes in the RNA library
#'   with at least one read.
#'   
#'   The method fits a two-component negative binomial mixture model for each
#'   hashtag. The type of mixture model used can be specified by \code{model}.
#'   "naive" fits a standard mixture model. "reg" fits a regression mixture
#'   model using the given number of detected genes (\code{rna}) as covariate
#'   in the regression model. "regpos" uses a regression model only for the
#'   positive but not for the negative component. If model is set to "auto",
#'   all three models are fitted and the model with the lowest posterior
#'   classification error probability summed over all droplets is selected.
#'   Details are stored in the slot \code{modelSelection} of the returned
#'   object. In most real HTO datasets, regression mixture models outperform
#'   the naive mixture model.
#'   
#'   The \code{demuxmix} method consists of 3 steps, which can be tuned
#'   by the respective parameters. The default settings work well for a wide
#'   range of datasets and usually do not need to be adapted unless any issues
#'   arise during model fitting and quality control. An exception is the
#'   acceptance probability \code{pAcpt}, which may be set to smaller or
#'   larger value depending on the desired trade-off between number of
#'   unclassified/discarded droplets and expected error rate. Steps 1 and 2
#'   are executed for each HTO separately; step 3 classifies the droplets based
#'   on the results from all HTOs. Therefore, parameters affecting steps 1 and 2
#'   (incl. \code{model}) can be specified for each HTO using a vector with
#'   one element per HTO. Shorter vectors will be extended.
#'   
#'   \enumerate{
#'     \item Preprocessing (\code{k.hto}, \code{k.rna}). Droplets are clustered
#'     into a negative and a positive group based on the HTO counts using
#'     k-means. Droplets in the positive group with HTO counts larger than the
#'     0.75 quantile + \code{k.hto} times the IQR of the HTO counts in the
#'     positive group are marked as outliers. Outliers are still classified but
#'     will not be used to fit the mixture model for this HTO in step 2. If
#'     the parameter \code{rna} is given and the \code{model} is "reg" or
#'     "regpos", all droplets (both groups) with number of detected genes larger
#'     than the 0.75 quantile + k.rna times the IQR are marked as outliers, too,
#'     since these cells could affect the fitting of the regression model
#'     negatively. If more than 5\% of the cells are marked as outliers,
#'     a warning message is printed and larger values for \code{k.hto}
#'     and \code{k.rna} might be preferable. If the model fit seems to be
#'     affected by a few large values (very high variance of the positive
#'     component), smaller values should be chosen.
#'     
#'     \item Model fitting (\code{model}, \code{alpha}, \code{beta},
#'     \code{correctTails}, \code{tol}, \code{maxIter}). An EM algorithm is
#'     used to fit the mixture model to the HTO counts which were not marked as
#'     outliers in step 1. \code{maxIter} defines the maximum number of
#'     iterations of the EM algorithm, and, if \code{model} is "reg", "regpos"
#'     or "auto", it also defines the maximum number of iterations to fit the
#'     negative binomial regression models within each EM iteration. \code{tol}
#'     defines the convergence criterion for the EM algorithm. The algorithm
#'     stops if \eqn{\Delta LL/LL \le} \code{tol}. After the mixture model has
#'     been fitted, the posterior probability that the i-th droplet is positive
#'     for the hashtag \eqn{P(C_i = pos)} is calculated. Depending on the given
#'     data, these probabilities can be inaccurate at the far tails of the
#'     mixture distribution. Specifically, a positive component with large
#'     variance can have a larger value close to zero than the negative
#'     component, if the negative component is narrow and shifted to the right
#'     due to background HTO reads. If \code{correctTails} is \code{TRUE}, the
#'     following two rules are applied to avoid false classifications at the
#'     far tails. First, if the i-th droplet is classified as positive based on
#'     the posterior probability, but the probability to detected more than the
#'     observed \eqn{y_i} HTO counts in a negative droplet is
#'     \eqn{P(Y \ge y_i | neg)} > \code{alpha}, then \eqn{P(C_i = pos)} is set
#'     to 0 (left tail). Second, if the i-th droplet is classified as negative,
#'     but \eqn{P(Y \le y_i | pos)} > \code{beta}, \eqn{P(C_i = pos)} is set to
#'     1 (right tail). For most datasets, these rules will not apply and it is
#'     recommended not to change these values. If \code{correctTails} is
#'     \code{FALSE}, posterior probabilities will not be altered, but potential
#'     problems at the tails will still be logged in the slot
#'     \code{tailException} of the returned object.
#'     
#'     \item Classification (\code{pAcpt}). The posterior probabilities
#'     obtained from the models fitted to each HTO separately are 
#'     used to calculate the most likely class for each cell. The following
#'     classes are considered: one class for each HTO (singlets), one
#'     class for each possible multiplet, and a negative class representing
#'     droplets negative for all HTOs (i.e. empty droplets or droplets
#'     containing only cell debris). Each droplet is assigned to the most
#'     likely class unless the probability is smaller than \code{pAcpt},
#'     in which case the droplet is assigned to the class "uncertain".
#'     Classification results can be accessed by running
#'     \code{\link{dmmClassify}} on an object returned by \code{demuxmix}. The
#'     acceptance probability can be changed after running \code{demuxmix} using
#'     \code{\link{pAcpt<-}}.
#'   }
#'   
#' @return \code{demuxmix} returns an object of class \code{\link{Demuxmix}}.
#'   Classification results can be extracted with \code{\link{dmmClassify}}.
#'   Various plot methods (see below) are available to assess the model fit.
#'  
#' @seealso \code{\link{dmmClassify}} to extract the classification results
#'   and \code{\link{summary}} to summarize the results.
#'   \code{\link{plotDmmHistogram}}, \code{\link{plotDmmScatter}},
#'   \code{\link{plotDmmPosteriorP}}, and \code{\link{dmmOverlap}} to assess
#'   the model fit.
#'   
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, model = "naive")
#' dmm
#' table(dmmClassify(dmm)$HTO, simdata$groundTruth)
#' 
#' dmmreg <- demuxmix(simdata$hto, rna = simdata$rna)
#' dmm
#' table(dmmClassify(dmmreg)$HTO, simdata$groundTruth)
#' summary(dmmreg)
#' 
#' pAcpt(dmmreg) <- 0.5
#' summary(dmmreg)
#' 
#' dmmOverlap(dmmreg)
#' \donttest{
#' plotDmmHistogram(dmmreg)
#' plotDmmScatter(dmmreg, hto="HTO_1")}
#' 
#' @aliases demuxmix,matrix,missing-method demuxmix,matrix,numeric-method
#'   demuxmix,Matrix,missing-method demuxmix,Matrix,numeric-method
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("demuxmix",
    function(hto, rna, pAcpt = 0.9^nrow(hto), model = "auto",
             alpha = 0.9, beta = 0.9, correctTails = TRUE,
             tol = 10^-5, maxIter = 100, k.hto = 1.5, k.rna = 1.5) {
        standardGeneric("demuxmix")
    },
    signature = c("hto", "rna")
)



#' Return classification results from a Demuxmix object
#' 
#' This method uses the posterior probabilities from the given demuxmix model
#' to assign each droplet to the most likely class, either a single HTO,
#' a combination of HTOs (multiplet) or the negative class (non-labeled
#' cells, empty droplets, cell debris). If the assignment cannot be made with
#' certainty above a defined threshold, the droplet is labeled as "uncertain".
#' 
#' @param object An object of class \code{\link{Demuxmix}}.
#' 
#' @details A droplet is labeled as "uncertain" if the posterior probability of
#'   the most likely class is smaller than the threshold \code{pAcpt}, which
#'   is stored in the given \code{\link{Demuxmix}} object. The acceptance
#'   probability \code{pAcpt} can be inspected and set to a different value
#'   by applying the getter/setter method \code{\link{pAcpt}} to the
#'   \code{\link{Demuxmix}} object before calling this method. The method
#'   \code{\link{summary}} is useful to inspect classification results and to
#'   estimate error rates for different values of \code{pAcpt}.
#' 
#' @return A \code{data.frame} with 3 columns and one row for each droplet
#'   in the dataset. The first column gives the class (HTO) the droplet has been
#'   assigned to. The second column contains the posterior probability. And the
#'   third column specifies the type of the assigned class, i.e., "singlet",
#'   "multiplet", "negative" or "uncertain".
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, rna = simdata$rna)
#' head(dmmClassify(dmm))
#' table(dmmClassify(dmm)$HTO, simdata$groundTruth)
#' 
#' pAcpt(dmm) <- 0.5
#' sum(dmmClassify(dmm)$HTO == "uncertain")
#' pAcpt(dmm) <- 0.9999
#' sum(dmmClassify(dmm)$HTO == "uncertain")
#' 
#' @aliases dmmClassify,Demuxmix-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("dmmClassify",
    function(object) {
        standardGeneric("dmmClassify")
    },
    signature = "object"
)



#' Calculate the intersection of two components of a mixture model
#' 
#' \code{dmmOverlap} sums over the probability mass intersected by the two
#' components of the given mixture model. The sum should be
#' close to 0 if the HTO labeling experiment was successful.
#'
#' @param object An object of class \code{\link{Demuxmix}}.
#' @param hto Optional vector specifying a subset of HTOs in \code{object} which
#'   should be used by this function.
#' @param tol The maximum acceptable error when calculating the
#'   area.
#'   
#' @details The probability mass shared between the negative and positive
#'   component is an informative quality metric for the labeling efficiency
#'   of the HTO. Values under 0.03 can be considered as good, values larger
#'   than 0.1 are problematic.
#'   
#'   The probability mass functions of the negative and positive component are
#'   not scaled by the estimated proportions of negative and positive droplets.
#'   Therefore, the result does not depend on the proportion of cells stained
#'   with the HTO and the returned value lies between 0 and 1. 
#'   
#'   The definition of the shared probability mass is not obvious for a
#'   regression mixture model since the distributions' means depend on the
#'   covariate, i.e., the number of detected genes in the RNA library. If a
#'   regression mixture model is given, this method calculates for each of the
#'   two components the weighted mean number of detected genes and uses these
#'   numbers to calculate the expectation value for the negative
#'   and positive component respectively.
#' 
#' @return A numeric vector with the shared probability mass for each HTO in
#'   the given \code{object}.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples 
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, model = "naive")
#' dmmOverlap(dmm)
#' 
#' dmmreg <- demuxmix(simdata$hto, rna = simdata$rna)
#' dmmOverlap(dmmreg)
#' dmmOverlap(dmmreg, hto = "HTO_1")
#' dmmOverlap(dmmreg, hto = 2)
#' 
#' @aliases dmmOverlap,Demuxmix,missing-method dmmOverlap,Demuxmix,ANY-method
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("dmmOverlap",
    function(object, hto, tol = 0.001) {
        standardGeneric("dmmOverlap")
    },
    signature = c("object", "hto")
)



#' Simulate HTO sequencing data
#' 
#' This method simulates HTO count data and corresponding numbers of detected
#' RNA features using the negative binomial distribution. The purpose of this
#' method is to provide simple example datasets for testing and documentation.
#'
#' @param class A \code{matrix} of type logical defining the number of
#'   HTOs, the number of droplets, and the droplets' class memberships, i.e.,
#'   which droplets contain cells that have been tagged with a certain HTO.
#'   Each row corresponds to one HTO and each column to a droplet. Negative
#'   droplets (all entries in the column are \code{FALSE}) and multiplets
#'   (more than one entry are \code{TRUE}) are allowed. If the matrix has row
#'   names, the names must be unique and are used as HTO names.
#' @param mu Vector of expectation values of the HTO counts if a
#'   droplet is positive for the HTO. Values are recycled if \code{mu} is
#'   shorter than number of HTOs defined by \code{class}.
#' @param theta Vector of dispersion parameters of the HTO counts if a
#'   droplet is positive for the HTO. Values are recycled if \code{theta} is
#'   shorter than number of HTOs defined by \code{class}.
#' @param muAmbient Vector of expectation values of the HTO counts if a
#'   droplet is negative for the HTO. Values are recycled if \code{mu} is
#'   shorter than number of HTOs defined by \code{class}.
#' @param thetaAmbient Vector of dispersion parameters of the HTO counts if a
#'   droplet is negative for the HTO. Values are recycled if \code{theta} is
#'   shorter than number of HTOs defined by \code{class}.
#' @param muRna Single expectation value for the number of detected RNA
#'   features.
#' @param thetaRna Single dispersion parameter for the number of detected RNA
#'   features.
#'
#' @details A vector \eqn{r} of detected RNA features (same length as columns
#'   in \code{class}) is simulated using \code{\link[stats]{rnbinom}} with
#'   \code{muRna} and \code{thetaRna} as parameters. HTO counts of positive
#'   droplets are then simulated using \code{\link[stats]{rnbinom}} with 
#'   \eqn{r} \code{mu}/\code{muRna} as expectation value and \code{theta}
#'   as dispersion. If a droplet is negative for the HTO,
#'   \eqn{r} \code{muAmbient}/\code{muRna} and \code{thetaAmbient}
#'   are used respectively.
#' 
#' @return A list with three elements: \code{hto} is a matrix of the same
#'   dimension as the given \code{class} matrix and contains the simulated HTO
#'   counts. \code{rna} is a vector of simulated detected number of genes (same
#'   length as \code{hto} has columns). \code{groundTruth} is a character
#'   vector encoding the class labels given by \code{class} as character
#'   strings for convenience.
#'
#' @seealso \code{\link{demuxmix}} 
#' 
#' @examples
#' set.seed(2642)
#' class <- rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                c(rep(FALSE, 200), rep(TRUE, 220)))
#' simdata <- dmmSimulateHto(class = class, mu = c(150, 300), theta = c(15, 20),
#'                           muAmbient = c(30, 30), thetaAmbient = c(10, 10),
#'                           muRna = 3000, thetaRna = 30)
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
    function(class, mu = 180, theta = 15,
             muAmbient = 30, thetaAmbient = 10,
             muRna = 3000, thetaRna = 30) {
        standardGeneric("dmmSimulateHto")
    },
    signature = "class"
)



#' Summarize classification results of a Demuxmix model
#' 
#' This method takes the demultiplexing results from an HTO experiment
#' returned by \code{\link{demuxmix}} and returns a \code{data.frame}
#' summarizing the classification results and expected error rates.
#' 
#' @param object An object of class \code{\link{Demuxmix}}.
#' @param ... Additional parameters (ignored).
#' 
#' @details Results are summarized for the individual HTOs, for all
#'   singlets combined, for all multiplets combined, and for the negative
#'   class. Relative frequencies are calculated after excluding the
#'   "uncertain" class. The estimated number of false positive droplets and the
#'   estimated FDR are based on several assumptions, one of which is the
#'   independence of the HTO counts from different hashtags. This
#'   assumption is unlikely for real data where all HTO counts
#'   are obtained from the same droplet. Usually, the positive
#'   correlation among HTOs causes an overestimation of multiplets and
#'   negative/empty droplets. Error rates are more accurate when regression
#'   mixture models are used since the number of detected genes explains some
#'   of the positive correlation between HTOs.
#' 
#' @return A \code{data.frame} with one row per class showing the number
#'   of droplets in the class (NumObs), the relative frequency of the class
#'   (RelFreq), the median probability with which a droplet was assigned to the
#'   class (MedProb), the estimated number of droplets falsely assigned to the
#'   class (ExpFPs), and the corresponding estimated false discovery rate (FDR).
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples 
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, rna = simdata$rna)
#' summary(dmm)
#' pAcpt(dmm) <- 0.05
#' summary(dmm)
#' 
#' @aliases summary,data.frame-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("summary")



#' Plotting a histogram with mixture probability mass function
#' 
#' This methods plots the mixture probability mass function with the negative
#' and positive component on top of a histogram of the HTO counts used to fit
#' the mixture model. The mixture model must be generated by
#' \code{\link{demuxmix}}.
#' 
#' @param object An object of class \code{\link{Demuxmix}}.
#' @param hto Optional vector specifying a subset of HTOs in \code{object} which
#'   should be used by this function.
#' @param quantile Quantile of the mixture distribution which is used as
#'   right limit of the plot's x axis.
#' @param binwidth Width of the bins of the histogram.
#' 
#' @details A histogram overlaid with the pmf is a standard tool to 
#'   assess the fit of a the mixture model and trivial for a naive mixture
#'   model. However, if a regression mixture model is given, the expectation
#'   values of the components are different for each droplet depending on the
#'   covariates (here the number of genes detected in the droplet). This method
#'   calculates the weighted mean number of detected genes in droplets in the
#'   positive and negative component, and then uses these numbers to
#'   calculate expectation values for an average droplet of the positive
#'   and negative component. The HTO counts shown in the histogram
#'   are adjusted to account for different numbers of detected genes by
#'   replacing the original HTO counts with the expected counts given the mean
#'   number of detected genes plus the residuals from the regression model. In
#'   other words, the effect of the number of detected genes was regressed out
#'   before plotting the HTO counts in the histogram.
#'   
#'   It may be useful to zoom into the plot to obtain a better view
#'   of the fit. To restrict the plot to a certain range on the x or y axis,
#'   the method \code{\link[ggplot2]{coord_cartesian}} from the \code{ggplot2}
#'   package should be used (see examples).
#'   
#' @return An object of class \code{ggplot} is returned, if only one HTO is
#'   plotted. If several HTOs are plotted simultaneously, a grid of plots is
#'   returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, simdata$rna)
#' \donttest{plotDmmHistogram(dmm)}
#' p <- plotDmmHistogram(dmm, hto = 1)
#' \donttest{p + ggplot2::coord_cartesian(xlim = c(25, 100), ylim = c(0, 0.01))}
#'
#' @aliases plotDmmHistogram,Demuxmix,missing-method
#'   plotDmmHistogram,Demuxmix,ANY-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmHistogram",
    function(object, hto, quantile = 0.95, binwidth = 5) {
        standardGeneric("plotDmmHistogram")
    },
    signature = c("object", "hto")
)



#' Plotting RNA features versus HTO counts
#' 
#' This methods plots the number of genes detected in a droplet versus the
#' number of sequenced HTOs. The posterior probability that the droplet is
#' positive for the HTO is indicated by a color gradient. Optionally, the
#' decision boundary with posterior probability 0.5 can be plotted. The mixture
#' model passed to this function must be a regression mixture model
#' generated by \code{\link{demuxmix}}.
#' 
#' @param object An object of class \code{\link{Demuxmix}}.
#' @param hto Optional vector specifying a subset of HTOs in \code{object} which
#'   should be used by this function.
#' @param log Logical value indicating whether both HTO counts and number
#'   of detected genes should be log transformed.
#' @param pointsize Numeric value specifying the size of the points.
#' @param plotDecBoundary Logical value indicating whether the decision
#'   boundary should be added to the plot.
#' @param tol Numeric value between 0 and 1 specifying the error tolerance of
#'   the decision boundary, i.e., a point on the plotted line has a posterior
#'   probability within 0.5 +- \code{tol}. Only used of \code{plotDecBoundary}
#'   is \code{true}.
#' 
#' @details The scatterplot produced by this method is helpful to assess
#'   the relation between the number of detected genes and the number
#'   of HTO counts obtained for a droplet. A positive association is usually
#'   visible for the positive cells (i.e., droplets with cells treated with the
#'   oligo-labeled antibodies). The association is often weak/absent in the
#'   droplets negative for the HTO.
#'   This method can only be applied to regression mixture models and not
#'   to naive mixture models. To see whether a \code{\link{Demuxmix}} object
#'   contains regression mixture models, type \code{show(object)} to display
#'   the type of model used for each HTO.
#'   
#' @return An object of class \code{ggplot} is returned, if only one HTO is
#'   plotted. If several HTOs are plotted simultaneously, a grid of plots is
#'   returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmmreg <- demuxmix(simdata$hto, rna = simdata$rna, model = "reg")
#' \donttest{plotDmmScatter(dmmreg)}
#' \donttest{plotDmmScatter(dmmreg, hto = 1, log = FALSE)}
#'
#' @aliases plotDmmScatter,Demuxmix,missing-method
#'   plotDmmScatter,Demuxmix,ANY-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmScatter",
    function(object, hto, log = TRUE, pointsize = 1.2,
             plotDecBoundary = TRUE, tol = 0.01) {
        standardGeneric("plotDmmScatter")
    },
    signature = c("object", "hto")
)



#' Plotting a histogram of posterior probabilities
#' 
#' This methods plots a histogram of posterior probabilities obtained
#' from the given mixture model. The posterior probabilities indicate whether
#' the droplet likely contains a cell labeled by the respective HTO.
#' The mixture model passed to this function must be generated by
#' \code{\link{demuxmix}}.
#'
#' @param object An object of class \code{\link{Demuxmix}}.
#' @param hto Optional vector specifying a subset of HTOs in \code{object} which
#'   should be used by this function.
#' @param bins The number of bins of the histogram.
#' 
#' @details The histogram visualizes how well the positive droplets can be
#'   separated from the negative droplets. Ideally, the histogram shows many
#'   droplets with a posterior probability very close to 0 and many droplets
#'   close to 1, but no or very few droplets with probabilities somewhere in
#'   between. The histogram can be useful for guiding the selection of the
#'   acceptance probability \code{pAcpt}.
#'   
#' @return An object of class \code{ggplot} is returned, if only one HTO is
#'   plotted. If several HTOs are plotted simultaneously, a grid of plots is
#'   returned.
#' 
#' @seealso \code{\link{demuxmix}}
#' 
#' @examples
#' set.seed(2642)
#' simdata <- dmmSimulateHto(class = rbind(c(rep(TRUE, 220), rep(FALSE, 200)),
#'                                         c(rep(FALSE, 200), rep(TRUE, 220))))
#' 
#' dmm <- demuxmix(simdata$hto, model = "naive")
#' \donttest{plotDmmPosteriorP(dmm)}
#' 
#' dmmreg <- demuxmix(simdata$hto, rna = simdata$rna, model = "auto")
#' \donttest{plotDmmPosteriorP(dmmreg)}
#' \donttest{plotDmmPosteriorP(dmmreg, hto = 1)}
#'
#' @aliases plotDmmPosteriorP,Demuxmix,missing-method
#'   plotDmmPosteriorP,Demuxmix,ANY-method
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("plotDmmPosteriorP",
    function(object, hto, bins = 50) {
        standardGeneric("plotDmmPosteriorP")
    },
    signature = c("object", "hto")
)



#' @export
setGeneric("pAcpt",
    function(object) {
        standardGeneric("pAcpt")
    },
    signature = "object"
)

#' @export
setGeneric("pAcpt<-",
    function(object, value) {
        standardGeneric("pAcpt<-")
    },
    signature = c("object", "value")
)



# Internal methods for NaiveMixModel and RegMixModel
setGeneric("getPosteriorProbability",
    function(model) {
        standardGeneric("getPosteriorProbability")
    },
    signature = "model"
)


setGeneric("getHto",
    function(model, standardize = FALSE) {
        standardGeneric("getHto")
    },
    signature = "model"
)

setGeneric("getMu1",
    function(model, standardize = FALSE) {
        standardGeneric("getMu1")
    },
    signature = "model"
)

setGeneric("getMu2",
    function(model, standardize = FALSE) {
        standardGeneric("getMu2")
    },
    signature = "model"
)

setGeneric("getTheta1",
    function(model) {
        standardGeneric("getTheta1")
    },
    signature = "model"
)

setGeneric("getTheta2",
    function(model) {
        standardGeneric("getTheta2")
    },
    signature = "model"
)

setGeneric("getPi",
    function(model) {
        standardGeneric("getPi")
    },
    signature = "model"
)