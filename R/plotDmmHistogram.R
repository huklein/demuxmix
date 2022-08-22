#' @importFrom stats dnbinom qnbinom density
#' @importFrom ggplot2 ggplot geom_histogram stat_function xlab ylab coord_cartesian aes stat
.plotDmmHistogram <- function(model, quantile = 0.95, binwidth = 5) {
    hto <- getHto(model, standardize = TRUE)
    mu1 <- getMu1(model, standardize = TRUE)
    mu2 <- getMu2(model, standardize = TRUE)
    theta1 <- getTheta1(model)
    theta2 <- getTheta2(model)
    pi <- getPi(model)

    if (is(model, "NaiveMixModel")) {
        xlab <- paste(model@htoId, "HTO counts")
    } else {
        xlab <- paste(model@htoId, "HTO counts (adjusted for RNA)")
    }

    comp1 <- function(x) {
        return(dnbinom(round(x), mu = mu1, size = theta1) * pi[1])
    }
    comp2 <- function(x) {
        return(dnbinom(round(x), mu = mu2, size = theta2) * pi[2])
    }
    mixture <- function(x) {
        return(comp1(x) + comp2(x))
    }

    xmax <- max(
        qnbinom(quantile, mu = mu1, size = theta1),
        qnbinom(quantile, mu = mu2, size = theta2)
    )
    df <- data.frame(hto = hto)

    p <- ggplot(df, aes(x = hto)) +
        geom_histogram(aes(y = stat(density)), alpha = 0.4, binwidth = binwidth) +
        stat_function(fun = mixture, lwd = 1, n = xmax, col = "black") +
        stat_function(fun = comp1, lwd = 1, n = xmax, col = "dodgerblue") +
        stat_function(fun = comp2, lwd = 1, n = xmax, col = "firebrick2") +
        xlab(xlab) +
        ylab("Density") +
        coord_cartesian(xlim = c(0, xmax))

    return(p)
}


#' @importFrom gridExtra grid.arrange
setMethod("plotDmmHistogram",
    signature = c(object = "Demuxmix", hto = "missing"),
    function(object, hto, quantile = 0.95, binwidth = 5) {
        plots <- lapply(object@models, .plotDmmHistogram,
            quantile = quantile, binwidth = binwidth
        )
        if (length(plots) == 1) {
            return(plots[[1]])
        } else {
            return(do.call(grid.arrange, plots))
        }
    }
)

setMethod("plotDmmHistogram",
    signature = c(object = "Demuxmix", hto = "ANY"),
    function(object, hto, quantile = 0.95, binwidth = 5) {
        if (is.numeric(hto) & any(hto > length(object@models))) {
            stop("Invalid HTO identifier.")
        }
        if (is.character(hto) & any(!is.element(hto, names(object@models)))) {
            stop("Invalid HTO identifier.")
        }
        plots <- lapply(object@models[hto], .plotDmmHistogram,
            quantile = quantile, binwidth = binwidth
        )
        if (length(plots) == 1) {
            return(plots[[1]])
        } else {
            return(do.call(grid.arrange, plots))
        }
    }
)