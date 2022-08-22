test_that("Model fitting and classification with demuxmix", {
    set.seed(45274)
    class <- rbind(
        c(rep(TRUE, 220), rep(FALSE, 200)),
        c(rep(FALSE, 200), rep(TRUE, 220))
    )

    mu <- c(150, 300)
    theta <- c(15, 20)
    muAmbient <- c(30, 30)
    thetaAmbient <- c(10, 10)
    simData <- list()

    # Simulate data with pos. association between number of genes and HTO counts
    simData$simDep <- dmmSimulateHto(
        class = class,
        mu = mu,
        theta = theta,
        muAmbient = muAmbient,
        thetaAmbient = thetaAmbient,
        muRna = 3000,
        thetaRna = 30
    )

    # Simulate data without association between number of genes and HTO counts
    hto <- rbind(
        c(
            rnbinom(n = 220, mu = mu[1], size = theta[1]),
            rnbinom(n = 200, mu = muAmbient[1], size = thetaAmbient[1])
        ),
        c(
            rnbinom(n = 200, mu = muAmbient[2], size = thetaAmbient[2]),
            rnbinom(n = 220, mu = mu[2], size = theta[2])
        )
    )
    rownames(hto) <- c("HTO_1", "HTO_2")
    simData$simIndep <- list(
        hto = hto,
        rna = simData$simDep$rna,
        groundTruth = simData$simDep$groundTruth
    )

    # Fit models and test expectation values and classifications
    for (model in c("naive", "regpos", "reg")) {
        for (data in names(simData)) {
            dmm <- demuxmix(simData[[data]]$hto,
                rna = simData[[data]]$rna,
                model = model
            )
            classLabels <- dmmClassify(dmm)

            # mean of negative component of 1st HTO approximately muAmbient[1]
            expect_equal(demuxmix:::getMu1(dmm@models[[1]], standardize = TRUE),
                muAmbient[1],
                tolerance = 0.2
            )
            # mean of positive component of 1st HTO approximately mu[1]
            expect_equal(demuxmix:::getMu2(dmm@models[[1]], standardize = TRUE),
                mu[1],
                tolerance = 0.1
            )
            # mean of negative component of 2nd HTO approximately muAmbient[2]
            expect_equal(demuxmix:::getMu1(dmm@models[[2]], standardize = TRUE),
                muAmbient[2],
                tolerance = 0.2
            )
            # mean of positive component of 2nd HTO approximately mu[2]
            expect_equal(demuxmix:::getMu2(dmm@models[[2]], standardize = TRUE),
                mu[2],
                tolerance = 0.1
            )
            # not more than 10% classified as "uncertain" in these datasets
            expect_lt(
                sum(classLabels$Type == "uncertain") / nrow(classLabels),
                0.1
            )
            # more than 95% of singlets and doublets classified correctly
            valid <- classLabels$Type != "uncertain"
            expect_gt(sum(classLabels$HTO[valid] ==
                simData[[data]]$groundTruth[valid]) / sum(valid), 0.95)
        }
    }
})