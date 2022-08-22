test_that("dmmOverlap calculates shared probablity mass accurately", {
    mm <- demuxmix:::NaiveMixModel(
        mu1 = 100,
        mu2 = 150,
        theta1 = 20,
        theta2 = 5,
        pi = c(0.5, 0.5),
        hto = c(
            rnbinom(100, mu = 100, size = 20),
            rnbinom(100, mu = 150, size = 5)
        ),
        htoId = "myHTO",
        parameters = list(tol = 10^-5, maxIter = 100),
        log = data.frame(),
        converged = TRUE
    )

    dmm <- demuxmix:::Demuxmix(
        models = list(mm),
        outliers = matrix(NA),
        clusterInit = matrix(NA),
        posteriorProb = matrix(c(rep(0, 100), rep(1, 100)),
            nrow = 1,
            dimnames = list("myHTO", NULL)
        ),
        tailException = matrix(NA),
        modelSelection = data.frame(),
        parameters = list(pAcpt = 0.9)
    )

    # dnbinom(10000, mu = 150, size = 20) # -> 0
    # sum(pmin(
    #     dnbinom(0:10000, mu = 100, size = 20),
    #     dnbinom(0:10000, mu = 150, size = 5)
    # ))
    #    0.5277511

    # tolerance not violated
    expect_equal(dmmOverlap(dmm, tol = 0.8), 0.5277511,
        tolerance = 0.8 / 0.5277511
    )
    expect_equal(dmmOverlap(dmm, tol = 0.1), 0.5277511,
        tolerance = 0.1 / 0.5277511
    )
    expect_equal(dmmOverlap(dmm, tol = 0.001), 0.5277511,
        tolerance = 0.001 / 0.5277511
    )

    # mixture components exchangeable
    dmmSwap <- dmm
    dmmSwap@models[[1]]@mu1 <- dmm@models[[1]]@mu2
    dmmSwap@models[[1]]@mu2 <- dmm@models[[1]]@mu1
    dmmSwap@models[[1]]@theta1 <- dmm@models[[1]]@theta2
    dmmSwap@models[[1]]@theta2 <- dmm@models[[1]]@theta1
    expect_equal(dmmOverlap(dmm), dmmOverlap(dmmSwap))
})