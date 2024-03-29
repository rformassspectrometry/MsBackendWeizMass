test_that("MsBackendWeizMass class works", {
    obj <- new("MsBackendWeizMass", dbcon = dbc)
    expect_true(validObject(obj))
    obj <- new("MsBackendWeizMass")
    expect_true(validObject(obj))

    show(obj)
})

test_that("backendInitialize,MsBackendWeizMass works", {
    expect_error(backendInitialize(MsBackendWeizMass()), "required")

    be <- backendInitialize(MsBackendWeizMass(), dbcon = dbc)
    expect_true(length(be@spectraIds) > 0)
    expect_true(length(be@spectraVariables) > 0)
    expect_output(show(be), "MsBackendWeizMass")
})

test_that("peaksVariables,MsBackendWeizMass works", {
    be <- backendInitialize(MsBackendWeizMass(), dbcon = dbc)
    res <- peaksVariables(be)
    expect_equal(res, c("mz", "intensity", "relative_intensity",
                        "peak_annotation"))

    be <- MsBackendWeizMass()
    res <- peaksVariables(be)
    expect_equal(res, character())
})

test_that("length,MsBackendWeizMass works", {
    expect_equal(length(MsBackendWeizMass()), 0L)
    be <- backendInitialize(MsBackendWeizMass(), dbc)
    expect_equal(length(be), 2L)
})

test_that("dataStorage,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    expect_equal(dataStorage(be), character())

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    expect_equal(dataStorage(be), rep("<WeizMass>", length(be)))
})

test_that("peaksVariables,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    expect_equal(peaksVariables(be), character())

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    expect_equal(peaksVariables(be),
                 c("mz", "intensity", "relative_intensity", "peak_annotation"))
})

test_that("peaksData,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- peaksData(be)
    expect_true(is.list(res))

    be <- backendInitialize(be, dbc)
    res <- peaksData(be)
    expect_true(is.list(res))
    expect_true(length(res) == length(be))
    expect_true(is.matrix(res[[1]]))
    expect_true(is.matrix(res[[2]]))

    ## duplicated spectra.
    be2 <- be[c(2, 1, 1, 2)]
    res2 <- peaksData(be2)
    expect_equal(res2[[1]], res[[2]])
    expect_equal(res2[[2]], res[[1]])
    expect_equal(res2[[3]], res[[1]])
    expect_equal(res2[[4]], res[[2]])
})

test_that("spectraVariables,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- spectraVariables(be)
    expect_equal(res, names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- spectraVariables(be)
    expect_true(length(res) > length(names(Spectra:::.SPECTRA_DATA_COLUMNS)))

    be$a <- 1
    res_2 <- spectraVariables(be)
    expect_true(length(res_2) == (length(res) + 1))
})

test_that("spectraData,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)

    ##
    be <- backendInitialize(be, dbc)
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 2)
    expect_true(all(c("rtime", "mz", "intensity", "formula", "exactmass",
                      "smiles", "inchikey", "common_name", "iupac_name",
                      "peak_annotation") %in% colnames(res)))

    res <- spectraData(be, columns = c("common_name", "iupac_name"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 2)
    expect_equal(colnames(res), c("common_name", "iupac_name"))
    expect_true(is(res$common_name, "CharacterList"))
    expect_true(is(res$iupac_name, "CharacterList"))

    ## Add an additional spectra variable and require spectra data to be
    ## equal with that additional variable
    be_2 <- be
    be_2$other_column <- c("A", "B")
    res <- spectraData(be)
    res_2 <- spectraData(be_2)
    expect_equal(sort(c(colnames(res), "other_column")), sort(colnames(res_2)))

    ## Different order
    be_2 <- be[c(2, 1, 1)]
    res_2 <- spectraData(be_2)
    rownames(res) <- NULL
    rownames(res_2) <- NULL
    expect_equal(res_2[1, ], res[2, ])
    expect_equal(res_2[2, ], res[1, ])
    expect_equal(res_2[3, ], res[1, ])
})

test_that("$,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- be$rtime
    expect_true(is.numeric(res))
    expect_true(length(res) == 0)

    res <- be$mz
    expect_true(is(res, "NumericList"))

    expect_error(be$other, "not available")

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    expect_true(length(be$mz) == 2)
})

test_that("$<-,MsBackendWeizMass works", {
    be <- backendInitialize(MsBackendWeizMass(), dbc)
    ## Test adding new column
    be$new_col <- "a"
    expect_equal(be$new_col, rep("a", length(be)))

    ## Test replacing column
    be$new_col <- seq_len(length(be))
    expect_equal(be$new_col, seq_len(length(be)))

    be$rtime <- seq_len(length(be))
    expect_equal(be$rtime, seq_len(length(be)))

    be$authors <- "a"
    expect_equal(be$authors, rep("a", length(be)))

    ## Test errors with m/z etc.
    expect_error(be$rtime <- c(1, 2, 3), "length 1 or")
    expect_error(be$mz <- 1:3, "Replacing m/z and intensity")
})

test_that("acquisitionNum,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- acquisitionNum(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- acquisitionNum(be)
    expect_equal(res, rep(NA_integer_, length(be)))
})

test_that("centroided,centroided<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- centroided(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- centroided(be)
    expect_equal(res, rep(NA, length(be)))

    centroided(be) <- FALSE
    res <- centroided(be)
    expect_equal(res, rep(FALSE, length(be)))

    expect_error(centroided(be) <- 3, "logical")
})

test_that("collisionEnergy,collisionEnergy<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- collisionEnergy(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- collisionEnergy(be)
    expect_true(is.numeric(res))

    be$collisionEnergy <- 1.2
    res <- collisionEnergy(be)
    expect_equal(res, rep(1.2, length(be)))

    expect_error(collisionEnergy(be) <- "a", "numeric")
})

test_that("dataOrigin, dataOrigin<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- dataOrigin(be)
    expect_equal(res, character())

    be <- backendInitialize(be, dbc)
    res <- dataOrigin(be)
    expect_true(is.character(res))

    dataOrigin(be) <- "b"
    res <- dataOrigin(be)
    expect_equal(res, rep("b", length(be)))

    expect_error(dataOrigin(be) <- 1, "character")
})

test_that("selectSpectraVariables,MsBackendWeizMass works", {
    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- selectSpectraVariables(be, spectraVariables = spectraVariables(be))
    expect_equal(spectraVariables(be), spectraVariables(res))

    ## errors
    expect_error(selectSpectraVariables(be, c("rtime", "other")), "available")

    res <- selectSpectraVariables(be, c("rtime", "mz", "intensity"))
    expect_equal(sort(spectraVariables(res)),
                 sort(unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS), "rtime",
                               "mz", "intensity"))))
    res$new_col <- "b"
    expect_equal(sort(spectraVariables(res)),
                 sort(unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS), "rtime",
                               "mz", "intensity", "new_col"))))
    res_2 <- selectSpectraVariables(res, spectraVariables(res))
    expect_equal(spectraVariables(res), spectraVariables(res_2))
    res <- selectSpectraVariables(res, c("rtime"))
    expect_equal(
        spectraVariables(res),
        unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS)), "rtime"))
})

test_that("[,MsBackendWeizMass works", {
    be <- backendInitialize(MsBackendWeizMass(), dbc)

    res <- be[2:1]
    expect_true(length(res) == 2)
    expect_equal(res@spectraIds, be@spectraIds[2:1])
    expect_equal(res@spectraIds, res$spectrumId)
    expect_equal(be$mz[2:1], res$mz)

    res <- be[c(2, 1, 2, 2, 1)]
    expect_equal(mz(res)[[1]], mz(be)[[2]])
    expect_equal(mz(res)[[2]], mz(be)[[1]])
    expect_equal(mz(res)[[3]], mz(be)[[2]])
    expect_equal(mz(res)[[4]], mz(be)[[2]])
    expect_equal(mz(res)[[5]], mz(be)[[1]])
    expect_equal(intensity(res)[[1]], intensity(be)[[2]])
    expect_equal(intensity(res)[[2]], intensity(be)[[1]])
    expect_equal(intensity(res)[[3]], intensity(be)[[2]])
    expect_equal(intensity(res)[[4]], intensity(be)[[2]])
    expect_equal(intensity(res)[[5]], intensity(be)[[1]])
    expect_equal(res$spectrumId[1], be$spectrumId[2])
    expect_equal(res$spectrumId[2], be$spectrumId[1])
    expect_equal(res$spectrumId[3], be$spectrumId[2])
    expect_equal(res$spectrumId[4], be$spectrumId[2])
    expect_equal(res$spectrumId[5], be$spectrumId[1])
})

test_that("lengths,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    expect_equal(lengths(be), integer())

    be <- backendInitialize(be, dbc)
    res <- lengths(be)
    expect_true(length(res) > 0)
    expect_true(all(res > 0))
})

test_that("intensity,intensity<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- intensity(be)
    expect_equal(res, IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, dbc)
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(intensity(be) <- 3, "Can not")

    ## duplicated spectra
    be2 <- be[c(2, 1, 2, 1, 2, 1)]
    res2 <- intensity(be2)
    expect_equal(res[[2]], res2[[1]])
    expect_equal(res[[1]], res2[[2]])
    expect_equal(res[[2]], res2[[3]])
    expect_equal(res[[1]], res2[[4]])
    expect_equal(res[[2]], res2[[5]])
    expect_equal(res[[1]], res2[[6]])
})

test_that("mz,mz<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- mz(be)
    expect_equal(res, IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, dbc)
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(mz(be) <- 3, "Can not")

    ## duplicated spectra.
    be2 <- be[c(2, 1, 2, 1, 2, 1)]
    res2 <- mz(be2)
    expect_equal(res[[2]], res2[[1]])
    expect_equal(res[[1]], res2[[2]])
    expect_equal(res[[2]], res2[[3]])
    expect_equal(res[[1]], res2[[4]])
    expect_equal(res[[2]], res2[[5]])
    expect_equal(res[[1]], res2[[6]])
})

test_that("ionCount,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- ionCount(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- ionCount(be)
    expect_true(is.numeric(res))
    expect_equal(length(res), length(be))
    expect_true(all(res > 0))
})

test_that("isEmpty,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- isEmpty(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- isEmpty(be)
    expect_true(all(!res))
})

test_that("isolationWindowLowerMz,&<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- isolationWindowLowerMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowLowerMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowLowerMz(be) <- seq_len(length(be))
    res <- isolationWindowLowerMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowLowerMz(be) <- 1:3, "length 1")
    expect_error(isolationWindowLowerMz(be) <- "a", "numeric")
})

test_that("isolationWindowTargetMz,&<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- isolationWindowTargetMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowTargetMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowTargetMz(be) <- seq_len(length(be))
    res <- isolationWindowTargetMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowTargetMz(be) <- 1:3, "length 1")
    expect_error(isolationWindowTargetMz(be) <- "a", "numeric")
})

test_that("isolationWindowUpperMz,&<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- isolationWindowUpperMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowUpperMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowUpperMz(be) <- seq_len(length(be))
    res <- isolationWindowUpperMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowUpperMz(be) <- 1:3, "length 1")
    expect_error(isolationWindowUpperMz(be) <- "a", "numeric")
})

test_that("polarity,&<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- polarity(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- polarity(be)
    expect_true(all(res %in% c(0L, 1L)))

    polarity(be) <- rep(c(1, 0), length(be)/2)
    res <- polarity(be)
    expect_equal(res, rep(c(1, 0), length(be)/2))

    expect_error(polarity(be) <- 1:3, "length 1")
    expect_error(polarity(be) <- "a", "integer")
})

test_that("precursorCharge,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- precursorCharge(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- precursorCharge(be)
    expect_equal(res, rep(NA_integer_, length(be)))

    be$precursorCharge <- 1L
    res <- precursorCharge(be)
    expect_equal(res, rep(1L, length(be)))
})

test_that("precursorIntensity,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- precursorIntensity(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- precursorIntensity(be)
    expect_true(is.numeric(res))

    nmbrs <- abs(rnorm(length(be)))
    be$precursorIntensity <- nmbrs
    res <- precursorIntensity(be)
    expect_equal(res, nmbrs)
})

test_that("precursorMz,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- precursorMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- precursorMz(be)
    expect_true(is.numeric(res))

    be$precursorMz <- 12.21
    res <- precursorMz(be)
    expect_equal(res, rep(12.21, length(be)))
})

test_that("reset,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- reset(be)
    expect_equal(res, be)

    be <- backendInitialize(be, dbc)
    res <- reset(be)
    expect_equal(res, be)

    be$new_col <- "c"
    be <- be[c(2, 1)]
    res <- reset(be)
    expect_true(!any(spectraVariables(res) == "new_col"))
})

test_that("rtime,rtime<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- rtime(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- rtime(be)
    expect_true(!any(is.na(res)))

    rtime(be) <- 1.4
    res <- rtime(be)
    expect_equal(res, rep(1.4, length(be)))

    expect_error(rtime(be) <- 1:3, "length 1")
    expect_error(rtime(be) <- "a", "numeric")
})

test_that("scanIndex,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- scanIndex(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- scanIndex(be)
    expect_equal(res, rep(NA_integer_, length(be)))

    be$scanIndex <- seq_len(length(be))
    res <- scanIndex(be)
    expect_equal(res, seq_len(length(be)))
})

test_that("smoothed,smoothed<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- smoothed(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- smoothed(be)
    expect_equal(res, rep(NA, length(be)))

    smoothed(be) <- rep(c(TRUE, FALSE), length(be)/2)
    res <- smoothed(be)
    expect_equal(res, rep(c(TRUE, FALSE), length(be)/2))

    expect_error(smoothed(be) <- c(TRUE, FALSE, TRUE), "length 1")
    expect_error(smoothed(be) <- "a", "logical")
})

test_that("spectraData<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    expect_error(spectraData(be) <- spectraData(be), "not support")
})

test_that("spectraNames,spectraNames<-,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- spectraNames(be)
    expect_equal(res, character())

    be <- backendInitialize(be, dbc)
    res <- spectraNames(be)
    expect_true(is.character(res))
    expect_true(length(res) == length(be))
    expect_equal(res, as.character(be@spectraIds))

    expect_error(spectraNames(be) <- "a", "not support")
})

test_that("tic,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- tic(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- tic(be)
    expect_true(is.numeric(res))
    expect_true(all(is.na(res)))

    be$totIonCurrent <- 1.12
    res <- tic(be)
    expect_equal(res, rep(1.12, length(be)))

    res <- tic(be, initial = FALSE)
    expect_true(is.numeric(res))
    expect_true(length(res) == length(be))
    expect_true(all(res > 0))
})

test_that("msLevel,MsBackendWeizMass works", {
    be <- MsBackendWeizMass()
    res <- msLevel(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- msLevel(be)
    expect_true(is.integer(res))
    expect_true(length(res) == length(be))
    expect_true(all(res %in% c(NA_integer_, 2L)))
})

test_that("precursorMz,Spectra works", {
    sps <- Spectra(backendInitialize(MsBackendWeizMass(), dbc))
    precursors_cached <- precursorMz(sps)
    sps@backend@localData$precursorMz <- NULL
    precursors_direct <- precursorMz(sps)
    expect_equal(precursors_cached, precursors_direct)
})
