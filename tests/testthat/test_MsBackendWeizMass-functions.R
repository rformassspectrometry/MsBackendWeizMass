test_that(".valid_dbcon works", {
    expect_true(length(.valid_dbcon(dbc)) == 0)
    expect_true(length(.valid_dbcon(4)) == 1)
})

test_that("MsBackendWeizMass works", {
    res <- MsBackendWeizMass()
    expect_true(validObject(res))
    expect_true(is(res, "MsBackendWeizMass"))
})

test_that(".fetch_peaks_sql works", {
    be <- MsBackendWeizMass()
    res <- .fetch_peaks_sql(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- .fetch_peaks_sql(be, columns = c("mz", "intensity"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) > 0)
    expect_true(is.numeric(res$mz))
    expect_true(is.numeric(res$intensity))

    res <- .fetch_peaks_sql(be, columns = peaksVariables(be))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity",
                                  "relative_intensity", "peak_annotation"))
    expect_true(nrow(res) > 0)
    expect_true(is.numeric(res$mz))
    expect_true(is.numeric(res$intensity))
    expect_true(is.character(res$peak_annotation))
})

test_that(".fetch_peaks_sql_order_mz works", {
    be <- MsBackendWeizMass()
    res <- .fetch_peaks_sql_order_mz(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- .fetch_peaks_sql_order_mz(be, columns = c("mz", "intensity"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) > 0)
    expect_true(is.numeric(res$mz))
    expect_true(is.numeric(res$intensity))
    expect_false(is.unsorted(res$mz))

    res <- .fetch_peaks_sql_order_mz(be, columns = peaksVariables(be))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity",
                                  "relative_intensity", "peak_annotation"))
    expect_true(nrow(res) > 0)
    expect_true(is.numeric(res$mz))
    expect_true(is.numeric(res$intensity))
    expect_true(is.character(res$peak_annotation))
    expect_false(is.unsorted(res$mz))
})

test_that(".peaks_data works", {
    be <- MsBackendWeizMass()
    res <- .peaks_data(be)
    expect_true(length(res) == 0)

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- .peaks_data(be)
    expect_true(is.list(res))
    expect_true(length(res) == length(be))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("mz", "intensity"))
    expect_equal(colnames(res[[2L]]), c("mz", "intensity"))
    expect_true(is.numeric(res[[1L]]))
    expect_true(is.numeric(res[[2L]]))

    be2 <- be[c(2, 1, 2)]
    res2 <- .peaks_data(be2)
    expect_true(length(res2) == length(be2))
    expect_equal(res[[1L]], res2[[2L]])
    expect_equal(res[[2L]], res2[[1L]])
    expect_equal(res[[2L]], res2[[3L]])

    res <- .peaks_data(be, columns = peaksVariables(be))
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_true(is.character(res[[1L]]))
    expect_true(is.character(res[[2L]]))
    expect_equal(colnames(res[[1L]]),
                 c("mz", "intensity", "relative_intensity", "peak_annotation"))
    expect_equal(colnames(res[[2L]]),
                 c("mz", "intensity", "relative_intensity", "peak_annotation"))
})

test_that(".fetch_spectra_data_sql works", {
    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- .fetch_spectra_data_sql(be, c("rtime"))
    expect_true(nrow(res) > 0)
    expect_true(colnames(res) == "rtime")

    res <- .fetch_spectra_data_sql(be, c("rtime", "spectrumId"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("rtime", "spectrumId"))

    res <- .fetch_spectra_data_sql(be, columns = "polarity")
    expect_equal(res$polarity, c(1L, 0L))

    ## synonym and compound_name
    res <- .fetch_spectra_data_sql(
        be, columns = c("rtime", "inchikey",
                        "common_name"))
    expect_equal(colnames(res),
                 c("rtime", "inchikey", "common_name"))
    expect_true(is.list(res$common_name))


    res <- .fetch_spectra_data_sql(
        be, columns = c("precursorMz",
                        "common_name", "iupac_name"))
    expect_equal(colnames(res),
                 c("precursorMz", "common_name", "iupac_name"))
    expect_true(is.list(res$common_name))
    expect_true(is.list(res$iupac_name))

    res <- .fetch_spectra_data_sql(be, columns = c("iupac_name"))
    expect_equal(colnames(res), "iupac_name")
    expect_true(is.list(res$iupac_name))

    res <- .fetch_spectra_data_sql(
        be, columns = c("dataOrigin", "ION", "adduct",
                        "exactmass", "formula", "spectrumId"))
    expect_equal(res$ION, c("MSE", "MSE"))
    expect_equal(res$adduct, c("[M]+", "[M-H]-"))
    expect_true(is.numeric(res$exactmass))
    expect_true(is.character(res$formula))
    expect_equal(res$spectrumId, be@spectraIds)
})

test_that(".spectra_data_weizmass_sql works", {
    be <- MsBackendWeizMass()
    res <- .spectra_data_weizmass_sql(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0L)
    expect_equal(colnames(res), names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    ## Full data.
    expect_error(.spectra_data_weizmass_sql(be, columns = "other"), "available")
    res <- .spectra_data_weizmass_sql(be)
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), spectraVariables(be))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))
    expect_true(is.numeric(res$rtime))

    ## A single column
    res <- .spectra_data_weizmass_sql(be, columns = "rtime")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "rtime")
    expect_true(all(!is.na(res$rtime)))

    res <- .spectra_data_weizmass_sql(be, "msLevel")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "msLevel")
    expect_true(all(res$msLevel %in% c(NA, 2L)))

    res <- .spectra_data_weizmass_sql(be, "mz")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "mz")
    expect_true(is(res$mz, "NumericList"))

    ## A combination of core and db
    res <- .spectra_data_weizmass_sql(be, columns = c("msLevel", "polarity"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("msLevel", "polarity"))
    expect_true(all(res$polarity %in% c(0L, 1L)))
    expect_true(all(is.na(res$msLevel)))

    res <- .spectra_data_weizmass_sql(
        be, columns = c("mz", "polarity", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("mz", "polarity", "intensity"))
    expect_true(all(res$polarity %in% c(0L, 1L)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    ## With local data.
    be@localData <- data.frame(rtime = seq_len(length(be)))
    res <- .spectra_data_weizmass_sql(
        be, columns = c("mz", "rtime", "spectrumId"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(res$rtime, seq_len(length(be)))

    be@localData$authors <- rep("A", length(be))
    res <- .spectra_data_weizmass_sql(be, columns = c("rtime", "authors"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(res$authors, rep("A", length(be)))
})

test_that(".join_query works", {
    ## Without any tables.
    be <- MsBackendWeizMass()
    res <- .join_query(be, .map_spectraVariables_to_sql(
                               c("compound_id", "spectrumId", "inchikey")))
    expect_equal(res, "SPECTRUM")

    be <- backendInitialize(MsBackendWeizMass(), dbc)
    res <- .join_query(be, .map_spectraVariables_to_sql(
                               c("compound_id", "spectrumId", "inchikey")))
    expect_equal(res, "SPECTRUM JOIN RECORD ON (SPECTRUM.ID=RECORD.ID)")

    res <- .join_query(
        be, .map_spectraVariables_to_sql(c("spectrumId", "rtime")))
    expect_equal(res, "SPECTRUM")
})
