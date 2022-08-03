#' @rdname MsBackendWeizMass
#'
#' @export MsBackendWeizMass
MsBackendWeizMass <- function() {
    if (!requireNamespace("DBI", quietly = TRUE))
        stop("'MsBackendWeizMass' requires package 'DBI'. Please ",
             "install with 'install.packages(\"DBI\")'")
    new("MsBackendWeizMass")
}

#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_dbcon <- function(x) {
    if (length(x)) {
        if (!inherits(x, "DBIConnection"))
            return("'dbcon' is expected to be a connection to a database")
        tables <- dbListTables(x)
        if (!all(c("RECORD", "SPECTRUM", "PEAK") %in% tables))
            return(paste0("Database lacks some required tables. Is 'dbcon' a",
                          " connection to a WeizMass database?"))
    }
    NULL
}

#' Returns the spectra data, from the database and eventually filling with
#' *core* spectra variables, if they are not available in the database.
#'
#' The data can be either:
#' - in the database.
#' - in the local data (if new variables were added with $name <-).
#' - core spectra variables - if they are not in the database they have to be
#'   initialized with `NA` and the correct data type.
#'
#' @return a `data.frame` - always, even if only with a single column.
#'
#' @importFrom IRanges NumericList CharacterList
#'
#' @importFrom S4Vectors extractCOLS
#'
#' @importFrom S4Vectors make_zero_col_DFrame
#'
#' @importFrom methods as callNextMethod getMethod
#'
#' @importMethodsFrom Spectra spectraVariables
#'
#' @author Johannes Rainer
#'
#' @noRd
.spectra_data_weizmass_sql <- function(x, columns = spectraVariables(x)) {
    res <- getMethod("spectraData", "MsBackendCached")(x, columns = columns)
    if (is.null(res))
        res <- make_zero_col_DFrame(length(x))
    ## Define what needs to be still retrieved.
    db_cols <- intersect(columns, x@spectraVariables)
    db_cols <- db_cols[!db_cols %in% c("mz", "intensity",
                                       "relative_intensity", "peak_annotation",
                                       colnames(res))]
    mz_cols <- intersect(columns, c("mz", "intensity",
                                    "relative_intensity", "peak_annotation"))
    if (length(db_cols)) {
        res <- cbind(
            res, as(.fetch_spectra_data_sql(x, columns = db_cols), "DataFrame"))
        if (any(colnames(res) == "common_name"))
            res$common_name <- CharacterList(res$common_name, compress = FALSE)
        if (any(colnames(res) == "iupac_name"))
            res$iupac_name <- CharacterList(res$iupac_name, compress = FALSE)
    }
    ## Get m/z and intensity values
    if (length(mz_cols)) {
        pks <- .fetch_peaks_sql(x, columns = mz_cols)
        f <- as.factor(pks$spectrum_id)
        if (any(mz_cols == "mz")) {
            mzs <- unname(
                split(as.numeric(pks$mz), f)[as.character(x@spectraIds)])
            res$mz <- NumericList(mzs, compress = FALSE)
        }
        if (any(mz_cols == "intensity")) {
            ints <- unname(
                split(as.numeric(pks$intensity), f)[as.character(x@spectraIds)])
            res$intensity <- NumericList(ints, compress = FALSE)
        }
        if (any(mz_cols == "relative_intensity")) {
            ints <- unname(split(as.numeric(pks$relative_intensity), f)[
                as.character(x@spectraIds)])
            res$relative_intensity <- NumericList(ints, compress = FALSE)
        }
        if (any(mz_cols == "peak_annotation")) {
            ints <- unname(split(pks$peak_annotation, f)[
                as.character(x@spectraIds)])
            res$peak_annotation <- CharacterList(ints, compress = FALSE)
        }
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
    extractCOLS(res, columns)
}

#' @importFrom DBI dbSendQuery dbBind dbFetch dbClearResult
#'
#' @noRd
.fetch_peaks_sql <- function(x, columns = c("mz", "intensity")) {
    if (length(x@dbcon)) {
        cols <- .map_spectraVariables_to_sql(columns)
        res <- dbGetQuery(
            x@dbcon,
            paste0("select SPECTRUM_ID,", paste(cols, collapse = ","),
                   " from PEAK where SPECTRUM_ID in (",
                   paste0(unique(x@spectraIds), collapse = ","),")"))
        colnames(res) <- c("spectrum_id", columns)
        if (any(columns == "mz"))
            res$mz <- as.numeric(res$mz)
        if (any(columns == "intensity"))
            res$intensity <- as.numeric(res$intensity)
        if (any(columns == "relative_intensity"))
            res$relative_intensity <- as.numeric(res$relative_intensity)
        res
    } else {
        data.frame(spectrum_id = character(), mz = numeric(),
                   intensity = numeric(), relative_intensity = numeric(),
                   peak_annotation = character())[, c("spectrum_id", columns)]
    }
}

#' Fetches the m/z and intensity values from the database and returns a list
#' of two column matrices (m/z, intensity). The function ensures that the data
#' is returned in the same order than x@spectraIds (also allowing duplicated
#' entries).
#'
#' @param x `MsBackendWeizMass`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_data <- function(x, columns = c("mz", "intensity")) {
    p <- .fetch_peaks_sql(x, columns = columns)
    p <- unname(split.data.frame(
        p, as.factor(p$spectrum_id))[as.character(x@spectraIds)])
    emat <- matrix(ncol = length(columns), nrow = 0,
                   dimnames = list(character(), columns))
    idx <- seq(2, (length(columns) + 1L))
    if (length(idx) == 1) {
        lapply(p, function(z) {
            if (nrow(z))
                matrix(z[, idx], dimnames = list(c(), columns))
            else emat
        })
    } else {
        lapply(p, function(z) {
            if (nrow(z))
                as.matrix(z[, idx], rownames.force = FALSE)
            else emat
        })
    }
}

.columns_sql <- c(
    spectrumId = "SPECTRUM_ID",
    polarity = "ION_MODE",
    adduct = "PRECURSOR_ION",
    rtime_ci = "RT_CI",
    rtime = "RT",
    dataOrigin = "SAMPLE_FILE",
    formula = "FORMULA",
    exactmass = "EXACT_MASS",
    smiles = "SMILES",
    inchikey = "INCHIKEY",
    precursorMz = "PRECURSOR_MZ",
    compound_id = "SPECTRUM.ID",
    mz = "MZ",
    relative_intensity = "RELATIVE_INTENSITY",
    intensity = "INTENSITY",
    peak_annotation = "ANNOTATION",
    common_name = "COMMON_NAME",
    iupac_name = "IUPAC_NAME",
    instrument = "INSTRUMENT"
)

.map_spectraVariables_to_sql <- function(x) {
    for (i in seq_along(.columns_sql))
        x[x == names(.columns_sql)[i]] <- .columns_sql[i]
    x
}

.map_sql_to_spectraVariables <- function(x) {
    for (i in seq_along(.columns_sql))
        x[x == .columns_sql[i]] <- names(.columns_sql[i])
    x
}

#' Simple helper that creates a join query depending on the provided columns.
#'
#' @param x `Spectra`.
#'
#' @param columns `character` with the column names.
#'
#' @noRd
.join_query <- function(x, columns) {
    res <- "SPECTRUM"
    if (any(columns %in% x@.tables$RECORD))
        res <- paste0(res, " JOIN RECORD ON (SPECTRUM.ID",
                      "=RECORD.ID)")
    res
}

.fetch_spectra_data_sql <- function(x, columns = c("spectrum_id")) {
    orig_columns <- columns
    if (any(columns %in% c("common_name", "iupac_name"))) {
        columns <- columns[!columns %in% c("common_name", "iupac_name")]
        columns <- unique(c(columns, "compound_id"))
    }
    sql_columns <-
        unique(c("SPECTRUM_ID", .map_spectraVariables_to_sql(columns)))
    ## That turns out to be faster than dbBind if we use a field in the
    ## database that is unique (such as spectrum_id).
    res <- dbGetQuery(
        x@dbcon,
        paste0("SELECT ", paste(sql_columns, collapse = ","), " FROM ",
               .join_query(x, sql_columns), " WHERE SPECTRUM_ID IN (",
               paste0(unique(x@spectraIds), collapse = ", ") ,")"))
    idx <- match(x@spectraIds, res$SPECTRUM_ID)
    res <- res[idx[!is.na(idx)], , drop = FALSE]
    rownames(res) <- NULL
    colnames(res) <- .map_sql_to_spectraVariables(colnames(res))
    colnames(res)[colnames(res) == "ID"] <- "compound_id"
    if (any(columns == "polarity")) {
        pol <- rep(NA_integer_, nrow(res))
        pol[res$polarity == "POSITIVE"] <- 1L
        pol[res$polarity == "NEGATIVE"] <- 0L
        res$polarity <- pol
    }
    ## So far we're not dealing with multiple precursor m/z here!
    if (any(columns == "precursorMz"))
        suppressWarnings(
            res$precursorMz <- as.numeric(res$precursorMz))
    ## manage synonym and compound_name. Need a second query for that.
    if (any(orig_columns %in% c("common_name", "iupac_name"))) {
        cmps <- dbGetQuery(
            x@dbcon,
            paste0("SELECT * FROM CH_NAME WHERE CH_NAME.ID IN (",
                   paste0("'", unique(res$compound_id), "'",
                          collapse = ","), ")"))
        if (any(orig_columns == "common_name"))
            res$common_name <- split(
                cmps$COMMON_NAME,
                as.factor(cmps$ID))[as.character(res$compound_id)]
        if (any(orig_columns == "iupac_name"))
            res$iupac_name <- split(
                cmps$IUPAC_NAME,
                as.factor(cmps$ID))[as.character(res$compound_id)]
    }
    res[, orig_columns, drop = FALSE]
}
