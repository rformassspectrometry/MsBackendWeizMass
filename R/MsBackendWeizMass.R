#' @title MS backend accessing the WeizMass MySQL database
#'
#' @aliases MsBackendWeizMass-class compounds
#'
#' @description
#'
#' The `MsBackendWeizMass` provides access to WeizMass mass spectrometry
#' libraries by directly accessing its MySQL/MariaDb database. In addition the
#' backend supports adding new spectra variables to the object, or to *locally*
#' change spectra variables (without changing the original values in the
#' database).
#'
#' Note that `MsBackendWeizMass` requires access to a WeizMass MySQL/MariaDB
#' database.
#'
#' Also, some of the fields in the WeizMass database are not directly compatible
#' with `Spectra`, as the data is stored as text instead of numeric. The
#' precursor m/z values are for example stored as character in the database, but
#' are converted to numeric during the data access. Thus, for spectra with
#' non-numeric values stored in that field an `NA` is reported.
#'
#' @param dbcon For `backendInitialize,MsBackendWeizMass`: SQL database
#'     connection to the WeizMass database.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param drop For `[`: not considered.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`).
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param name name of the variable to replace for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param object Object extending `MsBackendWeizMass`.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackendWeizMass`.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Supported Backend functions:
#'
#' The following functions are supported by the `MsBackendWeizMass`.
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `backendInitialize`: initialises the backend by retrieving the IDs of all
#'   spectra in the database. Parameter `dbcon` with the connection to the
#'   WeizMass MySQL database is required.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of spectra
#'   in `object` with the *data origin* of each spectrum. This could e.g. be
#'   the mzML file from which the data was read.
#'
#' - `dataStorage`: returns `"<WeizMass>"` for all spectra.
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`. Note that
#'   the collision energy description from WeizMass are provided as spectra
#'   variable `"collisionEnergyText"`.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl` th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `length`: returns the number of spectra in the object.
#'
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `peaksData` returns a `list` with the spectras' peak data. The length of
#'   the list is equal to the number of spectra in `object`. Each element of
#'   the list is a `matrix` with columns defined by parameter `columns` which
#'   defaults to `columns = c("mz", "intensity")` but any of
#'   `peaksVariables(object)` would be supported.
#'   Note that if `columns` contains `"peak_annotation"`, the whole matrix will
#'   be of type `character` (i.e. even the m/z and intensity values will be
#'   provided as text). See examples below for details. For an empty spectrum,
#'   a `matrix` with 0 rows is returned.
#'
#' - `peaksVariables` returns a `character` with the provided peaks variables
#'   (i.e. data available for each individual mass peak). These can be used in
#'   `peaksData` to retrieve the specified values.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   2 and above spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `reset`: restores the backend to its original state, i.e. deletes all
#'   locally modified data and reinitializes the backend to the full data
#'   available in the database.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum (in seconds). `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the
#'   spectrum as reported in the mzML file.
#'
#' - `selectSpectraVariables`: reduces the information within the backend to
#'   the selected spectra variables.
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`: gets general spectrum metadata (annotation, also called
#'   header).  `spectraData` returns a `DataFrame`. Note that replacing the
#'   spectra data with `spectraData<-` is not supported.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object`.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`. This should return **all** spectra variables which
#'   are present in `object`, also `"mz"` and `"intensity"` (which are by
#'   default not returned by the `spectraVariables,Spectra` method).
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' @section Not supported Backend functions:
#'
#' The following functions are not supported by the `MsBackendWeizMass` since
#' the original data can not be changed.
#'
#' `backendMerge`, `export`, `filterDataStorage`, `filterPrecursorScan`,
#' `peaksData<-`, `filterAcquisitionNum`, `intensity<-`, `mz<-`, `precScanNum`,
#' `spectraData<-`, `spectraNames<-`.
#'
#' @section Retrieving compound annotations for spectra:
#'
#' While compound annotations are also provided *via* the `spectraVariables` of
#' the backend, it would also be possible to use the `compounds` function on
#' a `Spectra` object (that uses a `MsBackendWeizMass` backend) to retrieve
#' compound annotations for the specific spectra.
#'
#' @name MsBackendWeizMass
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @references
#'
#' Shahaf N., Rogachev I., Heinig U., Meir S, Malitsky S, Battat M. et al.
#' (2016). The WEIZMASS spectra library for high-confidence metabolite
#' identification. Nature Communications 7:12423. \doi{10.1038/ncomms12423}.
#'
#' @md
#'
#' @exportClass MsBackendWeizMass
#'
#' @examples
#'
#' ## Create a connection to a database with WeizMass data - in the present
#' ## example we connect to a tiny SQLite database bundled in this package.
#' library(RSQLite)
#' con <- dbConnect(SQLite(), system.file("sqlite", "weizmassv2.sqlite",
#'     package = "MsBackendWeizMass"))
#'
#' ## Given that we have the connection to a WeizMass database we can
#' ## initialize the backend:
#' be <- backendInitialize(MsBackendWeizMass(), dbcon = con)
#' be
#'
#' ## List available peak variables
#' peaksVariables(be)
#'
#' ## Get peaks data; by default only m/z and intensity values are returned
#' peaksData(be)
#'
#' ## Get peaks data including peak annotations; note that now for each
#' ## spectrum a character matrix is returned!
#' res <- peaksData(be, columns = c("mz", "intensity", "peak_annotation"))
#' res[[1L]]
#'
#' ## Get the m/z values for all spectra
#' mz(be)
#'
#' ## annotations for the invidual peaks can be retrieved with
#' be$peak_annotation
#'
#' ## List available spectra variables
#' spectraVariables(be)
#'
#' ## Access MS level
#' msLevel(be)
#' be$msLevel
#'
#' ## Access m/z values
#' be$mz
#'
#' ## Access the full spectra data (including m/z and intensity values)
#' spectraData(be)
#'
#' ## Add a new spectra variable
#' be$new_variable <- "b"
#' be$new_variable
#'
#' ## Subset the backend
#' be_sub <- be[c(2, 1)]
#'
#' spectraNames(be)
#' spectraNames(be_sub)
NULL

#' @importClassesFrom DBI DBIConnection
setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importClassesFrom Spectra MsBackendCached
#'
#' @importClassesFrom S4Vectors DataFrame
setClass(
    "MsBackendWeizMass",
    contains = "MsBackendCached",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "integer",
        .tables = "list"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = integer(),
        .tables = list(),
        readonly = TRUE, version = "0.2"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackendWeizMass", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @importFrom DBI dbGetQuery
#'
#' @rdname MsBackendWeizMass
setMethod("backendInitialize", "MsBackendWeizMass", function(object,
                                                             dbcon, ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsBackendWeizMass'")
    msg <- .valid_dbcon(dbcon)
    object@dbcon <- dbcon
    if (length(msg))
        stop(msg)

    res <- dbGetQuery(
        dbcon, "SELECT SPECTRUM_ID, PRECURSOR_MZ FROM SPECTRUM")
    object@spectraIds <- as.integer(res[, "SPECTRUM_ID"])
    object@.tables <- list(
        SPECTRUM = colnames(
            dbGetQuery(dbcon, "SELECT * FROM SPECTRUM LIMIT 0")),
        RECORD = colnames(
            dbGetQuery(dbcon, "SELECT * FROM RECORD LIMIT 0")),
        CH_NAME = colnames(
            dbGetQuery(dbcon, "SELECT * FROM CH_NAME LIMIT 0")),
        PEAK = colnames(dbGetQuery(dbcon, "SELECT * FROM PEAK LIMIT 0")))
    ## Initialize cached backend\
    cns <- unique(unlist(object@.tables))
    cns[cns == "ID"] <- "SPECTRUM.ID"
    object <- callNextMethod(
        object, nspectra = length(object@spectraIds),
        spectraVariables = c(.map_sql_to_spectraVariables(cns),
                             "precursor_mz_text"))
    object@localData$precursor_mz_text <- res[, "PRECURSOR_MZ"]
    suppressWarnings(object@localData$precursorMz <-
                         as.numeric(res[, "PRECURSOR_MZ"]))
    validObject(object)
    object
})

#' @importMethodsFrom Spectra peaksVariables
#'
#' @exportMethod peaksVariables
#'
#' @rdname MsBackendWeizMass
setMethod("peaksVariables", "MsBackendWeizMass", function(object) {
    res <- MsBackendWeizMass:::.map_sql_to_spectraVariables(object@.tables$PEAK)
    res[res != "spectrumId"]
})

#' @importMethodsFrom Spectra peaksData
#'
#' @exportMethod peaksData
#'
#' @rdname MsBackendWeizMass
setMethod("peaksData", "MsBackendWeizMass",
          function(object, columns = c("mz", "intensity")) {
        .peaks_data(object, columns = columns)
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname MsBackendWeizMass
setMethod("dataStorage", "MsBackendWeizMass", function(object) {
    rep("<WeizMass>", length(object))
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendWeizMass
setReplaceMethod("intensity", "MsBackendWeizMass", function(object, value) {
    stop("Can not replace original intensity values in WeizMass.")
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendWeizMass
setReplaceMethod("mz", "MsBackendWeizMass", function(object, value) {
    stop("Can not replace original data in WeizMass.")
})

#' @exportMethod reset
#'
#' @importMethodsFrom Spectra reset backendInitialize
#'
#' @rdname MsBackendWeizMass
setMethod("reset", "MsBackendWeizMass", function(object) {
    message("Restoring original data ...", appendLF = FALSE)
    if (is(object@dbcon, "DBIConnection"))
        object <- backendInitialize(object, object@dbcon)
    message("DONE")
    object
})

#' @importMethodsFrom Spectra spectraData
#'
#' @exportMethod spectraData
#'
#' @rdname MsBackendWeizMass
setMethod(
    "spectraData", "MsBackendWeizMass",
    function(object, columns = spectraVariables(object)) {
        .spectra_data_weizmass_sql(object, columns = columns)
    })

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackendWeizMass
setMethod("spectraNames", "MsBackendWeizMass", function(object) {
    as.character(object@spectraIds)
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackendWeizMass
setReplaceMethod("spectraNames", "MsBackendWeizMass",
                 function(object, value) {
                     stop(class(object)[1],
                          " does not support replacing spectra names (IDs).")
})

#' @exportMethod tic
#'
#' @importMethodsFrom ProtGenerics tic
#'
#' @importFrom Spectra intensity
#'
#' @importFrom MsCoreUtils vapply1d
#'
#' @rdname MsBackendWeizMass
setMethod("tic", "MsBackendWeizMass", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@localData) == "totIonCurrent"))
            object@localData[, "totIonCurrent"]
        else rep(NA_real_, times = length(object))
    } else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @exportMethod [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @importFrom methods slot<-
#'
#' @importFrom S4Vectors extractROWS
#'
#' @rdname MsBackendWeizMass
setMethod("[", "MsBackendWeizMass", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    slot(x, "spectraIds", check = FALSE) <- x@spectraIds[i]
    x <- callNextMethod(x, i = i)
    x
})

#' @rdname MsBackendWeizMass
#'
#' @export
setReplaceMethod("$", "MsBackendWeizMass", function(x, name, value) {
    if (name %in% c("spectrum_id"))
        stop("Spectra IDs can not be changed.", call. = FALSE)
    callNextMethod()
})


#' @rdname MsBackendWeizMass
#'
#' @importMethodsFrom Spectra precScanNum
#' @export
setMethod("precScanNum", "MsBackendWeizMass", function(object) {
    warning("precursor scan numbers not available")
    rep(NA_integer_, length(object))
})
