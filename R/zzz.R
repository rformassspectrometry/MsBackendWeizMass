.onLoad <- function(libname, pkgname) {
    options(.PEAKS_FUN = .fetch_peaks_sql_order_mz)
    if (length(getOption("NO_ORDER_MZ")))
        options(.PEAKS_FUN = .fetch_peaks_sql)
}
