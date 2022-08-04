library(testthat)
library(MsBackendWeizMass)

library(RSQLite)
dbc <- dbConnect(SQLite(), system.file("sqlite", "weizmassv2.sqlite",
                                       package = "MsBackendWeizMass"))

test_check("MsBackendWeizMass")

## Run additional tests from Spectra:
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")

be <- MsBackendWeizMass()
be <- backendInitialize(be, dbc)
res <- test_file(paste0(test_suite, "/test_spectra_variables.R"),
                 reporter = check_reporter(), stop_on_failure = TRUE)
