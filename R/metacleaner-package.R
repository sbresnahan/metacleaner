#' metacleaner: Automated Curation of DNA Barcode Reference Databases
#'
#' DNA barcode reference databases built from public repositories such as
#' NCBI GenBank contain mislabeled and contaminant sequences. metacleaner
#' finds and removes them.
#'
#' Most users only need three functions:
#' * [clean_reference()] runs the whole cleaning workflow;
#' * [write_results()] saves the cleaned database and a report;
#' * [build_reference_db()] downloads reference or contaminant databases from
#'   NCBI.
#'
#' Start with `vignette("metacleaner")`.
#'
#' @keywords internal
#' @useDynLib metacleaner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
