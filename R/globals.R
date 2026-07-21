#' Internal package imports and global-variable declarations
#'
#' These declarations prevent false-positive R CMD check notes for
#' non-standard evaluation columns used in dplyr pipelines and spatial helpers.
#'
#' @keywords internal
#' @importFrom CoordinateCleaner clean_coordinates
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom utils tail
"_PACKAGE"

utils::globalVariables(
  c(
    ".bf_rowid",
    ".bf_sinas_location_id",
    ".data",
    "X",
    "Y",
    "category",
    "decimalLatitude",
    "decimalLongitude",
    "zip_path"
  )
)
