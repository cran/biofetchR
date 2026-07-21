#' Report optional biofetchR dependencies
#'
#' Reports whether optional packages used by extended biofetchR workflows are
#' available. This function reports package availability only. Optional
#' packages should be installed by the user outside package examples, tests
#' and vignettes.
#'
#' @return A data frame with one row per optional package and columns:
#'   `package`, `available` and `purpose`.
#'
#' @examples
#' install_optional_deps()
#'
#' @export
install_optional_deps <- function() {
  optional <- data.frame(
    package = c(
      "geosphere",
      "ggplot2",
      "httr",
      "lwgeom",
      "mapme.biodiversity",
      "osmdata",
      "readxl",
      "rnaturalearth",
      "stringi",
      "stringr",
      "terra",
      "tidyr",
      "worrms"
    ),
    purpose = c(
      "Distance calculations and spatial utilities",
      "Plotting examples and visual summaries",
      "HTTP requests for optional web workflows",
      "Geometry repair and advanced spatial operations",
      "Optional biodiversity context layers",
      "OpenStreetMap-based contextual overlays",
      "Reading spreadsheet inputs",
      "Natural Earth contextual overlays",
      "String encoding utilities",
      "String manipulation utilities",
      "Raster context workflows",
      "Data reshaping utilities",
      "World Register of Marine Species workflows"
    ),
    stringsAsFactors = FALSE
  )

  optional$available <- vapply(
    optional$package,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )

  optional
}

