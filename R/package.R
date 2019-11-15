#' Tools for faster raster manipulation
#'
#' \code{fraster} provides tools for performing faster raster manipulations and
#' calculations.
#'
#' @details
#'
#' @author D. Scott Rinnan
#'
#' @seealso \code{\link[raster]{raster-package}}, \code{\link[sf]{sf}}, \code{\link[sp]{sp}}
#'
#' @docType package
#' @name fraster-package
## @useDynLib fraster
#' @import raster
#' @import fasterize
#' @importClassesFrom sp SpatialPolygons SpatialPolygonsDataFrame SpatialPoints SpatialPointsDataFrame
#' @importFrom sf st_sf st_sfc
#' @import methods
## @importFrom Rcpp sourceCpp
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nSupport women in science: 500womenscientists.org")
}

.onUnload <- function (libpath) {
  closeAllConnections()
}
