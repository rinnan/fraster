#' Fast rasterize
#'
#' Expands the fasterize() function to sp compatibility.
#'
#'
#' @param x Spatial* object from sf or sp classes
#' @param y Raster* object
#' @param field numeric or character. The value(s) to be transferred. This can
#' be a single number, or a vector of numbers that has the same length as the
#' number of spatial features (points, lines, polygons). If x is a Spatial*DataFrame,
#' this can be the column name of the variable to be transferred. If missing, the
#' attribute index is used (i.e. numbers from 1 to the number of features). You
#' can also provide a vector with the same length as the number of spatial features,
#' or a matrix where the number of rows matches the number of spatial features
#' @param fun function or character. To determine what values to assign to cells
#' that are covered by multiple spatial features. You can use functions such as
#' min, max, or mean, or one of the following character values: 'first', 'last',
#' 'count'. The default value is 'last'. In the case of SpatialLines*, 'length'
#' is also allowed (currently for planar coordinate systems only).
#'
#' If x represents points, fun must accept a na.rm argument, either explicitly or
#' through 'dots'. This means that fun=length fails, but fun=function(x,...)length(x)
#' works, although it ignores the na.rm argument. To use the na.rm argument you can
#' use a function like this:
#' fun=function(x, na.rm)if (na.rm) length(na.omit(x)) else (length(x), or use a
#' function that removes NA values in all cases, like this function to compute the
#' number of unique values per grid cell "richness":
#' fun=function(x, ...) {length(unique(na.omit(x)))} . If you want to count the
#' number of points in each grid cell, you can use fun='count' or
#' fun=function(x,...){length(x)}.
#'
#' You can also pass multiple functions using a statement like
#' fun=function(x, ...) c(length(x),mean(x)), in which case the returned object
#' is a RasterBrick (multiple layers).
#' @param background numeric. Value to put in the cells that are not covered by
#' any of the features of x. Default is NA
#' @param by character.  The name of a column in `sf` by which to aggregate
#' layers.  If set, fasterize will return a RasterBrick with as many layers
#' as unique values of the `by` column.
#' @param ... additional arguments (none currently implemented)
#'
#' @examples
#'
#' @return A raster of the same size, extent, resolution and projection as the
#' provided raster template.
#'
#' @references Wylie, C., Romney, G., Evans, D., & Erdahl, A. (1967).
#'   Half-tone perspective drawings by computer. Proceedings of the November
#'   14-16, 1967, Fall Joint Computer Conference. AFIPS '67 (Fall).
#'   <https://dx.doi.org/10.1145/1465611.1465619>
#
#' @export
#' @name frasterize

setGeneric("frasterize", function(x, y, ...){
  standardGeneric("frasterize")})

#' @rdname frasterize
setMethod('frasterize', signature(x = 'sf', y = 'Raster'),
          function(x, y, field = NULL, fun = "last", background = NA, by = NULL) {
            .Call('_fasterize_fasterize', PACKAGE = 'fasterize', x, y, field, fun, background, by)
          }
)

#' @rdname frasterize
setMethod('frasterize', signature(x = 'Spatial', y = 'Raster'),
          function(x, y, field = NULL, fun = "last", background = NA, by = NULL) {
            temp <- sf::st_as_sf(x)
            .Call('_fasterize_fasterize', PACKAGE = 'fasterize', temp, y, field, fun, background, by)
          }
)
