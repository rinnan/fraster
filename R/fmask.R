#' Fast masking of values in a Raster object
#'
#' Create a new Raster* object that has the same values as x, except for the cells
#' that are NA (or other maskvalue) in a 'mask'. These cells become NA (or other
#' updatevalue). The mask can be either another Raster* object of the same extent
#' and resolution, or a Spatial* object (e.g. SpatialPolygons) in which case all
#' cells that are not covered by the Spatial object are set to updatevalue. You
#' can use inverse=TRUE to set the cells that are not NA (or other maskvalue) in
#' the mask, or not covered by the Spatial* object, to NA (or other updatvalue).
#'
#' @param x Raster* object
#' @param mask Raster* object or a Spatial* object from sf or sp classes
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
#' @name fmask

setGeneric("fmask", function(x, mask, ...){
  standardGeneric("fmask")})

#' @rdname fmask
setMethod('fmask', signature(x='Raster', mask='sf'),
          function(x, mask, filename="", inverse=FALSE, updatevalue=NA, updateNA=FALSE, ...){
            #if (inherits(mask, 'sf')) {
              m <- fasterize::fasterize(mask, x)
            #} else {
            #  m <- rasterize(mask, x, 1, silent=TRUE)
            #}
            mask(x, m, filename=filename, inverse=inverse, maskvalue=NA, updatevalue=updatevalue, ...)
          } )

#' @rdname fmask
setMethod('fmask', signature(x='Raster', mask='Spatial'),
          function(x, mask, filename="", inverse=FALSE, updatevalue=NA, updateNA=FALSE, ...){
            #if (inherits(mask, 'SpatialPolygons')) {
              m <- .fasterize(mask, x, values=rep(1,length(mask)))
            #} else {
            #  m <- rasterize(mask, x, 1, silent=TRUE)
            #}
            mask(x, m, filename=filename, inverse=inverse, maskvalue=NA, updatevalue=updatevalue, ...)
          } )

#' @rdname fmask
setMethod('fmask', signature(x='RasterLayer', mask='RasterLayer'),
          function(x, mask, filename="", inverse=FALSE, maskvalue=NA, updatevalue=NA, updateNA=FALSE, ...){

            maskvalue <- maskvalue[1]
            updatevalue <- updatevalue[1]

            compareRaster(x, mask)
            out <- .copyWithProperties(x)

            if ( canProcessInMemory(x, 3)) {

              x <- getValues(x)
              mask <- getValues(mask)
              if (is.na(maskvalue)) {
                if (updateNA) {
                  if (inverse) {
                    x[!is.na(mask)] <- updatevalue
                  } else {
                    x[is.na(mask)] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[!is.na(mask) & !is.na(x)] <- updatevalue
                  } else {
                    x[is.na(mask) & !is.na(x)] <- updatevalue
                  }
                }
              } else {
                if (updateNA) {
                  if (inverse) {
                    x[mask != maskvalue] <- updatevalue
                  } else {
                    x[mask == maskvalue] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[mask != maskvalue & !is.na(x)] <- updatevalue
                  } else {
                    x[mask == maskvalue & !is.na(x)] <- updatevalue
                  }
                }
              }
              x <- setValues(out, x)
              if (filename != '') {
                x <- writeRaster(x, filename, ...)
              }
              return(x)

            } else {

              if (filename=='') {
                filename <- rasterTmpFile()
              }

              out <- writeStart(out, filename=filename, ...)
              tr <- blockSize(out)
              pb <- pbCreate(tr$n, label='mask', ...)

              if (is.na(updatevalue)) {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m)] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m)] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m != maskvalue] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m==maskvalue] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }
              } else {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (updateNA) {
                    if (inverse) {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m != maskvalue] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    } else {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m==maskvalue] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    }
                  } else {
                    if (inverse) {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m != maskvalue & !is.na(v)] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    } else {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m==maskvalue & !is.na(v)] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    }
                  }
                }
              }


              pbClose(pb)
              out <- writeStop(out)
              return(out)
            }
          }
)

#' @rdname fmask
setMethod('fmask', signature(x='RasterStackBrick', mask='RasterLayer'),
          function(x, mask, filename="", inverse=FALSE, maskvalue=NA, updatevalue=NA, updateNA=FALSE, ...){

            compareRaster(x, mask)
            maskvalue <- maskvalue[1]
            updatevalue <- updatevalue[1]

            out <- .copyWithProperties(x)


            if (canProcessInMemory(x, nlayers(x)+4)) {

              x <- getValues(x)

              if (is.na(maskvalue)) {
                if (updateNA) {
                  if (inverse) {
                    x[!is.na(getValues(mask))] <- updatevalue
                  } else {
                    x[is.na(getValues(mask))] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[!is.na(getValues(mask)) & !is.na(x)] <- updatevalue
                  } else {
                    x[is.na(getValues(mask)) & !is.na(x)] <- updatevalue
                  }
                }
              } else {
                if (updateNA) {
                  if (inverse) {
                    x[getValues(mask) != maskvalue] <- updatevalue
                  } else {
                    x[getValues(mask) == maskvalue] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[getValues(mask) != maskvalue & !is.na(x)] <- updatevalue
                  } else {
                    x[getValues(mask) == maskvalue & !is.na(x)] <- updatevalue
                  }
                }
              }
              out <- setValues(out, x)
              if (filename != '') {
                out <- writeRaster(out, filename, ...)
              }
              return(out)

            } else {

              if ( filename=='') {
                filename <- rasterTmpFile()
              }

              out <- writeStart(out, filename=filename, ...)
              tr <- blockSize(out)
              pb <- pbCreate(tr$n, label='mask', ...)

              if (is.na(updatevalue)) {

                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m)] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m)] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m != maskvalue] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m == maskvalue] <- NA
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }


              } else {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {

                  if (updateNA) {
                    if (inverse) {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m != maskvalue] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    } else {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m == maskvalue] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    }

                  } else {

                    if (inverse) {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m != maskvalue & !is.na(v)] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    } else {
                      for (i in 1:tr$n) {
                        v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                        m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                        v[m == maskvalue & !is.na(v)] <- updatevalue
                        out <- writeValues(out, v, tr$row[i])
                        pbStep(pb, i)
                      }
                    }
                  }
                }
              }
              pbClose(pb)
              out <- writeStop(out)
              return(out)
            }
          }
)

#' @rdname fmask
setMethod('fmask', signature(x='RasterLayer', mask='RasterStackBrick'),
          function(x, mask, filename="", inverse=FALSE, maskvalue=NA, updatevalue=NA, updateNA=FALSE, ...){

            compareRaster(x, mask)

            out <- brick(mask, values=FALSE)
            maskvalue <- maskvalue[1]
            updatevalue <- updatevalue[1]

            if (canProcessInMemory(mask, nlayers(x)*2+2)) {

              x <- getValues(x)
              x <- matrix(rep(x, nlayers(out)), ncol=nlayers(out))

              if (updateNA) {

                if (is.na(maskvalue)) {
                  if (inverse) {
                    x[!is.na(getValues(mask))] <- updatevalue
                  } else {
                    x[is.na(getValues(mask))] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[getValues(mask)!=maskvalue] <- updatevalue
                  } else {
                    x[getValues(mask)==maskvalue] <- updatevalue
                  }
                }
              } else {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    x[!is.na(getValues(mask)) & !is.na(x)] <- updatevalue
                  } else {
                    x[is.na(getValues(mask)) & !is.na(x)] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[getValues(mask)!=maskvalue & !is.na(x)] <- updatevalue
                  } else {
                    x[getValues(mask)==maskvalue & !is.na(x)] <- updatevalue
                  }
                }
              }
              out <- setValues(out, x)
              if (filename != '') {
                out <- writeRaster(out, filename, ...)
              }
              return(out)

            } else {

              if ( filename=='') { filename <- rasterTmpFile() }
              out <- writeStart(out, filename=filename, ...)
              tr <- blockSize(out)
              pb <- pbCreate(tr$n, label='mask', ...)

              if (updateNA) {

                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m != maskvalue] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m == maskvalue] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }
              } else {

                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[!is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[is.na(m) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m != maskvalue & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      v <- matrix(rep(v, nlayers(out)), ncol=nlayers(out))
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[m == maskvalue & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }
              }

              pbClose(pb)
              out <- writeStop(out)
              return(out)
            }
          }
)


#' @rdname fmask
setMethod('fmask', signature(x='RasterStackBrick', mask='RasterStackBrick'),
          function(x, mask, filename="", inverse=FALSE, maskvalue=NA, updatevalue=NA, updateNA=FALSE, ...){


            nlx <- nlayers(x)
            nlk <- nlayers(mask)
            if ( nlx != nlk ) {
              if (nlx == 1) {
                x <- raster(x, 1)
                return(mask(x, mask, filename=filename, inverse=inverse, maskvalue=maskvalue, updatevalue=updatevalue, ...))
              }
              if (nlk == 1) {
                mask <- raster(mask, 1)
                return(mask(x, mask, filename=filename, inverse=inverse, maskvalue=maskvalue, updatevalue=updatevalue, ...))
              }

              if (! ((nlx > nlk) & (nlx %% nlk == 0)) ) {
                stop('number of layers of x and mask must be the same,\nor one of the two should be 1, or the number of layers of x\nshould be divisible by the number of layers of mask')
              }
            }

            updatevalue <- updatevalue[1]
            maskvalue <- maskvalue[1]

            compareRaster(x, mask)
            out <- brick(x, values=FALSE)
            ln <- names(x)
            names(out) <- ln

            if (canProcessInMemory(x, nlx*2)) {

              x <- getValues(x)

              if (updateNA) {

                if (is.na(maskvalue)) {
                  if (inverse) {
                    x[!is.na(as.vector(getValues(mask)))] <- updatevalue
                  } else {
                    x[is.na(as.vector(getValues(mask)))] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[as.vector(getValues(mask)) != maskvalue] <- updatevalue
                  } else {
                    x[as.vector(getValues(mask)) == maskvalue] <- updatevalue
                  }
                }


              } else {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    x[!is.na(as.vector(getValues(mask))) & !is.na(x)] <- updatevalue
                  } else {
                    x[is.na(as.vector(getValues(mask))) & !is.na(x)] <- updatevalue
                  }
                } else {
                  if (inverse) {
                    x[as.vector(getValues(mask)) != maskvalue  & !is.na(x)] <- updatevalue
                  } else {
                    x[as.vector(getValues(mask)) == maskvalue & !is.na(x)] <- updatevalue
                  }
                }
              }
              out <- setValues(out, x)
              if (filename != '') {
                out <- writeRaster(out, filename, ...)
                names(out) <- ln
              }
              return(out)

            } else {

              if ( filename=='') { filename <- rasterTmpFile() }

              out <- writeStart(out, filename=filename, ...)
              tr <- blockSize(out)
              pb <- pbCreate(tr$n, label='mask', ...)

              if (updateNA) {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(!is.na(m))] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(is.na(m))] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(m != maskvalue)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(m == maskvalue)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }
              } else {
                if (is.na(maskvalue)) {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(!is.na(m)) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(is.na(m)) &  !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                } else {
                  if (inverse) {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(m != maskvalue) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  } else {
                    for (i in 1:tr$n) {
                      v <- getValues( x, row=tr$row[i], nrows=tr$nrows[i] )
                      m <- getValues( mask, row=tr$row[i], nrows=tr$nrows[i] )
                      v[as.vector(m == maskvalue) & !is.na(v)] <- updatevalue
                      out <- writeValues(out, v, tr$row[i])
                      pbStep(pb, i)
                    }
                  }
                }
              }
              pbClose(pb)
              out <- writeStop(out)
              names(out) <- ln
              return(out)
            }
          }
)
