#' @alisases Organism.dplyr-depricated
#'
#' @title Deprecated functionality
#'
#' @description
#'
#' @name Depricated
NULL

#' @description \code{BasicFilter} will initiate a \code{AnnotationFilter}
#' object.  Use \code{AnnotationFilter} instead.
#'
#' @rdname Deprecated
BasicFilter <- function(){
	.Depricated("AnnotationFilter")
	
	return(new("AnnotationFilter"))
}

setMethod("initialize", "Foo", 


#' @description \code{GRangesFilter} will initiate a \code
