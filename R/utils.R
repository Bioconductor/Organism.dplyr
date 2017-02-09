#' Utilities used in examples, vignettes, and tests
#'
#' These functions are primarily for illustrating
#' funtionality. \code{hg38light()} and \code{mm10light()} provide
#' access to trimmed-down versions of Organism.dplyr data based
#' derived from the TxDb.Hsapiens.UCSC.hg38.knownGene and
#' TxDb.Mmusculus.UCSC.mm10.ensGene data bases.
#'
#' @return character(1) file path to the trimmed-down data base
#' @examples
#' hg38light()
#' mm10light()
#' @rdname utils
#' @export
hg38light <- function() {
    system.file(
        package="Organism.dplyr", "extdata", "light.hg38.knownGene.sqlite"
    )
}

#' @rdname utils
#' @export
mm10light <- function() {
    system.file(
        package="Organism.dplyr", "extdata", "light.mm10.ensGene.sqlite"
    )
}
