#' Functions for filter objects
#'
#' These functions are used to create filters for genomic extrators.
#'
#' All these filters except \code{GRangesFilter()} extend \code{BasicFilter}
#' class. Each filter takes value(s) from the corresponding column. For example,
#' \code{AccnumFilter()} takes value of accession number(s), which come from
#' column \code{accnum}.
#'
#' \code{GRangesFilter()} takes a \code{GRanges} object as filter, and returns
#' genomic extractors (\code{genes}, \code{transcripts}, etc.) that are
#' partially overlapping with the region.
#'
#' \code{possibleFilters()} lists all available filters for \code{src_organism}
#' object.
#'
#' @aliases AccnumFilter AliasFilter Cds_chromFilter Cds_idFilter Cds_nameFilter
#'     Cds_strandFilter EnsemblFilter EnsemblprotFilter EnsembltransFilter
#'     EntrezFilter EnzymeFilter EvidenceFilter EvidenceallFilter
#'     Exon_chromFilter Exon_idFilter Exon_nameFilter Exon_rankFilter
#'     Exon_strandFilter FlybaseFilter Flybase_cgFilter Flybase_protFilter 
#'     Gene_chromFilter Gene_strandFilter GenenameFilter GoFilter GoallFilter
#'     IpiFilter MapFilter MgiFilter OmimFilter OntologyFilter OntologyallFilter
#'     PfamFilter PmidFilter PrositeFilter RefseqFilter SymbolFilter
#'     Tx_chromFilter Tx_idFilter Tx_nameFilter Tx_strandFilter Tx_typeFilter
#'     UnigeneFilter UniprotFilter WormbaseFilter ZfinFilter Cds_startFilter
#'     Cds_endFilter Exon_startFilter Exon_endFilter Gene_startFilter
#'     Gene_endFilter Tx_startFilter Tx_endFilter
#'
#' @usage AccnumFilter(value, condition = "==")
#' AliasFilter(value, condition = "==")
#' Cds_chromFilter(value, condition = "==")
#' Cds_idFilter(value, condition = "==")
#' Cds_nameFilter(value, condition = "==")
#' Cds_strandFilter(value, condition = "==")
#' EnsemblFilter(value, condition = "==")
#' EnsemblprotFilter(value, condition = "==")
#' EnsembltransFilter(value, condition = "==")
#' EntrezFilter(value, condition = "==")
#' EnzymeFilter(value, condition = "==")
#' EvidenceFilter(value, condition = "==")
#' EvidenceallFilter(value, condition = "==")
#' Exon_chromFilter(value, condition = "==")
#' Exon_idFilter(value, condition = "==")
#' Exon_nameFilter(value, condition = "==")
#' Exon_rankFilter(value, condition = "==")
#' Exon_strandFilter(value, condition = "==")
#' FlybaseFilter(value, condition = "==")
#' Flybase_cgFilter(value, condition = "==")
#' Flybase_protFilter(value, condition = "==")
#' Gene_chromFilter(value, condition = "==")
#' Gene_strandFilter(value, condition = "==")
#' GenenameFilter(value, condition = "==")
#' GoFilter(value, condition = "==")
#' GoallFilter(value, condition = "==")
#' IpiFilter(value, condition = "==")
#' MapFilter(value, condition = "==")
#' MgiFilter(value, condition = "==")
#' OmimFilter(value, condition = "==")
#' OntologyFilter(value, condition = "==")
#' OntologyallFilter(value, condition = "==")
#' PfamFilter(value, condition = "==")
#' PmidFilter(value, condition = "==")
#' PrositeFilter(value, condition = "==")
#' RefseqFilter(value, condition = "==")
#' SymbolFilter(value, condition = "==")
#' Tx_chromFilter(value, condition = "==")
#' Tx_idFilter(value, condition = "==")
#' Tx_nameFilter(value, condition = "==")
#' Tx_strandFilter(value, condition = "==")
#' Tx_typeFilter(value, condition = "==")
#' UnigeneFilter(value, condition = "==")
#' UniprotFilter(value, condition = "==")
#' WormbaseFilter(value, condition = "==")
#' ZfinFilter(value, condition = "==")
#' Cds_startFilter(value, condition = "==")
#' Cds_endFilter(value, condition = "==")
#' Exon_startFilter(value, condition = "==")
#' Exon_endFilter(value, condition = "==")
#' Gene_startFilter(value, condition = "==")
#' Gene_endFilter(value, condition = "==")
#' Tx_startFilter(value, condition = "==")
#' Tx_endFilter(value, condition = "==")
#'
#' @param value Value of the filter. For \code{GRangesFilter} vaule should be a
#'     \code{GRanges} object.
#'
#' @param condition The condition to be used in filter for genomic extractors,
#'     one of "==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=".  For
#'     character values "==", "!=", "startsWith" and "endsWith" are allowed, for
#'     numeric values (\code{Cds_startFilter}, \code{Cds_endFilter},
#'     \code{Exon_startFilter}, \code{Exon_endFilter}, \code{Gene_startFilter},
#'     \code{Gene_endFilter}, \code{Tx_startFilter} and \code{Tx_endFilte}),
#'     "==", "!=", ">", ">=", "<" and "<=". Default condition is "==".
#'
#' @seealso \code{\link{src_organism}} for creating a \code{src_organism}
#'     object.
#'
#'     \code{\link[Organism.dplyr]{transcripts_tbl}} for generic functions
#'      to extract genomic features from a \code{src_organism} object.
#'
#'      \code{\link[Organism.dplyr]{select,src_organism-method}} for "select"
#'     interface on \code{src_organism} objects.
#'
#' @examples
#' ## filter by ensembl
#' EnsemblFilter("ENSG00000171862")
#'
#' ## filter by gene symbol start with "BRAC"
#' SymbolFilter("BRCA", "startsWith")
#'
#' ## filter by GRanges
#' GRangesFilter(as("chr10:87869000-87876000", "GRanges"))
#'
#' ## filter by transcript start position
#' Tx_startFilter(87863438,">")
#'
#' @export AccnumFilter AliasFilter Cds_chromFilter Cds_idFilter Cds_nameFilter
#' @export Cds_strandFilter EnsemblFilter EnsemblprotFilter EnsembltransFilter
#' @export EntrezFilter EnzymeFilter EvidenceFilter EvidenceallFilter
#' @export Exon_chromFilter Exon_idFilter Exon_nameFilter Exon_rankFilter
#' @export Exon_strandFilter FlybaseFilter Flybase_cgFilter Flybase_protFilter
#' @export Gene_chromFilter Gene_strandFilter GenenameFilter GoFilter
#' @export GoallFilter IpiFilter MapFilter MgiFilter OmimFilter OntologyFilter
#' @export OntologyallFilter PfamFilter PmidFilter PrositeFilter RefseqFilter
#' @export SymbolFilter Tx_chromFilter Tx_idFilter Tx_nameFilter Tx_strandFilter
#' @export Tx_typeFilter UnigeneFilter UniprotFilter WormbaseFilter ZfinFilter
#' @export Cds_startFilter Cds_endFilter Exon_startFilter Exon_endFilter
#' @export Gene_startFilter Gene_endFilter Tx_startFilter Tx_endFilter
#' @rdname BasicFilter
#' @importFrom methods new setClass slot
#' @export
setClass("BasicFilter",
         representation(
             "VIRTUAL",
             field="character",
             condition="character",
             value="ANY",
             .valueIsCharacter="logical"
         ),
         prototype=list(
             condition= "==",
             value=character(),
             .valueIsCharacter=TRUE
         )
)

#' @rdname BasicFilter
#' @export
setClass("GRangesFilter",
         slots = list(field="character",
                      value="GRanges"))

#' @rdname BasicFilter
#' @export
GRangesFilter <- function(value){
    new("GRangesFilter",
        field="granges",
        value=value)
}

setValidity("GRangesFilter", function(object) {
    value <- .value(object)
    txt <- character()
    if (!is(value, "GRanges"))
        txt <- c(txt, "'value' must be 'GRanges' object")
    if (length(txt)) txt else TRUE
})

setValidity("BasicFilter", function(object) {
    value <- .value(object)
    condition <- .condition(object)
    isCharacter <- .isCharacter(object)
    txt <- character()
    if (length(condition) != 1L)
        txt <- c(txt, "'condition' must be length 1")
    if (!condition %in% .OPS)
        txt <- c(txt,
                 sprintf("'condition' must be one of %s",
                         paste("'", .OPS, "'", collapse=", ")))
    if (isCharacter && !is.character(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take character value"))
    if (!isCharacter && (!is.integer(value)) || is.na(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take integer value"))
    if (condition  %in% c("startsWith", "endsWith", ">", "<", ">=", "<=") &&
        length(value) > 1L)
        txt <- c(txt,
                 paste0("'value' must be length 1 when condition is '",
                        condition, "'"))
    if (condition  %in% c("startsWith", "endsWith") && !isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with character value"))
    if (condition  %in% c(">", "<", ">=", "<=") && isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with integer value"))
    if (length(txt)) txt else TRUE
})

.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CHAR_FIELDS <- c(
    "accnum", "alias", "cds_chrom", "cds_id", "cds_name", "cds_strand",
    "ensembl", "ensemblprot", "ensembltrans", "entrez", "enzyme", "evidence",
    "evidenceall", "exon_chrom", "exon_id", "exon_name", "exon_rank", 
    "exon_strand", "flybase", "flybase_cg", "flybase_prot", "gene_chrom",
    "gene_strand", "genename", "go", "goall", "ipi", "map", "mgi", "omim",
    "ontology", "ontologyall", "pfam", "pmid", "prosite", "refseq", "symbol",
    "tx_chrom", "tx_id", "tx_name", "tx_strand", "tx_type", "unigene",
    "uniprot", "wormbase", "zfin")

.INT_FIELDS <- c(
    "cds_start", "cds_end", "exon_start", "exon_end", "gene_start", "gene_end",
    "tx_start", "tx_end")

.filterFactory <- function(field, isCharacter) {
    class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
    .valueIsCharacter <- isCharacter
    .as.value <-
        if (.valueIsCharacter) {
            as.character
        } else {
            function(x) {
                stopifnot(is.numeric(x))
                as.integer(x)
            }
        }
    function(value, condition="==") {
        new(class,
            field=field,
            condition=condition,
            value=.as.value(value),
            .valueIsCharacter=.valueIsCharacter)
    }
}

## create filter functions
local({
    for (field in c(.CHAR_FIELDS, .INT_FIELDS))  {
        class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
        setClass(class, contains="BasicFilter", where=topenv())
        if (field %in% .CHAR_FIELDS)
            assign(class, .filterFactory(field, TRUE), envir=topenv())
        else
            assign(class, .filterFactory(field, FALSE), envir=topenv())
    }
})

.field <- function(x) x@field

.condition <- function(x) x@condition

.value <- function(x) x@value

.isCharacter <- function(x) x@.valueIsCharacter


#' @param object A \code{src_organism} object
#'
#' @rdname BasicFilter
#' @export
setMethod("show", "BasicFilter", function(object){
    cat("class:", class(object), "\n")
    cat("condition:", .condition(object), "\n")
    cat("value:", .value(object), "\n")
})

.fields <- function(filter) {
    vapply(filter, .field, character(1))
}

.convertFilter <- function(filter) {
    field <- .field(filter)
    value <- .value(filter)
    condition <- .condition(filter)

    op <- switch(
        condition,
        "==" = if (length(value) == 1) "==" else "%in%",
        "!=" = if (length(value) == 1) "!=" else "%in%",
        "startsWith" = "%like%",
        "endsWith" = "%like%"
    )

    if (condition %in% c("==", "!="))
        value <- paste0("'", value, "'", collapse=", ")

    if (!is.null(op) && op %in% c("==", "!="))
        sprintf("%s %s %s", field, op, value)
    else if ((condition == "==") && op == "%in%")
        sprintf("%s %s c(%s)", field, op, value)
    else if ((condition == "!=") && op == "%in%")
        sprintf("!%s %s c(%s)", field, op, value)
    else if (condition == "startsWith")
        sprintf("%s %s '%s%%'", field, op, value)
    else if (condition == "endsWith")
        sprintf("%s %s '%%%s'", field, op, value)
    else if (condition %in% c(">", "<", ">=", "<=")) {
        sprintf("%s %s %s", field, condition, as.integer(value))
    }
}

#' @rdname BasicFilter
#' @export
supportedFilters <- function() {
    paste0(sub("([a-z])", "\\U\\1", c(.CHAR_FIELDS, .INT_FIELDS),
               perl=TRUE), "Filter")
}
