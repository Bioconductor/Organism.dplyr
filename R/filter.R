#' @export AccnumFilter AliasFilter Cds_chromFilter Cds_idFilter Cds_nameFilter
#' @export Cds_strandFilter EnsemblFilter EnsemblprotFilter EnsembltransFilter
#' @export EntrezFilter EnzymeFilter EvidenceFilter EvidenceallFilter
#' @export Exon_chromFilter Exon_idFilter Exon_nameFilter Exon_rankFilter
#' @export Exon_strandFilter Gene_chromFilter Gene_strandFilter GenenameFilter
#' @export GoFilter GoallFilter IpiFilter MapFilter OmimFilter OntologyFilter
#' @export OntologyallFilter PfamFilter PmidFilter PrositeFilter RefseqFilter
#' @export SymbolFilter Tx_chromFilter Tx_idFilter Tx_nameFilter
#' @export Tx_strandFilter Tx_typeFilter UnigeneFilter UniprotFilter
#' @export Cds_startFilter Cds_endFilter Exon_startFilter Exon_endFilter
#' @export Gene_startFilter Gene_endFilter Tx_startFilter Tx_endFilter

## transcripts_tbl(src, filter=list(
##     SymbolFilter("BRCA", "startsWith"),
##     EntrezFilter(672)))
## 
## transcripts_tbl(src, filter=list(
##     SymbolFilter(c("PTEN", "BRCA1")), 
##     SymbolFilter("BRCA1","!="), 
##     Tx_startFilter(87863438,">"), 
##     Tx_endFilter(87933487, "<")))

#' @rdname BasicFilter
#' @importFrom methods new setClass slot
#' @export
setClass("BasicFilter",
         representation(
             "VIRTUAL",
             field="character",
             condition="character",
             value="character",
             .valueIsCharacter="logical"
         ),
         prototype=list(
             condition= "==",
             value=character(),
             .valueIsCharacter=TRUE
         )
)

#' Functions for filter objects
#' 
#' These functions are used to create filter objects for genomic extrators. 
#' 
#' All these filters except \code{GRangesFilter} extend \code{BasicFilter}. 
#' 
#' @aliases AccnumFilter AliasFilter Cds_chromFilter Cds_idFilter Cds_nameFilter
#'     Cds_strandFilter EnsemblFilter
#' 
#' @usage AccnumFilter(value, condition = "==")
#' AliasFilter(value, condition = "==")
#' Cds_chromFilter(value, condition = "==")
#' Cds_idFilter(value, condition = "==")
#' Cds_nameFilter(value, condition = "==")
#' Cds_strandFilter(value, condition = "==")
#' EnsemblFilter(value, condition = "==")
#' 
#' @param value value of the filter
#' 
#' @param condition one of "==", "!=", "startsWith", ">", "<", ">=", "<="
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
    if (!isCharacter && !is.integer(value)) 
        txt <- c(txt, 
                 paste0("'", class(object),
                        "' can only take integer value"))
    if (condition  %in% c(">", "<", ">=", "<=") && length(value) > 1L)
        txt <- c(txt, 
                 paste0("'value' must be length 1 when condition is '", 
                        condition, "'"))
    if (condition  %in% c(">", "<", ">=", "<=") && isCharacter)
        txt <- c(txt, 
                 paste0("'", condition,
                        "' can only work with integer value"))
    if (length(txt)) txt else TRUE
})

.OPS <- c("==", "!=", "startsWith", ">", "<", ">=", "<=")
.CHAR_FIELDS <- c("accnum", "alias", "cds_chrom", "cds_id", "cds_name", 
                  "cds_strand", "ensembl", "ensemblprot", "ensembltrans", 
                  "entrez", "enzyme", "evidence", "evidenceall", "exon_chrom", 
                  "exon_id", "exon_name", "exon_rank", "exon_strand", 
                  "gene_chrom", "gene_strand", "genename", "go", "goall", 
                  "ipi", "map", "omim", "ontology", "ontologyall", "pfam", 
                  "pmid", "prosite", "refseq", "symbol", "tx_chrom", "tx_id", 
                  "tx_name", "tx_strand", "tx_type", "unigene", "uniprot")
.INT_FIELDS <- c("cds_start", "cds_end", "exon_start", "exon_end", 
                 "gene_start", "gene_end", "tx_start", "tx_end")

.filterFactory <- function(field, isCharacter) {
    class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
    .valueIsCharacter <- isCharacter
    function(value, condition="==") {
        new(class,
            field=field,
            condition=condition,
            value=as.character(value),
            .valueIsCharacter=.valueIsCharacter)
    }
}


.onLoad <- function(libname, pkgname) {
    for (field in .CHAR_FIELDS) {
        class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
        setClass(class, contains="BasicFilter")
        assign(class, .filterFactory(field, TRUE),
               envir=topenv(parent.frame()))
    }
    for (field in .INT_FIELDS) {
        class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
        setClass(class, contains="BasicFilter")
        assign(class, .filterFactory(field, FALSE),
               envir=topenv(parent.frame()))
    }
}


.field <- function(x) x@field
.condition <- function(x) x@condition
.value <- function(x) {
    if (is(x, "BasicFilter"))
        if (.isCharacter(x)) x@value else as.integer(x@value)
    else
        x@value
}
.isCharacter <- function(x) x@.valueIsCharacter



#' @param object A \code{src_organism} object
#' @examples 
#' ## filter by symbol start with "BRCA"
#' SymbolFilter("BRCA", "startsWith")
#' 
#' @rdname transcripts_tbl
#' @export
setMethod("show", "BasicFilter", function(object){
    cat("Object of class:", class(object), "\n")
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
        "startsWith" = "%like%"
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
        sprintf("%s %s %s", field, op, paste0("'", value, "%", "'"))
    else if (condition %in% c(">", "<", ">=", "<=")) {
        sprintf("%s %s %s", field, condition, as.integer(value))
    }
}

 