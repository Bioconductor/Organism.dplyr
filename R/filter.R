#' Filtering src_organism objects
#'
#' These functions create filters to be used by the "select"
#' interface to \code{src_organism} objects.
#'
#' All filters except \code{GRangesFilter()} takes value(s) from
#' corresponding fields in the data base. For example,
#' \code{AccnumFilter()} takes values of accession number(s), which
#' come from field \code{accnum}. See \code{keytypes()} and
#' \code{keys()} for possible values.
#'
#' \code{GRangesFilter()} takes a \code{GRanges} object as filter, and returns
#' genomic extractors (\code{genes}, \code{transcripts}, etc.) that are
#' partially overlapping with the region.
#'
#' \code{supportedFilters()} lists all available filters for
#' \code{src_organism} object.
#'
#' @aliases AccnumFilter AliasFilter CdsChromFilter CdsIdFilter
#'     CdsNameFilter CdsStrandFilter EnsemblFilter EnsemblprotFilter
#'     EnsembltransFilter EntrezFilter EnzymeFilter EvidenceFilter
#'     EvidenceallFilter ExonChromFilter ExonIdFilter ExonNameFilter
#'     ExonRankFilter ExonStrandFilter FlybaseFilter FlybaseCgFilter
#'     FlybaseProtFilter GeneChromFilter GeneStrandFilter GenenameFilter
#'     GoFilter GoallFilter IpiFilter MapFilter MgiFilter OmimFilter
#'     OntologyFilter OntologyallFilter PfamFilter PmidFilter PrositeFilter
#'     RefseqFilter SymbolFilter TxChromFilter TxIdFilter TxNameFilter
#'     TxStrandFilter TxTypeFilter UnigeneFilter UniprotFilter WormbaseFilter
#'     ZfinFilter CdsStartFilter CdsEndFilter ExonStartFilter ExonEndFilter
#'     GeneStartFilter GeneEndFilter TxStartFilter TxEndFilter
#'
#' @usage AccnumFilter(value, condition = "==")
#' AliasFilter(value, condition = "==")
#' CdsChromFilter(value, condition = "==")
#' CdsIdFilter(value, condition = "==")
#' CdsNameFilter(value, condition = "==")
#' CdsStrandFilter(value, condition = "==")
#' EnsemblFilter(value, condition = "==")
#' EnsemblprotFilter(value, condition = "==")
#' EnsembltransFilter(value, condition = "==")
#' EntrezFilter(value, condition = "==")
#' EnzymeFilter(value, condition = "==")
#' EvidenceFilter(value, condition = "==")
#' EvidenceallFilter(value, condition = "==")
#' ExonChromFilter(value, condition = "==")
#' ExonIdFilter(value, condition = "==")
#' ExonNameFilter(value, condition = "==")
#' ExonRankFilter(value, condition = "==")
#' ExonStrandFilter(value, condition = "==")
#' FlybaseFilter(value, condition = "==")
#' FlybaseCgFilter(value, condition = "==")
#' FlybaseProtFilter(value, condition = "==")
#' GeneChromFilter(value, condition = "==")
#' GeneStrandFilter(value, condition = "==")
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
#' TxChromFilter(value, condition = "==")
#' TxIdFilter(value, condition = "==")
#' TxNameFilter(value, condition = "==")
#' TxStrandFilter(value, condition = "==")
#' TxTypeFilter(value, condition = "==")
#' UnigeneFilter(value, condition = "==")
#' UniprotFilter(value, condition = "==")
#' WormbaseFilter(value, condition = "==")
#' ZfinFilter(value, condition = "==")
#' CdsStartFilter(value, condition = "==")
#' CdsEndFilter(value, condition = "==")
#' ExonStartFilter(value, condition = "==")
#' ExonEndFilter(value, condition = "==")
#' GeneStartFilter(value, condition = "==")
#' GeneEndFilter(value, condition = "==")
#' TxStartFilter(value, condition = "==")
#' TxEndFilter(value, condition = "==")
#'
#' @param value Value of the filter. For \code{GRangesFilter} value should be a
#'     \code{GRanges} object.
#'
#' @param condition The condition to be used in filter for genomic
#'     extractors, one of "==", "!=", "startsWith", "endsWith", ">",
#'     "<", ">=", "<=".  For character values "==", "!=", "startsWith"
#'     and "endsWith" are allowed, for numeric values
#'     (\code{CdsStartFilter}, \code{CdsEndFilter},
#'     \code{ExonStartFilter}, \code{ExonEndFilter},
#'     \code{GeneStartFilter}, \code{GeneEndFilter},
#'     \code{TxStartFilter} and \code{TxEndFilter}), "==", "!=", ">",
#'     ">=", "<" and "<=". Default condition is "==".
#'
#' @return A Filter object showing class, value and condition of the filter
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
#' @author Yubo Cheng.
#'
#' @examples
#' src <- src_organism(dbpath=hg38light())
#' keytypes(src)
#' head(keys(src, "ensembl"))
#'
#' ## filter by ensembl
#' EnsemblFilter("ENSG00000171862")
#'
#' ## filter by gene symbol start with "BRAC"
#' SymbolFilter("BRCA", "startsWith")
#'
#' ## filter by GRanges
#' GRangesFilter(GenomicRanges::GRanges("chr10:87869000-87876000"))
#'
#' ## filter by transcript start position
#' TxStartFilter(87863438, ">")
#'
#' @export AccnumFilter AliasFilter CdsChromFilter CdsIdFilter CdsNameFilter
#' @export CdsStrandFilter EnsemblFilter EnsemblprotFilter EnsembltransFilter
#' @export EntrezFilter EnzymeFilter EvidenceFilter EvidenceallFilter
#' @export ExonChromFilter ExonIdFilter ExonNameFilter ExonRankFilter
#' @export ExonStrandFilter FlybaseFilter FlybaseCgFilter FlybaseProtFilter
#' @export GeneChromFilter GeneStrandFilter GenenameFilter GoFilter
#' @export GoallFilter IpiFilter MapFilter MgiFilter OmimFilter OntologyFilter
#' @export OntologyallFilter PfamFilter PmidFilter PrositeFilter RefseqFilter
#' @export SymbolFilter TxChromFilter TxIdFilter TxNameFilter TxStrandFilter
#' @export TxTypeFilter UnigeneFilter UniprotFilter WormbaseFilter ZfinFilter
#' @export CdsStartFilter CdsEndFilter ExonStartFilter ExonEndFilter
#' @export GeneStartFilter GeneEndFilter TxStartFilter TxEndFilter
#' @rdname filter
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

#' @rdname filter
#' @export
setClass("GRangesFilter",
         slots = list(field="character",
                      value="GRanges"))

#' @rdname filter
#' @export
GRangesFilter <- function(value) {
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

#' @rdname filter
#' @exportMethod show
setMethod("show", "GRangesFilter",
    function(object)
{
    cat("class:", class(object),
        "\nvalue:\n")
    print(.value(object))
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
    "accnum", "alias", "cds_chrom", "cds_name", "cds_strand",
    "ensembl", "ensemblprot", "ensembltrans", "entrez", "enzyme",
    "evidence", "evidenceall", "exon_chrom", "exon_name",
    "exon_strand", "flybase", "flybase_cg", "flybase_prot",
    "gene_chrom", "gene_strand", "genename", "go", "goall", "ipi",
    "map", "mgi", "omim", "ontology", "ontologyall", "pfam", "pmid",
    "prosite", "refseq", "symbol", "tx_chrom", "tx_name", "tx_strand",
    "tx_type", "unigene", "uniprot", "wormbase", "zfin")

.INT_FIELDS <- c(
    "cds_id", "cds_start", "cds_end", "exon_id", "exon_start",
    "exon_end", "exon_rank", "gene_start", "gene_end", "tx_id",
    "tx_start", "tx_end")

.fieldToClass <- function(field) {
    class <- sub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, "Filter")
}

.filterFactory <- function(field, class, .valueIsCharacter) {
    force(field); force(class)          # watch for lazy evaluation
    as.value <-
        if (.valueIsCharacter) {
            as.character
        } else {
            function(x) {
                stopifnot(is.numeric(x))
                as.integer(x)
            }
        }
    function(value, condition = "==") {
        value <- as.value(value)
        new(class, field=field, condition=condition,
            value=value, .valueIsCharacter=.valueIsCharacter)
    }
}

## create filter functions
local({
    field <- c(.CHAR_FIELDS, .INT_FIELDS)
    class <- .fieldToClass(field)
    for (i in seq_along(field)) {
        setClass(class[[i]], contains="BasicFilter", where=topenv())
        assign(class[[i]],
               .filterFactory(
                   field[[i]], class[[i]], field[[i]] %in% .CHAR_FIELDS
               ),
               envir=topenv())
    }
})

.field <- function(x) x@field

.condition <- function(x) x@condition

.value <- function(x) x@value

.isCharacter <- function(x) x@.valueIsCharacter

#' @param object A \code{BasicFilter} or \code{GRangesFilter} object
#'
#' @importFrom methods show
#' @rdname filter
#' @exportMethod show
setMethod("show", "BasicFilter",
    function(object)
{
    cat("class:", class(object),
        "\ncondition:", .condition(object),
        "\nvalue:", .value(object), "\n")
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

#' @rdname filter
#' @export
supportedFilters <- function() {
    .fieldToClass(c(.CHAR_FIELDS, .INT_FIELDS))
}
