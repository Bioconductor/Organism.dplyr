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
#'     CharacterFilter-class IntegerFilter-class show,CharacterFilter-method
#'     show,IntegerFilter-method
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
#' @export EnzymeFilter EvidenceFilter EvidenceallFilter
#' @export ExonChromFilter
#' @export ExonStrandFilter FlybaseFilter FlybaseCgFilter FlybaseProtFilter
#' @export GeneChromFilter GeneStrandFilter GoFilter
#' @export GoallFilter IpiFilter MapFilter MgiFilter OmimFilter OntologyFilter
#' @export OntologyallFilter PfamFilter PmidFilter PrositeFilter RefseqFilter
#' @export TxChromFilter TxStrandFilter
#' @export TxTypeFilter UnigeneFilter WormbaseFilter ZfinFilter
#' @rdname filter
#' @importFrom methods new setClass slot setMethod setValidity
#' @importFrom AnnotationFilter AnnotationFilter GRangesFilter field value
#'      condition
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

setMethod("initialize", "BasicFilter", function(.Object) {
    .Deprecated("AnnotationFilter")
})

setValidity("BasicFilter", function(object) {
    value <- value(object)
    condition <- condition(object)
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
                 paste0("'", class(object), "' can only take character value"))
    if (!isCharacter && (!is.integer(value)) || is.na(value))
        txt <- c(txt,
                  paste0("'", class(object), "' can only take integer value"))
    if (condition  %in% c("startsWith", "endsWith", ">", "<", ">=", "<=") &&
        length(value) > 1L)
        txt <- c(txt,
                 paste0("'value' must be length 1 when condition is '",
                         condition, "'"))
    if (condition  %in% c("startsWith", "endsWith") && !isCharacter)
        txt <- c(txt,
                 paste0("'", condition, "' can only work with character value"))
    if (condition  %in% c(">", "<", ">=", "<=") && isCharacter)
        txt <- c(txt,
                 paste0("'", condition, "' can only work with integer value"))
    if (length(txt)) txt else TRUE
})

.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CONDITION <- list(
    IntegerFilter = c("==", "!=", ">", "<", ">=", "<="),
    CharacterFilter = c("==", "!=", "startsWith", "endsWith", "contains"),
    GRangesFilter = c("any", "start", "end", "within", "equal")
)

.FIELD <- list(
    CharacterFilter = c(
        "accnum", "alias", "cds_chrom", "cds_name", "cds_strand",
        "ensembl", "ensemblprot", "ensembltrans", "entrez", "enzyme",
        "evidence", "evidenceall", "exon_chrom", "exon_name",
        "exon_strand", "flybase", "flybase_cg", "flybase_prot",
        "gene_chrom", "gene_strand", "genename", "go", "goall", "ipi",
        "map", "mgi", "omim", "ontology", "ontologyall", "pfam", "pmid",
        "prosite", "refseq", "symbol", "tx_chrom", "tx_name", "tx_strand",
        "tx_type", "unigene", "uniprot", "wormbase", "zfin"),
    IntegerFilter = c(
        "cds_id", "cds_start", "cds_end", "exon_id", "exon_start",
        "exon_end", "exon_rank", "gene_start", "gene_end", "tx_id",
        "tx_start", "tx_end")
)

.valid_condition <- function(condition, class){
    txt <- character()

    test0 <- length(condition) == 1L
    if (!test0)
        txt <- c(txt, "'condition' must be length 1")

    test1 <- test0 && (condition %in% .CONDITION[[class]])
    if (!test1) {
        value <- paste(sQuote(.CONDITION[[class]]), collapse=" ")
        txt <- c(txt, paste0("'", condition, "' must be in ", value))
    }

    if(length(txt)) txt else TRUE
}

.fieldToClass <- function(field) {
    class <- sub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, "Filter")
}

.filterFactory <- function(field, class) {
    force(field); force(class)          # watch for lazy evaluation
    as.value <-
        if (field %in% .FIELD[["CharacterFilter"]]) {
            as.character
        } else {
            function(x) {
                stopifnot(is.numeric(x))
                as.integer(x)
            }
        }
    function(value, condition = "==") {
        value <- as.value(value)
        condition <- as.character(condition)
        new(class, field=field, condition=condition, value=value)
    }
}

## create filter functions not already implemented in AnnotationFilter
.filter_init <- function() {
    makeClass <- function(contains){
        fields <- .FIELD[[contains]]
        supported <- as.character(supportedFilters()[,2])
        fields <- fields[!(fields %in% supported)]
        classes <- .fieldToClass(fields)
        for (i in seq_along(fields)) {
            setClass(classes[[i]], contains=contains, where=topenv())
            assign(
               classes[[i]],
               .filterFactory(fields[[i]], classes[[i]]),
               envir=topenv()
            )
        }
    }
    for (contains in names(.FIELD))
        makeClass(contains)
}

#' @param object A \code{BasicFilter} or \code{GRangesFilter} object
#'
#' @importFrom methods show
#' @rdname filter
#' @exportMethod show
setMethod("show", "BasicFilter",
    function(object)
{
    cat("class:", class(object),
        "\ncondition:", condition(object),
        "\nvalue:", value(object), "\n")
})

.fields <- function(object) {
    res <- lapply(object, function(x) {
            if(is(x, "AnnotationFilter"))
                field(x)
            else
                .fields(x)
        })
    unlist(res)
}

.supportedFilters <- function() {
    df <- data.frame(
        filter = c(.fieldToClass(unlist(.FIELD, use.names=FALSE)),
            "GRangesFilter"),
        field = c(unlist(.FIELD, use.names=FALSE), "granges")
    )
    df[order(df[,1]),]
}

#' @rdname filter
#' @export
setMethod("supportedFilters", "src_organism", function(object){
    .supportedFilters()
})
