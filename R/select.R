.getFields <- function(src) {
    fields <- lapply(src_tbls(src), function(table) colnames(tbl(src, table)))
    unique(unlist(fields))
}

#' @importFrom AnnotationDbi keytypes
#' @rdname select-src_organism-method
#' @export
setMethod("keytypes", "src_organism", function(x) {
    .getFields(x)
})


#' @importFrom AnnotationDbi columns
#' @rdname select-src_organism-method
#' @export
setMethod("columns", "src_organism", function(x) {
    .getFields(x)
})

.findTable <- function(x, field) {
    if (field == x$schema) {
        table <- "id"
    } else if (field == "tx_id") {
        table <- "ranges_tx"
    } else {
        tbls <- src_tbls(x)
        for (i in tbls) {
            if (field %in% colnames(tbl(x, i))) {
                table <- i
                break()
            }
        }
    }
    table
}

.keys <- function (x, keytype) {
    testForValidKeytype(x, keytype)
    
    table <- tbl(x, .findTable(x, keytype))
    res <- table %>% select_(keytype) %>% collect()
    
    if (length(res) == 0) {
        stop(paste(keytype, 
                   "is not a supported keytype.", " Please use the keytypes", 
                   "method to identify viable keytypes"))
    }
    
    res[[keytype]]
}

#' @param keytype specifies the kind of keys that will be returned. By default 
#'     keys will return the keys for schema of the \code{\link{src_organism}} 
#'     object.
#' 
#' @importFrom AnnotationDbi keys testForValidKeytype
#' @rdname select-src_organism-method
#' @export
setMethod("keys", "src_organism", function(x, keytype, ...) {
    if (missing(keytype)) 
        keytype <- x$schema
    AnnotationDbi:::smartKeys(x = x, keytype = keytype, ..., 
                              FUN = .keys)
})


.filterByKeys <- function(x, keys, keytype) {
    table <- tbl(x, .findTable(x, keytype))
    filter <- .tbl_filter(keytype, keys)
    
    fields <- colnames(table) 
    fields <- fields[fields %in% c(x$schema, "tx_id", "exon_rank")]
    fields <- c(fields, keytype)
    
    table <- table %>% filter_(filter) 
    do.call(select_, c(list(table), as.list(fields)))
}

.selectColumns <- function(x, table, keytype, cnames) {
    maintbl <- .findTable(x, keytype)
    tbls <- setdiff(src_tbls(x), maintbl)
    fields <- setdiff(cnames, keytype)
    
    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || length(keep) == 0)
            next
        table <- inner_join(table, tbl(x, i))
        fields <- setdiff(fields, keep)
    }
    do.call(select_, c(list(table), as.list(cnames)))
}

.select <- function (x, keys, columns, keytype) {
    if (missing(keys)) {
        stop("'keys' must be a character vector")
    }
    if (missing(columns)) {
        stop("'columns' must be a character vector")
    }
    testSelectArgs(x, keys = keys, cols = columns, keytype = keytype, 
                   skipValidKeysTest = FALSE)
    if (is.na(keys(x, keytype)[1]) & length(keys(x, keytype)) == 
        1) {
        stop(paste("There do not appear to be any keys", 
                   "for the keytype you have specified."))
    }
    cnames <- unique(c(keytype, columns))
    table <- .filterByKeys(x, keys, keytype)
    .selectColumns(x, table, keytype, cnames)
}

#' Using the "select" interface on src_organism objects
#' 
#' select, columns and keys can be used together to extract data for a 
#' \code{\link{src_organism}} object.
#' 
#' @param x a src_organism object
#' 
#' @param keys the keys to select records for from the database. All possible 
#'     keys are returned by using the \code{keys} method.
#' 
#' @param columns the columns or kinds of things that can be retrieved from  
#'     the database. As with keys, all possible columns are returned by using  
#'     the \code{columns} method.
#' 
#' @seealso \code{\link{AnnotationDb-class}} for more descriptsion of methods 
#'     \code{select}, \code{keytypes}, \code{keys} and \code{columns}.
#'     
#'     \code{\link{src_organism}} for creating a \code{src_organism} 
#'     object.
#'     
#'     \code{\link[Organism.dplyr]{transcripts_tbl}} for generic functions to 
#'     extract genomic features from a \code{src_organism} object.
#' 
#' @importFrom AnnotationDbi select testSelectArgs
#' 
#' @examples
#' src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' keytype <- "symbol"
#' keys <- c("PTEN", "BRCA1")
#' columns <- c("entrez", "tx_id", "tx_name","exon_id")
#' select(src, keys, columns, keytype)
#' 
#' @export
setMethod("select", "src_organism", function (x, keys, columns, keytype) 
{
    .select(x, keys, columns, keytype)
})
