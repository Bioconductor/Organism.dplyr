.getFields <-
    function(x)
{
    fields <- lapply(src_tbls(x), function(table) colnames(tbl(x, table)))
    unique(unlist(fields, use.names=FALSE))
}

#' Using the "select" interface on src_organism objects
#' 
#' \code{select}, \code{columns} and \code{keys} can be used together to 
#' extract data from a \code{\link{src_organism}} object.
#' 
#' \code{keytypes()}: discover which keytypes can be passed to keytype 
#'     argument of methods \code{select} or \code{keys}.
#' 
#' \code{keys()}: returns keys for the \code{src_organism} object. By default 
#'     it returns the primary keys for the database, and returns the keys from 
#'     that keytype when the keytype argument is used.
#' 
#' \code{columns()}: discover which kinds of data can be returned for the 
#'     \code{src_organism} object.
#' 
#' \code{select()}: retrieves the data as a \code{tibble} based on parameters 
#'     for selected keys columns and keytype arguments. If requested columns 
#'     that have multiple matches for the keys, `select()` will return a 
#'     \code{tibble} with one row for each possible match. 
#' 
#' \code{mapIds()}: gets the mapped ids (column) for a set of keys that are of 
#'     a particular keytype. Usually returned as a named character vector.
#' 
#' @param x a \code{src_organism} object
#' 
#' @param keys the keys to select records for from the database. All possible 
#'     keys are returned by using the \code{keys} method.
#' 
#' @param columns the columns or kinds of things that can be retrieved
#'     from the database. As with keys, all possible columns are
#'     returned by using the \code{columns} method.
#'     
#' @return \code{keys}, \code{columns} and \code{keytypes} each
#'     returns a character vector of possible values. \code{select}
#'     returns a \code{tibble}.
#' 
#' @seealso \code{\link{AnnotationDb-class}} for more descriptsion of
#'     methods \code{select}, \code{keytypes}, \code{keys} and
#'     \code{columns}.
#'     
#'     \code{\link{src_organism}} for creating a \code{src_organism}
#'     object.
#'     
#'     \code{\link[Organism.dplyr]{transcripts_tbl}} for generic
#'     functions to extract genomic features from a
#'     \code{src_organism} object.
#'     
#' @importFrom AnnotationDbi keytypes
#' @rdname select-src_organism-method
#' 
#' @examples
#' src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' 
#' ## keytypes
#' keytypes(src)
#' 
#' @export
setMethod("keytypes", "src_organism", .getFields)

#' @importFrom AnnotationDbi columns
#' @rdname select-src_organism-method
#' 
#' @examples 
#' ## columns
#' columns(src)
#' 
#' @export
setMethod("columns", "src_organism", .getFields)

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
    res <- table %>% select_(keytype) %>% collect(n=Inf)
    
    if (length(res) == 0L)
        stop("'", keytype, "' is not a supported keytype; see 'keytypes()'")
    
    as.character(res[[keytype]])
}

#' @param keytype specifies the kind of keys that will be returned. By
#'     default keys will return the keys for schema of the
#'     \code{src_organism} object.
#'     
#' @param ... other arguments. These include: 
#' 
#'     pattern: the pattern to match.
#'     
#'     column: the column to search on.
#'     
#'     fuzzy: TRUE or FALSE value. Use fuzzy matching? (this is used
#'            with pattern)
#' 
#' @importFrom AnnotationDbi keys testForValidKeytype
#' 
#' @examples 
#' ## keys
#' keys(src, "entrez")
#' 
#' @rdname select-src_organism-method
#' @export
setMethod("keys", "src_organism", function(x, keytype, ...) {
    if (missing(keytype)) 
        keytype <- x$schema
    AnnotationDbi:::smartKeys(x = x, keytype = keytype, ..., FUN = .keys)
})

.filterByKeys <- function(x, keys, keytype, cnames) {
    table <- tbl(x, .findTable(x, keytype))
    values <- paste0("'", keys, "'", collapse=", ")
    op <- if (length(keys) == 1) "==" else "%in%"
    filter <- sprintf("%s %s c(%s)", keytype, op, values)
    
    fields <- colnames(table) 
    keyfields <- fields[fields %in% c(x$schema, "tx_id", "exon_rank")]
    fields <- unique(c(keyfields, cnames[cnames %in% fields]))
    
    table <- table %>% filter_(filter) 
    do.call(select_, c(list(table), as.list(fields)))
}

.selectColumns <- function(x, table, keytype, cnames) {
    maintbl <- .findTable(x, keytype)
    tbls <- setdiff(src_tbls(x), maintbl)
    fields <- setdiff(cnames, keytype)
    fields <- setdiff(fields, colnames(table))
    
    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || length(keep) == 0)
            next

        if (all(c("tx_id", "exon_rank") %in% colnames(table)) && 
            all(c("tx_id", "exon_rank") %in% colnames(tbl(x, i)))) {
            table <- left_join(table, tbl(x, i), by = c("tx_id", "exon_rank"))
        } else if ("tx_id" %in% colnames(table) && 
                   "tx_id" %in% colnames(tbl(x, i))) {
            table <- left_join(table, tbl(x, i), by = "tx_id")
        } else {
            table <- left_join(table, tbl(x, i))
        }
        if ("entrez.x" %in% colnames(table))
            table <- rename(table, entrez = entrez.x)
        else if ("ensembl.x" %in% colnames(table))
            table <- rename(table, ensembl = ensembl.x)
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
        stop("no keys found for the keytype you specified.")
    }
    cnames <- unique(c(keytype, columns))
    table <- .filterByKeys(x, keys, keytype, cnames)
    .selectColumns(x, table, keytype, cnames)
}

 
#' @importFrom AnnotationDbi select testSelectArgs
#' 
#' @examples
#' keytype <- "symbol"
#' keys <- c("PTEN", "BRCA1")
#' columns <- c("entrez", "tx_id", "tx_name","exon_id")
#' 
#' ## select
#' select(src, keys, columns, keytype)
#' 
#' @export
setMethod("select", "src_organism", .select)

.mapIds_base <- 
    function(x, keys, column, keytype, 
        multiVals = c("first", "filter", "asNA", "list", "CharacterList")) 
{
    if (missing(multiVals)) 
        multiVals <- "first"
    if (!is.function(multiVals)) 
        match.arg(multiVals)
    if (length(keys) < 1) 
        stop(wmsg("mapIds must have at least one key to match against."))
    if (length(column) > 1) 
        stop(wmsg("mapIds can only use one column."))
    if (length(keytype) > 1) 
        stop(wmsg("mapIds can only use one keytype."))
    res <- select(x, 
                  keys = unique(keys), 
                  columns = column, 
                  keytype = keytype) %>% 
           collect(n=Inf) %>% as.data.frame()
    res <- split(res[[column]], f = res[[keytype]])[keys]
    .filter <- function(data) {
        idx <- elementNROWS(data) == 1
        unlist(data[idx])
    }
    .asNA <- function(data) {
        idx <- elementNROWS(data) > 1
        data[idx] <- NA
        unlist(data)
    }
    if (is.function(multiVals)) {
        sapply(res, FUN = multiVals)
    }
    else {
        switch(multiVals, 
               list = res, 
               filter = .filter(res), 
               asNA = .asNA(res), 
               CharacterList = as(res, "CharacterList"), 
               first = sapply(res, function(x) x[[1]])
               )
    }
}

#' @param column character(1) the column to search on, can only have a
#'     single element for the value
#' 
#' @param multiVals What should \code{mapIds} do when there are
#'     multiple values that could be returned. Options include:
#'     
#'     first: when there are multiple matches only the 1st thing that
#'            comes back will be returned. This is the default
#'            behavior.
#'            
#'     list: return a list object to the end user
#'     
#'     filter: remove all elements that contain multiple matches and
#'             will therefore return a shorter vector than what came
#'             in whenever some of the keys match more than one value
#'     
#'     asNA: return an NA value whenever there are multiple matches
#'     
#'     CharacterList: returns a SimpleCharacterList object
#'     
#'     FUN: can also supply a function to the multiVals argument for custom 
#'          behaviors. The function must take a single argument and return a 
#'          single value. This function will be applied to all the elements 
#'          and will serve a 'rule' that for which thing to keep when there 
#'          is more than one element. So for example this example function 
#'          will always grab the last element in each result: 
#'          \code{last <- function(x){x[[length(x)]]}}
#' 
#' @importFrom AnnotationDbi mapIds
#' @importFrom IRanges elementNROWS
#' 
#' @examples
#' ## mapIds
#' mapIds(src, keys, column = "tx_name", keytype)
#' 
#' @rdname select-src_organism-method
#' @export
setMethod("mapIds", "src_organism", 
function(x, keys, column, keytype, ..., multiVals) {
    .mapIds_base(x, keys, column, keytype, ..., multiVals = multiVals)
})
