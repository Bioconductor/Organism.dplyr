#' Create a dplyr view into an org package
#'
#' The view provides a convenient way to map between gene, transcript,
#' and protein identifiers.
#'
#' This function provides a way to map between identifiers. It does
#' not include details of GO evidence and ontology (see
#' \code{\link{tbl_go}}), and does not provide genomic coordinates
#' (see \code{\link{tbl_txdb}}.
#'
#' This function is meant to be a building block for
#' \code{\link{src_organism}}, which provides an integrated
#' presentation of identifiers and genomic coordinates.
#'
#'
#'
#' Create a dplyr view integrating org.* and TxDb.* information
#' 
#' @param org character(1) naming a \code{org.*} package (e.g.,
#'     \code{org.Hs.eg.db}) or \code{OrgDb} object instantiating the
#'     content of an \code{org.*} pacakge.
#'     
#' @param txdb character(1) naming a \code{TxDb.*} package (e.g.,
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#' 
#' @return A dplyr \code{tbl_sql} instance representing the data
#'     table.
#'     
#' @importFrom RSQLite dbSendQuery dbGetQuery
#' @importFrom AnnotationDbi dbconn dbfile
#' 
#' @examples
#' organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
#' src_tbls(organism)
#' tbl <- tbl(organism, "id")
#' tbl
#' tbl %>% dplyr::select(ensembl, symbol) %>% filter(symbol == "PTEN")
#' inner_join(tbl, tbl(organism, "txdb")) %>% filter(symbol == "PTEN") %>%
#'      dplyr::select(entrez, symbol, txid, txname, txstart, txend)
#' 
#' @export
src_organism <- function(org, txdb) {
    db <- loadNamespace(org)[[org]]
    con <- dbconn(db)
    
    .get_view(con, org, "organism")
    
    tbls <- dbGetQuery(src$con, "pragma database_list;")$name
    if (!("txdb" %in% tbls))
        .get_view(con, txdb, "txdb")
    
    src <- src_sql("sqlite", con, path=dbfile(con))
    class(src) = c("src_organism", class(src))
    src
}

.get_view <- function(con, db, tblname) {
    if (!(substr(db, 1, 3) == "org")) {
        db <- loadNamespace(db)[[db]]
        dbGetQuery(con, paste0("ATTACH '", dbfile(db), "' AS ", tblname))
    }
    
    fname <- system.file(package="Organism.dplyr", "schema",
                         paste0(tblname, ".sql"))
    schemas <- readLines(fname)
    grps <- cumsum(!nzchar(schemas)) + 1
    for (schema in split(schemas, grps)) {
        sql <- paste(schema, collapse="\n")
        dbSendQuery(con, sql)
    }
}

#' @importFrom RSQLite dbSendQuery dbListTables
#' @importFrom AnnotationDbi dbconn dbfile
#' @importFrom dplyr src_sql "%>%" tbl select_
#' 
#' @rdname src_organism
#' @export
select_.tbl_organism <- function(.data, ...) { 
    .data = NextMethod(.data, ...) 
    dplyr::distinct_(.data, ...) 
}

#' @importFrom dplyr src_tbls
#' @importFrom RSQLite dbGetQuery
#' @export
src_tbls.src_organism <- function(x) {
    sql <- "SELECT name FROM sqlite_temp_master WHERE type = 'view';"
    dbGetQuery(x$con, sql)$name
}

#' @importFrom dplyr tbl
#' @export
tbl.src_organism <- function(src, ...) {
    tbl <- NextMethod(src, ...)
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}
