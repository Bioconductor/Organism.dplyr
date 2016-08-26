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
#' @param org character(1) naming a \code{org.*} package (e.g.,
#'     \code{org.Hs.eg.db}) or \code{OrgDb} object instantiating the
#'     content of an \code{org.*} pacakge.
#'
#' @return A dplyr \code{tbl_sql} instance representing the data
#'     table.
#'
#' @examples
#' org = "org.Hs.eg.db"
#' idmap = tbl_org_idmap(org)
#' idmap
#' idmap %>% dplyr::select(ensembl, symbol) %>% filter(symbol == "PTEN")
#' symbols <- paste0("BRCA", 1:2)
#' idmap %>% dplyr::select(ensembl, symbol) %>% filter(symbol %in% symbols)
#' 
#' ipi = tbl_org_ipi(org)
#' inner_join(idmap, ipi) %>% filter(map == "17q21") %>%
#'      dplyr::select(entrez, symbol, ensembl, ipi)
#'
#' @importFrom RSQLite dbSendQuery dbListTables
#' @importFrom AnnotationDbi dbconn dbfile
#' @importFrom dplyr src_sql "%>%" tbl select_
#'
#' @rdname tbl_org
#' @export
#' @export
select_.tbl_organism <- function(.data, ...) { 
    .data = NextMethod(.data, ...) 
    dplyr::distinct_(.data, ...) 
}

#' Create a dplyr view of genomic coordinates
#'
#' @param txdb character(1) naming a \code{TxDb.*} package (e.g.,
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#' 
#' @examples
#' txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene"
#' tbl_txdb(txdb) %>% filter(entrez == "1") %>% 
#' dplyr::select(txid, txname, txstart, txend)
#' 
#' org = "org.Hs.eg.db"
#' idmap = tbl_org_idmap(org)
#' transcript = tbl_txdb(txdb)
#' inner_join(transcript, idmap, copy = TRUE) %>% filter(symbol == "PTEN") %>% 
#' dplyr::select(entrez, symbol, txid, txname, txstart, txend)
#' 
#' entrez <- idmap %>% filter(symbol == "PTEN") %>% dplyr::select(entrez, symbol)
#' inner_join(transcript, entrez, copy = TRUE) %>% filter(symbol == "PTEN") %>% 
#' dplyr::select(entrez, symbol, txid, txname, txstart, txend)
#' 
#' @export
tbl_txdb <- function(txdb) {
    .get_tbl(txdb, "txdb")
}

#' Create a dplyr view integrating org.* and TxDb.* information
#'
#' @inheritParams tbl_txdb
#' 
#' @export
src_organism <- function(org, txdb) {
    db <- org
    db <- loadNamespace(db)[[db]]
    conn = dbconn(db)

    fname <- system.file(package="Organism.dplyr", "schema", "organism.sql")
    schemas <- readLines(fname)
    grps <- cumsum(!nzchar(schemas)) + 1
    for (schema in split(schemas, grps)) {
        sql <- paste(schema, collapse="\n")
        dbSendQuery(conn, sql)
    }

    src <- src_sql("sqlite", conn, path=dbfile(db))
    class(src) = c("src_organism", class(src))
    src
}

#' @export
src_tbls.src_organism <- function(x) {
    sql <- "SELECT name FROM sqlite_temp_master WHERE type = 'view';"
    dbGetQuery(xx$con, sql)$name
}
