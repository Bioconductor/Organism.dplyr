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
#' tbl = tbl_org("org.Hs.eg.db")
#' tbl
#' tbl %>% select(ensembl, symbol) %>% filter(symbol == "PTEN")
#' symbols <- paste0("BRCA", 1:2)
#' tbl %>% select(ensembl, symbol) %>% filter(symbol %in% symbols)
#'
#' @importFrom RSQLite dbSendQuery dbListTables
#' @importFrom AnnotationDbi dbconn dbfile
#' @importFrom dplyr src_sql "%>%" tbl
#' 
#' @export
tbl_org_idmap <- function(org)
    .get_tbl(org, "idmap")

#' @export
tbl_org_genename <- function(org)
    .get_tbl(org, "genename")
    
.get_tbl <- function(org, tblname) {
    if (is.character(org))
        org <- loadNamespace(org)[[org]]
    conn = dbconn(org)
    if (!tblname %in% dbListTables(conn)) {
        fname <- paste0("org_", tblname, ".sql")
        schema <- system.file(package="Organism.dplyr", "schema", fname)
        sql <- paste(readLines(schema), collapse="\n")
        dbSendQuery(conn, sql)
    }
    src_sql("sqlite", conn, path=dbfile(org)) %>% tbl(tblname)
}

#' Create a dplyr view of GO details
#'
#' @inheritParams tbl_org
#' 
#' @export
tbl_go <- function(org) {}

#' Create a dplyr view of genomic coordinates
#'
#' @param txdb character(1) naming a \code{TxDb.*} package (e.g.,
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#' 
#' @export
tbl_txdb <- function(txdb) {}

#' Create a dplyr view integrating org.* and TxDb.* information
#'
#' @inheritParams tbl_org
#' @inheritParams tbl_txdb
#' 
#' @export
src_organism <- function(org, txdb) {}
