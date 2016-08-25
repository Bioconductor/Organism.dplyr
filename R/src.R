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
tbl_org_idmap <- function(org) 
{
    .get_tbl(org, "idmap")
}
    
#' @rdname tbl_org
#' @export
tbl_org_genename <- function(org)
{
    .get_tbl(org, "genename")
}

#' @rdname tbl_org
#' @export
tbl_org_pmid <- function(org)
{
    .get_tbl(org, "pmid")
}

#' @rdname tbl_org
#' @export
tbl_org_ensenbltrans <- function(org)
{
    .get_tbl(org, "ensenbltrans")
}

#' @rdname tbl_org
#' @export
tbl_org_protein <- function(org)
{
    .get_tbl(org, "protein")
}

#' @rdname tbl_org
#' @export
tbl_org_go <- function(org)
{
    .get_tbl(org, "view_go")
}

#' @rdname tbl_org
#' @export
tbl_org_go_all <- function(org)
{
    .get_tbl(org, "view_go_all")
}

#' @rdname tbl_org
#' @export
tbl_org_alias <- function(org)
{
    .get_tbl(org, "view_alias")
}

#' @rdname tbl_org
#' @export
tbl_org_ipi <- function(org)
{
    .get_tbl(org, "ipi")
}

.get_tbl <- function(db, tblname) {
    stopifnot(is.character(db), length(db) == 1L)
    stopifnot(is.character(tblnames), length(tblname) == 1)

    if (substr(db, 1, 3) == "org")
        fname <- paste0("org_", tblname, ".sql")
    else
        fname <- paste0(tblname, ".sql")
    db <- loadNamespace(db)[[db]]
    conn = dbconn(db)

    schema <- system.file(package="Organism.dplyr", "schema", fname)
    sql <- paste(readLines(schema), collapse="\n")
    dbSendQuery(conn, sql)

    tbl = src_sql("sqlite", conn, path=dbfile(db)) %>% tbl(tblname)
    class(tbl) = c("tbl_org", class(tbl))
    tbl
}

#' @rdname tbl_org
#' @export
tbl_go <- function(org) {}

#' @export
select_.tbl_org <- function(.data, ...) { 
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
    tbl = .get_tbl(txdb, "txdb")
    class(tbl) = c("tbl_org", class(tbl))
    tbl
}

#' Create a dplyr view integrating org.* and TxDb.* information
#'
#' @inheritParams tbl_org_idmap
#' @inheritParams tbl_txdb
#' 
#' @export
src_organism <- function(org, txdb) {}
