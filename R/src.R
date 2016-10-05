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
#' @importFrom AnnotationDbi dbconn dbfile taxonomyId
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' 
#' @examples
#' # human
#' library(org.Hs.eg.db)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' organism <- src_organism(org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene)
#' src_tbls(organism)
#' id <- tbl(organism, "id")
#' id
#' id %>% dplyr::select(ensembl, symbol) %>% filter(symbol == "PTEN")
#' inner_join(id, tbl(organism, "ranges_tx")) %>% filter(symbol == "PTEN") %>%
#'      dplyr::select(entrez, symbol, tx_id, tx_start, tx_end)
# 
# # mouse
# organism <- src_organism("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene")
# inner_join(tbl(organism, "id"), tbl(organism, "ranges_gene")) %>% 
#      filter(entrez == "11287") %>%
#      dplyr::select(entrez, symbol, tx_chrom, gene_start, gene_end)
# 
# # rat
# organism <- src_organism("org.Rn.eg.db", "TxDb.Rnorvegicus.UCSC.rn4.ensGene")
# id <- tbl(organism, "id")
# ranges_gene <- tbl(organism, "ranges_gene")
# inner_join(id, ranges_gene) %>% 
#     filter(ensembl == "ENSRNOG00000028896") %>%
#     dplyr::select(ensembl, symbol, chrom, start, end)
# 
# # transcripts()
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# transcripts(txdb, filter=list(tx_name="uc001aaa.3"))
# 
# human <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
# 
# tx <- inner_join(tbl(human, "ranges_tx"), tbl(human, "id_ranges")) %>%
#     dplyr::select(chrom, start, end, strand, txid, txname)
# 
# tx %>% filter(txname == "uc001aaa.3")
#' 
#' 
#' @export
src_organism <- function(org, txdb) {
    if (is.character(org))
        org <- loadNamespace(org)[[org]]
    stopifnot(is(org, "OrgDb"))

    if (is.character(txdb))
        txdb <- loadNamespace(txdb)[[txdb]]
    stopifnot(is(txdb, "TxDb"))

    stopifnot(identical(taxonomyId(org), taxonomyId(txdb)))

    org_meta <- metadata(org)
    org_schema <- org_meta$value[org_meta$name == "DBSCHEMA"]
    txdb_meta <- metadata(txdb)
    txdb_type <- txdb_meta$value[txdb_meta$name == "Type of Gene ID"]
    if (startsWith(tolower(txdb_type), "entrez"))
        txdb_schema <- "txdb_entrez"
    else if (startsWith(tolower(txdb_type), "ensembl"))
        txdb_schema <- "txdb_ensembl"
    else
        stop("unknown TxDb type: ", sQuote(txdb_type))

    .add_view(dbconn(org), org, org_schema)
    .add_view(dbconn(org), txdb, txdb_schema)
    
    src <- src_sql("sqlite", dbconn(org), path=dbfile(org))
    class(src) = c("src_organism", class(src))
    src
}

.add_view <- function(con, db, tblname) {
    tbls <- dbGetQuery(con, "pragma database_list;")$name
    if (tblname %in% tbls)
        return()
    if (is(db, "TxDb"))
        dbGetQuery(con, paste0("ATTACH '", dbfile(db), "' AS ", tblname))
    
    fname <- system.file(
        package="Organism.dplyr", "schema", paste0(tblname, ".sql"))
    if (!file.exists(fname))
        stop(sQuote(tblname), " schema not unsupported")

    schemas <- readLines(fname)
    grps <- cumsum(!nzchar(schemas)) + 1
    for (schema in split(schemas, grps)) {
        sql <- paste(schema, collapse="\n")
        dbSendQuery(con, sql)
    }
    if (is(db, "TxDb"))
        dbGetQuery(con, paste0("DETACH '", tblname, "'"))
}

#' @rdname src_organism
#' @importFrom GenomeInfoDb genomeBuilds
#' @importFrom utils read.csv tail
#' @export
src_ucsc <- function(organism, genome = NULL, id = NULL, verbose=TRUE) {
    stopifnot(is.character(organism), length(organism) == 1L)
    if (!missing(genome))
        stopifnot(is.character(genome), length(genome) == 1L)
    if (!missing(id))
        stopifnot(is.character(id), length(id) == 1L)

    filename <- system.file(
        package = "GenomeInfoDb", "extdata", "dataFiles",
        "genomeMappingTbl.csv")
    builds <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
    idx <- rowSums(builds[,1:2] == tolower(organism)) > 0
    if (!any(idx))
        stop("could not match organism ", organism,
             "; see GenomeInfoDb::listSpecies()")
    builds <- builds[idx,]

    species <- tail(builds$organism, 1L)
    twoletter <- sub("([A-z]).* ([a-z]).*", "\\U\\1\\L\\2", species, perl=TRUE)
    binomial <- sub("([A-z]).* ([[:alpha:]]+)", "\\U\\1\\L\\2", species,
                    perl=TRUE)
    org <- sprintf("org.%s.eg.db", twoletter)

    if (missing(genome))
        genome <- tail(builds$ucscID, 1L)

    if (missing(id)) {
        ids <- c("knownGene", "ensGene", "refGene")
        pkgs <- grep("TxDb", row.names(installed.packages()), value=TRUE)
        found <- FALSE
        for (id in ids) {
            txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, genome, id)
            if (txdb %in% pkgs) {
                found <- TRUE
                break
            }
        }
        if (!found)
            stop("could not guess TxDb package for:",
                 "\n  organism = ", organism,
                 "\n  genome = ", genome,
                 "\n  id = ", paste(ids, collapse=", "))
    }
    txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, genome, id)

    if (verbose)
        message("using ", org, ", ", txdb)
    src_organism(org, txdb)
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
    sql <- "SELECT name FROM sqlite_temp_master WHERE type IN ('view', 'table')"
    dbGetQuery(x$con, sql)$name
}

#' @importFrom dplyr tbl
#' @export
tbl.src_organism <- function(src, ...) {
    tbl <- NextMethod(src, ...)
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}
