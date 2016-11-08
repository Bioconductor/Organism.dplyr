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
#' @param dbpath path and file name where SQLite file will be accessed
#'      or created if not already exists.
#' 
#' @return A dplyr \code{tbl_sql} instance representing the data
#'     table.
#'     
#' @importFrom RSQLite dbSendQuery dbGetQuery dbConnect dbDisconnect SQLite
#' @importFrom AnnotationDbi dbconn dbfile taxonomyId
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom dplyr tbl
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
#' inner_join(id, tbl(organism, "ranges_tx")) %>% 
#'      filter(symbol == "PTEN") %>%
#'      dplyr::select(entrez, symbol, tx_id, tx_start, tx_end)
#' 
#' 
#' @export
src_organism <- function(org=NULL, txdb=NULL, dbpath=NULL) {
    if (is.null(org) && is.null(txdb)) {
        if (!(!is.null(dbpath) && file.exists(dbpath) && 
            file_ext(dbpath) == "sqlite"))
            stop("input valid sqlite file path")
    } else if (!is.null(org) && !is.null(txdb)) {
        ## check org, txdb
        if (is.character(org))
            org <- loadNamespace(org)[[org]]
        stopifnot(is(org, "OrgDb"))
        
        if (is.character(txdb))
            txdb <- loadNamespace(txdb)[[txdb]]
        stopifnot(is(txdb, "TxDb"))
        
        stopifnot(identical(taxonomyId(org), taxonomyId(txdb)))
        if ((!is.null(dbpath)) && file.exists(dbpath)) {
            ## FIXME: check metadata of org, txdb match tables in dbpath
        }
        
        org_meta <- metadata(org)
        org_schema <- org_meta$value[org_meta$name == "DBSCHEMA"]
        txdb_meta <- metadata(txdb)
        txdb_type <- txdb_meta$value[txdb_meta$name == "Type of Gene ID"]
        
        if (startsWith(tolower(txdb_type), "entrez")) {
            txdb_schema <- "txdb_entrez"
        } else if (startsWith(tolower(txdb_type), "ensembl")) {
            txdb_schema <- "txdb_ensembl"
        } else
            stop("unknown TxDb type: ", sQuote(txdb_type))
    } else
        stop("specify 'org' and 'txdb' or 'dbpath'")
    
    ## check dbpath
    if (is.null(dbpath)) {
        org_name <- file_path_sans_ext(basename(dbfile(org)))
        txdb_name <- file_path_sans_ext(basename(dbfile(txdb)))
        filename <- sprintf("%s_%s.sqlite", org_name, txdb_name)
        dbpath <- file.path(tempdir(), filename)
    }

    ## create db connection
    found <- file.exists(dbpath)
    con <- dbConnect(SQLite(), dbpath)
    if (!found) {
        message("creating 'src_organism' database...")
        withCallingHandlers({
            .add_view(con, org, org_schema)
            .add_view(con, txdb, txdb_schema)
        }, error=function(e) {
            dbDisconnect(con)         # clean-up
            file.remove(dbpath)
            stop(e)                   # signal condition to user
        })
    }
    
    src <- src_sql("sqlite", con, path=dbfile(con))
    src$schema <- tail(colnames(tbl(src, "ranges_gene")), 1)
    class(src) <- c("src_organism", class(src))
    src
}

.add_view <- function(con, db, tblname) {
    tbls <- dbGetQuery(con, "pragma database_list;")$name
    if (tblname %in% tbls)
        return()
    else {
        # ifelse(is(db, "OrgDb"), dbname <- "org", dbname <- tblname)
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
        dbGetQuery(con, paste0("DETACH '", tblname, "'"))
    }
}

# .add_view <- function(con, db, tblname) {
#     tbls <- dbGetQuery(con, "pragma database_list;")$name
#     if (tblname %in% tbls)
#         return()
#     if (is(db, "TxDb"))
#         dbGetQuery(con, paste0("ATTACH '", dbfile(db), "' AS ", tblname))
#     
#     fname <- system.file(
#         package="Organism.dplyr", "schema", paste0(tblname, ".sql"))
#     if (!file.exists(fname))
#         stop(sQuote(tblname), " schema not unsupported")
# 
#     schemas <- readLines(fname)
#     grps <- cumsum(!nzchar(schemas)) + 1
#     for (schema in split(schemas, grps)) {
#         sql <- paste(schema, collapse="\n")
#         dbSendQuery(con, sql)
#     }
#     if (is(db, "TxDb"))
#         dbGetQuery(con, paste0("DETACH '", tblname, "'"))
# }

#' @param organism
#' 
#' @param genome
#' 
#' @param id
#' 
#' @examples 
#' human <- src_ucsc("human")
#' 
#' @rdname src_organism
#' @importFrom GenomeInfoDb genomeBuilds
#' @importFrom utils read.csv tail installed.packages
#' @export
src_ucsc <- function(organism, genome = NULL, id = NULL, 
                     dbpath=NULL, verbose=TRUE) {
    stopifnot(is.character(organism), length(organism) == 1L)
    if (!missing(genome))
        stopifnot(is.character(genome), length(genome) == 1L)
    if (!missing(id))
        stopifnot(is.character(id), length(id) == 1L)
    
    # OrgDb
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
    twoletter <- 
        sub("([A-z]).* ([a-z]).*", "\\U\\1\\L\\2", species, perl=TRUE)
    binomial <- sub("([A-z]).* ([[:alpha:]]+)", "\\U\\1\\L\\2", species,
                    perl=TRUE)
    org <- sprintf("org.%s.eg.db", twoletter)
    
    # TxDb
    pkgs <- grep("TxDb", row.names(installed.packages()), value=TRUE)
    pkgs <- grep(binomial, pkgs, value=TRUE)
    if (length(pkgs) == 0)
        stop("No TxDb package available for organism: ", organism)
    
    found <- FALSE
    if (missing(id) & missing(genome)) {
        ids <- c("knownGene", "ensGene", "refGene")
        genome <- tail(builds$ucscID, 1L)
        for (id in ids) {
            txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, genome, id)
            if (txdb %in% pkgs) {
                found <- TRUE
                break
            }
        }
    } else if (missing(id)) {
        txdb <- .findPkg(pkgs, genome)
    } else if (missing(genome)) {
        txdb <- .findPkg(pkgs, id)
    } else {
        txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, genome, id)
        if (txdb %in% pkgs) 
            found <- TRUE
    }
    
    if (!is.null(txdb)) 
        found <- TRUE
    
    if (!found)
        stop("could not guess TxDb package for:",
             "\n  organism = ", organism,
             if (!missing(genome))
                 "\n  genome = ", genome,
             if (!missing(id))
                 "\n  genome = ", id)
    
    if (verbose)
        message("using ", org, ", ", txdb)
    src_organism(org, txdb, dbpath)
}

.findPkg <- function(pkgs, var) {
    txdb <- NULL
    pkgs <- grep(var, pkgs, value=TRUE)
    if (length(pkgs) != 0) {
        txdb = ifelse(length(pkgs) == 1, 
                      pkgs, 
                      names(tail(sapply(pkgs, function(pkg) 
                          unlist(strsplit(pkg, "[.]"))[4]), 1L)))
        }
    txdb
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
    sql <- "SELECT name FROM sqlite_master WHERE type IN ('view', 'table')"
    dbGetQuery(x$con, sql)$name
}

#' @importFrom dplyr tbl
#' @export
tbl.src_organism <- function(src, ...) {
    tbl <- NextMethod(src, ...)
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}

setOldClass("src_organism")
