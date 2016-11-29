#' Create a sqlite database from an org package and a txdb package
#'
#' The database provides a convenient way to map between gene, transcript,
#' and protein identifiers.
#'
#' This function is meant to be a building block for
#' \code{\link{src_organism}}, which provides an integrated
#' presentation of identifiers and genomic coordinates.
#'
#' Create a dplyr database integrating org.* and TxDb.* information
#'     
#' @param txdb character(1) naming a \code{TxDb.*} package (e.g.,
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#'     
#' @param dbpath path and file name where SQLite file will be accessed
#'      or created if not already exists.
#' 
#' @return A dplyr \code{src_sqlite} instance representing the data
#'     tables.
#' 
#' @rdname src_organism
#' 
#' @importFrom RSQLite dbGetQuery dbConnect dbDisconnect SQLite dbWriteTable
#' @importFrom DBI dbListTables
#' @importFrom AnnotationDbi dbconn dbfile 
#' @importFrom S4Vectors metadata
#' @importFrom methods is as
#' @importFrom tools file_ext 
#' @importFrom dplyr tbl 
#' @importFrom data.table setDT
#' @importFrom GenomeInfoDb as.data.frame
#' 
#' @examples
#' # create human sqlite database with TxDb.Hsapiens.UCSC.hg38.knownGene and 
#' # corresponding org.Hs.eg.db
#' organism <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' src_tbls(organism)
#' tbl(organism, "id")
#' 
#' @export
src_organism <- function(txdb=NULL, dbpath=NULL) {
    if (is.null(txdb)) {
        if (!(!is.null(dbpath) && file.exists(dbpath) && 
            file_ext(dbpath) == "sqlite"))
            stop("input valid sqlite file path or specify 'txdb'")
    } else {
        ## check org, txdb
        if (is.character(txdb))
            txdb <- loadNamespace(txdb)[[txdb]]
        stopifnot(is(txdb, "TxDb"))
        
        txdb_name <- basename(dbfile(txdb))
        org <- strsplit(txdb_name, "[.]")[[1]][2]
        org <-
            if (org != "Mmulatta") {
                sprintf("org.%s.eg.db", substr(org, 1, 2))
            } else "org.Mmu.eg.db"
        # org <- sprintf("org.%s.eg.db", substr(org, 1, 2))
        org <- loadNamespace(org)[[org]]
        
        org_meta <- metadata(org)
        org_schema <- org_meta$value[org_meta$name == "DBSCHEMA"]
        txdb_meta <- metadata(txdb)
        txdb_type <- txdb_meta$value[txdb_meta$name == "Type of Gene ID"]
        
        ## check metadata of org, txdb match tables in dbpath
        if ((!is.null(dbpath)) && file.exists(dbpath)) {
            con <- dbConnect(SQLite(), dbpath)
            src <- src_sql("sqlite", con, path=dbfile(con))
            txdb_md <- tbl(src, "metadata_txdb") %>% collect()
            if (!identical(txdb_meta$value, txdb_md$value))
                stop("'txdb', 'dbpath' sqlite schemas differ; use new dbpath?")
        }
        
        if (startsWith(tolower(txdb_type), "entrez")) {
            txdb_schema <- "txdb_entrez"
        } else if (startsWith(tolower(txdb_type), "ensembl")) {
            txdb_schema <- "txdb_ensembl"
        } else
            stop("unknown TxDb type: ", sQuote(txdb_type))
    } 
    
    ## check dbpath
    if (is.null(dbpath)) {
        dbpath <- file.path(tempdir(), txdb_name)
    }

    ## create db connection
    found <- file.exists(dbpath)
    con <- dbConnect(SQLite(), dbpath)
    if (!found) {
        message("creating 'src_organism' database...")
        withCallingHandlers({
            .add_view(con, org, org_schema)
            .add_view(con, txdb, txdb_schema)
            
            ## create seqinfo table
            seqinfo <- GenomicFeatures:::get_TxDb_seqinfo0(txdb)
            seqinfo <- as.data.frame(seqinfo)
            setDT(seqinfo, keep.rownames = TRUE)[]
            names(seqinfo) <- 
                c("seqnames", "seqlengths", "isCircular", "genome")
            dbWriteTable(con, "seqinfo", seqinfo)
            
            ## check ensembl length
            if (txdb_schema == "txdb_ensembl")
                .alter_ensembl_ids(con)
        }, error=function(e) {
            dbDisconnect(con)         # clean-up
            file.remove(dbpath)
            stop(e)                   # signal condition to user
        })
    }
    
    src <- src_sql("sqlite", con, path=dbfile(con))
    schema <- tail(colnames(tbl(src, "ranges_gene")), 1L)
    src$schema <-
        if (startsWith(schema, "entrez")) {
            "entrez"
        } else "ensembl"
    class(src) <- c("src_organism", class(src))
    src
}

.get_ensembl_id_len <- function(con, tablename) {
    dbGetQuery(con, 
               paste0("SELECT LENGTH(ensembl) FROM ", tablename, 
                      " WHERE ensembl is NOT NULL LIMIT 1")) [[1]]
}

.alter_table <- function(con, tablename, length) {
    dbGetQuery(con,
               paste0("ALTER TABLE ", tablename, 
                      " ADD COLUMN ensembl_original CHARACTER"))
    dbGetQuery(con,
               paste0("UPDATE ", tablename, 
                      " SET ensembl_original = ensembl"))
    dbGetQuery(con,
               paste0("UPDATE ", tablename, 
                      " SET ensembl = SUBSTR(ensembl, 1, ", length, ")"))
}

.alter_ensembl_ids <- function(con) {
    ## e.g., FBgn from org.Dm* has FBgn000xx, whereas
    ## TxDb.Dmelanogaster has FBgn000xx.1. Create a new column for
    ## original, upadate ensembl column to be like org.*
    orgens <- .get_ensembl_id_len(con, "id")
    txdbens <- .get_ensembl_id_len(con, "ranges_gene")
    if (orgens == txdbens)
        return()
    
    tbls <- dbListTables(con)
    if (orgens < txdbens) {
        tbls <- tbls[startsWith(tbls, "ranges")]
        for (tablename in tbls)
            .alter_table(con, tablename, orgens)
    } else {
        tbls <- tbls[startsWith(tbls, "id")]
        for (tablename in tbls)
            .alter_table(con, tablename, txdbens)
    }
}

.add_view <- function(con, db, tblname) {
    tbls <- dbGetQuery(con, "pragma database_list;")$name
    if (tblname %in% tbls)
        return()
    else {
        dbGetQuery(con, paste0("ATTACH '", dbfile(db), "' AS ", tblname))
        
        fname <- system.file(
            package="Organism.dplyr", "schema", paste0(tblname, ".sql"))
        if (!file.exists(fname))
            stop(sQuote(tblname), " schema not unsupported")
        
        schemas <- readLines(fname)
        grps <- cumsum(!nzchar(schemas)) + 1
        for (schema in split(schemas, grps)) {
            sql <- paste(schema, collapse="\n")
            dbGetQuery(con, sql)
        }
        dbGetQuery(con, paste0("DETACH '", tblname, "'"))
    }
}


#' @param organism organism name
#' 
#' @param genome genome name
#' 
#' @param id choose from "knownGene", "ensGene" and "refGene"
#' 
#' @param verbose logical. Should R report extra information on progress? 
#' Default is TRUE.
#' 
#' @examples 
#' # create human sqlite database using hg38 genome
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
    src_organism(txdb, dbpath)
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


#' @param .data A tbl.
#' 
#' @param ... Comma separated list of unquoted expressions. You can treat 
#' variable names like they are positions. Use positive values to select 
#' variables; use negative values to drop variables.
#' 
#' @importFrom dplyr src_sql "%>%" tbl select_
#' 
#' @rdname src_organism
#' @export
select_.tbl_organism <- function(.data, ...) { 
    .data = NextMethod(.data, ...) 
    dplyr::distinct_(.data, ...) 
}

#' @param src A src_organism object
#' 
#' @examples 
#' # Look at all available tables
#' src_tbls(organism)
#' 
#' @importFrom dplyr src_tbls
#' @importFrom RSQLite dbGetQuery
#' @rdname src_organism
#' @export
src_tbls.src_organism <- function(x) {
    sql <- "SELECT name FROM sqlite_master WHERE type IN ('view', 'table') 
            AND (SUBSTR(name,1,2) = 'id' OR SUBSTR(name,1,6) = 'ranges')"
    dbGetQuery(x$con, sql)$name
}

#' @examples 
#' # Look at data in table "id"
#' tbl(organism, "id")
#' 
#' @importFrom dplyr tbl
#' @rdname src_organism
#' @export
tbl.src_organism <- function(src, ...) {
    tbl <- NextMethod(src, ...)
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}

setOldClass("src_organism")
