#' @importFrom BiocFileCache BiocFileCache bfcquery bfcrpath bfcnew
#'     bfccache bfcremove

.src_organism_dbpath <-
    function(dbpath, txdb_name, txdb_meta, overwrite)
{
    bfc <- NULL
    if (is.null(dbpath))
        bfc <- BiocFileCache()

    if (is(bfc, "BiocFileCache")) {
        ## get dbpath from BiocFileCache
        nrec <- NROW(bfcquery(bfc, txdb_name, "rname", exact = TRUE))
        if (nrec == 0L) {
            dbpath <- bfcnew(bfc, txdb_name)
        } else if (nrec == 1L) {
            dbpath <- bfcrpath(bfc, txdb_name)
        } else {
            stop(
                "\n  'bfc' contains duplicate Organism.dplyr record names",
                "\n      bfccache(): '", bfccache(bfc), "'",
                "\n      rname: '", txdb_name, "'"
            )
        }
    }

    if (!file.exists(dbpath))
        return(dbpath)

    ## check metadata of org, txdb match tables in dbpath
    con <- dbConnect(SQLite(), dbpath)
    on.exit(dbDisconnect(con))
    src <- src_dbi(con)
    txdb_md <- tbl(src, "metadata_txdb") %>% collect()

    if (identical(txdb_meta$value, txdb_md$value)) {
        return(dbpath)
    } else if (!overwrite) {
        stop(
            "\n  'txdb', 'dbpath' sqlite schemas differ",
            "\n  use `overwrite = TRUE` or a different `dbpath = `"
        )
    }

    message("overwriting existing dbpath")
    if (is(bfc, "BiocFileCache")) {
        bfcremove(bfc, names(dbpath))
        dbpath <- bfcnew(bfc, txdb_name)
    } else {
        unlink(dbpath)
    }

    dbpath
}

#' Create a sqlite database from TxDb and corresponding Org packages
#'
#' The database provides a convenient way to map between gene, transcript,
#' and protein identifiers.
#'
#' \code{src_organism()} and \code{src_ucsc()} are meant to be a building block
#' for \code{\link{src_organism}}, which provides an integrated
#' presentation of identifiers and genomic coordinates.
#'
#' \code{src_organism()} creates a dplyr database integrating org.* and TxDb.*
#' information by given TxDb. And \code{src_ucsc()} creates the database by
#' given organism name, genome and/or id.
#'
#' supportedOrganisms() provides all supported organisms in this package with
#' corresponding OrgDb and TxDb.
#'
#' @param txdb character(1) naming a \code{TxDb.*} package (e.g.,
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or a \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#'
#' @param dbpath character(1) path or BiocFileCache instance
#'     representing the location where an Organism.dplyr SQLite
#'     database will be accessed or created. If no path is specified,
#'     the SQLite file is created in the default BiocFileCache()
#'     location.
#'
#' @param overwrite logical(1) overwrite an exisging `dbpath` contains
#'     an Organism.dplyr SQLite databse different from the version
#'     implied by `txdb`?
#'
#' @return \code{src_organism()} and \code{src_ucsc()} returns a dplyr
#'     \code{src_dbi} instance representing the data tables.
#'
#' @seealso \code{\link{dplyr}} for details about using \code{dplyr} to
#'     manipulate data.
#'
#'     \code{\link{transcripts_tbl}} for generic functions to extract
#'     genomic features from a \code{src_organism} object.
#'
#'     \code{\link{select,src_organism-method}} for "select" interface
#'     on \code{src_organism} objects.
#'
#' @author Yubo Cheng.
#'
#' @rdname src
#'
#' @importFrom  dplyr %>% arrange as.tbl collect compute desc distinct filter
#'     full_join inner_join is.tbl left_join mutate order_by rename right_join
#'     select_ src src_tbls summarise summarize tbl union union_all
#' @importFrom dbplyr build_sql src_sql src_dbi tbl_sql
#' @importFrom RSQLite dbGetQuery dbConnect dbDisconnect SQLite
#'     dbWriteTable dbListTables
#' @importFrom S4Vectors metadata
#' @importFrom methods is as
#' @importFrom tibble as_tibble
#' @importFrom tools file_ext
#' @importFrom AnnotationDbi dbfile
#' @importFrom GenomeInfoDb as.data.frame
#'
#' @examples
#' ## create human sqlite database with TxDb.Hsapiens.UCSC.hg38.knownGene and
#' ## corresponding org.Hs.eg.db
#' \dontrun{src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")}
#' src <- src_organism(dbpath=hg38light())
#'
#' ## query using dplyr
#' inner_join(tbl(src, "id"), tbl(src, "id_go")) %>%
#'      filter(symbol == "ADA") %>%
#'      dplyr::select(entrez, ensembl, symbol, go, evidence, ontology)
#'
#' @export
src_organism <- function(txdb=NULL, dbpath=NULL, overwrite=FALSE) {
    if (is.null(txdb)) {
        if (!(!is.null(dbpath) && file.exists(dbpath) &&
            file_ext(dbpath) == "sqlite"))
            stop("input valid sqlite file path or specify 'txdb'")
    } else {
        ## check org, txdb
        if (is.character(txdb) && length(txdb) == 1L)
            txdb <- loadNamespace(txdb)[[txdb]]
        stopifnot(is(txdb, "TxDb"))

        txdb_name <- paste0("dplyr.", basename(dbfile(txdb)))
        org <- strsplit(txdb_name, ".", fixed=TRUE)[[1]][3]
        org <-
            sprintf("org.%s.eg.db", substr(org, 1, 2 + (org == "Mmulatta")))
        org <- loadNamespace(org)[[org]]

        org_meta <- metadata(org)
        org_schema <- org_meta$value[org_meta$name == "DBSCHEMA"]
        txdb_meta <- metadata(txdb)
        txdb_type <- txdb_meta$value[txdb_meta$name == "Type of Gene ID"]

        dbpath <- .src_organism_dbpath(dbpath, txdb_name, txdb_meta, overwrite)

        if (startsWith(tolower(txdb_type), "entrez")) {
            txdb_schema <- "txdb_entrez"
        } else if (startsWith(tolower(txdb_type), "ensembl")) {
            txdb_schema <- "txdb_ensembl"
        } else
            stop("unknown TxDb type: ", sQuote(txdb_type))
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
            seqinfo$seqnames <- rownames(seqinfo)
            cols <- c("seqnames", "seqlengths", "isCircular", "genome")
            dbWriteTable(con, "seqinfo", seqinfo[, cols])

            ## check ensembl length
            if (txdb_schema == "txdb_ensembl")
                .alter_ensembl_ids(con)
        }, error=function(e) {
            dbDisconnect(con)         # clean-up
            file.remove(dbpath)
            stop(e)                   # signal condition to user
        })
    }

    src <- src_dbi(con)
    schema <- tail(colnames(tbl(src, "ranges_gene")), 1L)
    src$schema <-
        if (startsWith(schema, "entrez")) {
            "entrez"
        } else "ensembl"
    class(src) <- c("src_organism", class(src))
    ## TODO: Add support for using permanent database for table storage.
    ## Need to ensure no tables generated are kept inside permanent DB.
    if (FALSE) #file.access(src$dbpath, mode=2) == 0)
        src$db <- src$con
    else
        src$db <- dbConnect(RSQLite::SQLite(), "")
    src$dbpath <- src$db@dbname
    src
}

.get_ensembl_id_len <- function(con, tablename) {
    dbGetQuery(con,
               paste0("SELECT LENGTH(ensembl) FROM ", tablename,
                      " WHERE ensembl is NOT NULL LIMIT 1")) [[1]]
}

.alter_table <- function(con, tablename, length) {
    dbExecute(con,
              paste0("ALTER TABLE ", tablename,
                     " ADD COLUMN ensembl_original CHARACTER"))
    dbExecute(con,
              paste0("UPDATE ", tablename,
                     " SET ensembl_original = ensembl"))
    dbExecute(con,
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
        dbExecute(con, paste0("ATTACH '", dbfile(db), "' AS ", tblname))

        fname <- system.file(
            package="Organism.dplyr", "schema", paste0(tblname, ".sql"))
        if (!file.exists(fname))
            stop(sQuote(tblname), " schema not unsupported")

        schemas <- trimws(readLines(fname))
        grps <- cumsum(!nzchar(schemas)) + 1
        for (schema in split(schemas, grps)) {
            sql <- paste(schema, collapse="\n")
            dbExecute(con, sql)
        }
        dbExecute(con, paste0("DETACH '", tblname, "'"))
    }
}

.src_ucsc_builds <- function(organism) {
    filename <- system.file(
        package = "GenomeInfoDb", "extdata", "dataFiles",
        "genomeMappingTbl.csv"
    )
    builds <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
    idx <- rowSums(builds[,1:2] == tolower(organism)) > 0
    if (!any(idx))
        stop(
            "\n",
            "  could not match organism '", organism, "';\n",
            "  see 'commonName' field of 'GenomeInfoDb::listOrganisms()'"
        )
    builds[idx,]
}

.src_ucsc_organism <- function(builds) {
    unique(builds$organism)
}

.src_ucsc_binomial <- function(builds) {
    organism <- .src_ucsc_organism(builds)
     sub("([A-z]).* ([[:alpha:]]+)", "\\U\\1\\L\\2", organism, perl=TRUE)
}

.src_ucsc_org <- function(builds) {
    organism <- .src_ucsc_organism(builds)
    twoletter <- sub("([A-z]).* ([a-z]).*", "\\U\\1\\L\\2", organism, perl=TRUE)
    sprintf("org.%s.eg.db", twoletter)
}

.src_ucsc_ids <- function(builds) {
    ids <- unique(builds$ucscID)
    id_n <- as.integer(sub("^[^[:digit:]]*", "", ids))
    ids[order(id_n, decreasing = TRUE)]
}

.src_ucsc_txdb_packages <- function(organism, binomial) {
    pkgs <- grep("TxDb", row.names(installed.packages()), value=TRUE)
    pkgs <- grep(binomial, pkgs, value=TRUE)
    if (length(pkgs) == 0)
        stop(
            "\n",
            "  could not find installed TxDb package for '", organism, "'\n",
            "    binomial = ", binomial, "\n",
            "see 'BiocManager::available(\"TxDb\")'"
        )
    pkgs
}

.src_ucsc_missing_id_and_genome <- function(organism, builds, pkgs) {
    binomial <- .src_ucsc_binomial(builds)
    ids <- .src_ucsc_ids(builds) # all possible genomes
    resources <- c("knownGene", "ensGene", "refGene") # specified in docs

    found <- FALSE
    for (id in ids) {
        for (resource in resources) {
            txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, id, resource)
            if (txdb %in% pkgs) {
                found <- TRUE
                break
            }
        }
        if (found)
            break
    }

    if (!found)
        stop(
            "\n",
            "  could not find installed TxDb package for '", organism, "'\n",
            "    binomial = ", binomial, "\n",
            "    genomes = " , paste(ids, collapse = ", "), "\n",
            "    resources = ", paste(resources, collapse = ", "), "\n",
            "  see 'BiocManager::available(\"TxDb\")'"
        )

    txdb
}

#' @param organism organism or common name
#'
#' @param genome genome name
#'
#' @param id choose from "knownGene", "ensGene" and "refGene"
#'
#' @param verbose logical. Should R report extra information on progress?
#' Default is TRUE.
#'
#' @examples
#' ## create human sqlite database using hg38 genome
#' \dontrun{human <- src_ucsc("human")}
#'
#' @rdname src
#' @importFrom GenomeInfoDb genomeBuilds
#' @importFrom utils read.csv tail installed.packages
#' @export
src_ucsc <- function(organism, genome = NULL, id = NULL,
                     dbpath=NULL, verbose=TRUE) {
    stopifnot(is.character(organism), length(organism) == 1L, !is.na(organism))
    if (!missing(genome))
        stopifnot(is.character(genome), length(genome) == 1L)
    if (!missing(id))
        stopifnot(is.character(id), length(id) == 1L)

    ## OrgDb
    builds <- .src_ucsc_builds(organism)
    binomial <- .src_ucsc_binomial(builds)
    org <- .src_ucsc_org(builds)

    ## TxDb
    pkgs <- .src_ucsc_txdb_packages(organism, binomial)

    if (missing(id) && missing(genome)) {
        txdb <- .src_ucsc_missing_id_and_genome(organism, builds, pkgs)
    } else if (missing(id)) {
        txdb <- .findPkg(pkgs, genome)
    } else if (missing(genome)) {
        txdb <- .findPkg(pkgs, id)
    } else {
        txdb <- sprintf("TxDb.%s.UCSC.%s.%s", binomial, genome, id)
        if (!txdb %in% pkgs)
            txdb <- NULL
    }


    if (is.null(txdb))
        stop("could not guess TxDb package for:",
             "\n  organism = ", organism,
             if (!missing(genome))
                 "\n  genome = ", genome,
             if (!missing(id))
                 "\n  id = ", id)

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


#' @examples
#' ## all supported organisms with corresponding OrgDb and TxDb
#' supportedOrganisms()
#'
#' @rdname src
#' @export
supportedOrganisms <- function() {
    filename <- system.file(
        package = "Organism.dplyr", "extdata",
        "supportedOrganisms.csv")
    csv <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
    tibble::as_tibble(csv)
}

#' @param .data A tbl.
#'
#' @description DEPRECATED: Please use equivalent select() method.
#'
#' @param ... Comma separated list of unquoted expressions. You can treat
#' variable names like they are positions. Use positive values to select
#' variables; use negative values to drop variables.
#'
#' @rdname src
#' @export
select_.tbl_organism <- function(.data, ...) {
    .Deprecated("select")
    .data = NextMethod(.data, ...)
    dplyr::distinct(.data, ...)
}

#' @param x A src_organism object
#'
#' @examples
#' ## Look at all available tables
#' src_tbls(src)
#'
#' @importFrom RSQLite dbGetQuery dbExecute
#' @rdname src
#' @export
src_tbls.src_organism <- function(x, ...) {
    sql <- "SELECT name FROM sqlite_master WHERE type IN ('view', 'table')
            AND (SUBSTR(name,1,2) = 'id' OR SUBSTR(name,1,6) = 'ranges')"
    dbGetQuery(x$con, sql)$name
}

#' @param src An src_organism object
#'
#' @param .load_tbl_only a logic(1) that indicates whether only to load
#'  the table instead of also loading the pacakge in the temporary database.
#'  Default value is FALSE.
#'
#' @return A tibble of the requested table coming from the temporary database
#'  of the src_organism object.
#'
#' @examples
#' ## Look at data in table "id"
#' tbl(src, "id")
#'
#' ## Look at fields of one table
#' colnames(tbl(src, "id"))
#'
#' @rdname src
#' @importFrom DBI dbListTables dbWriteTable
#' @export
tbl.src_organism <- function(src, ..., .load_tbl_only = FALSE) {
    args <- list(...)
    if (length(args) == 0)
        stop("tbl name required.")
    table <- args[[1]]
    tbl <- NextMethod(src, ...)
    if (!.load_tbl_only && src$dbpath == "") {
        if (!(table %in% dbListTables(src$db))) {
            tbl <- tbl %>% collect(n=Inf)
            dbWriteTable(src$db, table, tbl)
        }
        tbl <- tbl(src$db, table)
    }
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}


setOldClass("src_organism")

.getSeqinfo <- function(x) {
    seqinfo <- mutate(tbl(x, "seqinfo") %>% collect(n=Inf),
                      seqnames = as.character(.data$seqnames),
                      seqlengths = as.integer(.data$seqlengths),
                      isCircular = as.logical(.data$isCircular),
                      genome = as.character(.data$genome))
    Seqinfo(seqnames=seqinfo[['seqnames']], seqlengths=seqinfo[['seqlengths']],
            isCircular=seqinfo[['isCircular']], genome=seqinfo[['genome']])
}

#' @examples
#' ## name of org package of src_organism object
#' orgPackageName(src)
#'
#' @importFrom AnnotationDbi orgPackageName
#' @importFrom dplyr pull
#' @rdname src
#' @exportMethod orgPackageName
setMethod("orgPackageName", "src_organism",
    function(x)
{
    org <-
        tbl(x, "metadata_txdb") %>%
        filter(.data$name == "Organism") %>%
        pull("value")
    supportedOrganisms() %>%
        filter(.data$organism == org) %>%
        pull("OrgDb") %>%
        unique()
})

#' @examples
#' ## seqinfo of src_organism object
#' seqinfo(src)
#'
#' @importFrom GenomeInfoDb Seqinfo seqinfo
#' @rdname src
#' @exportMethod seqinfo
setMethod("seqinfo", "src_organism",
    function(x)
{
    .getSeqinfo(x)
})
