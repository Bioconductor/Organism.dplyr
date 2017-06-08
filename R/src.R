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
#'     \code{TxDb.Hsapiens.UCSC.hg38.knownGene}) or \code{TxDb}
#'     object instantiating the content of a \code{TxDb.*} pacakge.
#'
#' @param dbpath path and file name where SQLite file will be accessed
#'      or created if not already exists.
#'
#' @return \code{src_organism()} and \code{src_ucsc()} returns a dplyr
#'     \code{src_sqlite} instance representing the data tables.
#'
#' @seealso \code{\link{dplyr}} for details about using \code{dplyr} to
#'     manipulate data.
#'
#'     \code{\link[Organism.dplyr]{transcripts_tbl}} for generic functions to
#'     extract genomic features from a \code{src_organism} object.
#'
#'     \code{\link[Organism.dplyr]{select,src_organism-method}} for "select"
#'     interface on \code{src_organism} objects.
#'
#' @author Yubo Cheng.
#'
#' @rdname src
#'
#' @importFrom  dplyr %>% arrange arrange_ as.tbl build_sql collect compute
#'     desc distinct distinct_ filter filter_ full_join group_by group_by_
#'     inner_join is.tbl left_join mutate mutate_ order_by rename rename_
#'     right_join select_ src src_sql src_sqlite src_tbls summarise summarise_
#'     summarize summarize_ tbl tbl_df tbl_sql union union_all
#' @importFrom RSQLite dbGetQuery dbConnect dbDisconnect SQLite
#'     dbWriteTable dbListTables
#' @importFrom S4Vectors metadata
#' @importFrom methods is as
#' @importFrom tools file_ext
#' @importFrom AnnotationDbi dbfile
#' @importFrom GenomeInfoDb as.data.frame
#' @importFrom BiocFileCache BiocFileCache bfccache
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
src_organism <- function(txdb=NULL, dbpath=NULL) {
    if (is.null(txdb)) {
        if (!(!is.null(dbpath) && file.exists(dbpath) &&
            file_ext(dbpath) == "sqlite"))
            stop("input valid sqlite file path or specify 'txdb'")
    } else {
        ## check org, txdb
        if (is.character(txdb) && length(txdb) == 1L)
            txdb <- loadNamespace(txdb)[[txdb]]
        else
            stop("input one valid 'txdb'")
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

        ## check dbpath
        if (is.null(dbpath)) {
            dbpath <- file.path(bfccache(BiocFileCache()), txdb_name)
        }

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
    tbl_df(csv)
}

#' @param .data A tbl.
#'
#' @param ... Comma separated list of unquoted expressions. You can treat
#' variable names like they are positions. Use positive values to select
#' variables; use negative values to drop variables.
#'
#' @rdname src
#' @export
select_.tbl_organism <- function(.data, ...) {
    .data = NextMethod(.data, ...)
    dplyr::distinct_(.data, ...)
}

#' @param x A src_organism object
#'
#' @examples
#' ## Look at all available tables
#' src_tbls(src)
#'
#' @importFrom RSQLite dbGetQuery
#' @rdname src
#' @export
src_tbls.src_organism <- function(x) {
    sql <- "SELECT name FROM sqlite_master WHERE type IN ('view', 'table')
            AND (SUBSTR(name,1,2) = 'id' OR SUBSTR(name,1,6) = 'ranges')"
    dbGetQuery(x$con, sql)$name
}

#' @param src A src_organism object
#'
#' @examples
#' ## Look at data in table "id"
#' tbl(src, "id")
#'
#' ## Look at fields of one table
#' colnames(tbl(src, "id"))
#'
#' @rdname src
#' @export
tbl.src_organism <- function(src, ...) {
    tbl <- NextMethod(src, ...)
    class(tbl) <- c("tbl_organism", class(tbl))
    tbl
}

setOldClass("src_organism")

.getSeqinfo <- function(x) {
    seqinfo <- mutate_(tbl(x, "seqinfo") %>% collect(n=Inf),
                       seqnames = ~ as.character(seqnames),
                       seqlengths = ~ as.integer(seqlengths),
                       isCircular = ~ as.logical(isCircular),
                       genome = ~ as.character(genome))
    Seqinfo(seqnames=seqinfo[[1]], seqlengths=seqinfo[[2]],
            isCircular=seqinfo[[3]], genome=seqinfo[[4]])
}

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
