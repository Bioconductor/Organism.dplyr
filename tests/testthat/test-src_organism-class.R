context("src_organism-class")

test_that("src_organism_constructor", {
    expect_error(src_organism())
    expect_error(src_organism("FOO"))

    if (interactive()) {
        message("creating expensive 'src_organism'")
        src <- src_organism(
            "TxDb.Hsapiens.UCSC.hg38.knownGene", overwrite = TRUE
        )
    } else {
        src <- src_organism(dbpath=hg38light())
    }
    expect_equal(is(src, "src_organism"), TRUE)
    expect_equal(length(src), 5)
    expect_equal(is(src$con, "SQLiteConnection"), TRUE)
    expect_equal(class(src$schema), "character")
    
    ## prevent overwrite
    expect_error(src_organism(
        "TxDb.Hsapiens.UCSC.hg19.knownGene",
        dbfile(src$con)
    ))
    
    ## test for expected tables
    tbls <- src_tbls(src)
    expect_equal(length(tbls[startsWith(tbls, "id")]), 7)
    expect_equal(length(tbls[startsWith(tbls, "ranges")]), 4)
    
    ## test for tables not empty
    for (table in tbls) 
        expect_true(dim(head(tbl(src, table)) %>% collect())[1] > 0)
    .deleteTempTables(src)
})

test_that(".src_ucsc_builds() and friends work", {
    expect_error(.src_ucsc_builds(""))
    expect_error(.src_ucsc_builds("humanz"))

    human <- .src_ucsc_builds("human")
    expect_true(NROW(human) > 1L)
    expect_identical(.src_ucsc_organism(human), "homo sapiens")
    expect_identical(.src_ucsc_binomial(human), "Hsapiens")
    expect_identical(.src_ucsc_org(human), "org.Hs.eg.db")
    expect_identical(.src_ucsc_ids(human), c("hg38", "hg19", "hs1"))

    mouse <- .src_ucsc_builds("mouse")
    expect_true(NROW(mouse) > 1L)
    expect_identical(.src_ucsc_organism(mouse), "mus musculus")
    expect_identical(.src_ucsc_binomial(mouse), "Mmusculus")
    expect_identical(.src_ucsc_org(mouse), "org.Mm.eg.db")
    expect_identical(.src_ucsc_ids(mouse), c("mm39", "mm10", "mm9"))
})

test_that(".src_ucsc_txdb_packages() works", {
    ## assumes that Suggests: packages are installed
    organism <- "human"
    builds <- .src_ucsc_builds(organism)
    binomial <- .src_ucsc_binomial(builds)
    pkgs <- .src_ucsc_txdb_packages(organism, binomial)
    expect_true("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% pkgs)

    organism <- "mouse"
    builds <- .src_ucsc_builds(organism)
    binomial <- .src_ucsc_binomial(builds)
    pkgs <- .src_ucsc_txdb_packages(organism, binomial)
    expect_true("TxDb.Mmusculus.UCSC.mm10.ensGene" %in% pkgs)
})

test_that(".src_ucsc_missing_id_and_genome() works", {
    organism <- "human"
    builds <- .src_ucsc_builds(organism)
    binomial <- .src_ucsc_binomial(builds)
    pkgs <- .src_ucsc_txdb_packages(organism, binomial)
    txdb <- .src_ucsc_missing_id_and_genome(organism, builds, pkgs)
    expect_true(startsWith(txdb, "TxDb.Hsapiens.UCSC"))

    organism <- "mouse"
    builds <- .src_ucsc_builds(organism)
    binomial <- .src_ucsc_binomial(builds)
    pkgs <- .src_ucsc_txdb_packages(organism, binomial)
    txdb <- .src_ucsc_missing_id_and_genome(organism, builds, pkgs)
    expect_true(startsWith(txdb, "TxDb.Mmusculus.UCSC"))
})

test_that("src_ucsc_constructor", {
    expect_error(src_ucsc())
    expect_error(src_ucsc("FOO"))

    if (interactive()) {
        message("creating expensive 'src_organism'")
        src <- src_ucsc("human")
    } else {
        src <- src_organism(dbpath=hg38light())
    }
    expect_equal(is(src, "src_organism"), TRUE)

    ## prevent overwrite
    expect_error(src_ucsc("human", genome = "hg37"))
    .deleteTempTables(src)
})

test_that("mouse", {
    ## test ensGene no filter
    suppressPackageStartupMessages({
        library(TxDb.Mmusculus.UCSC.mm10.ensGene)
    })
    txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
    src <- src_organism(dbpath=mm10light())
    expect_true(is(src, "src_organism"))
    .deleteTempTables(src)
    
    tx_src <- transcriptsBy(src)
    tx_txdb <- transcriptsBy(txdb)[names(tx_src)]
    expect_true(all(names(tx_src) %in% names(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    
    ## test ensGene with filter
    tx_src <- unlist(exonsBy(src, filter=AnnotationFilterList(TxIdFilter(2237))))
    tx_txdb <- unlist(exonsBy(txdb)["2237"])

    ## Removed `unlist()` from test as of 1.7.3 due to it not modifying the
    ## GRanges output.
    expect_equal(length(tx_src), length(tx_txdb))

    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
    .deleteTempTables(src)
})

