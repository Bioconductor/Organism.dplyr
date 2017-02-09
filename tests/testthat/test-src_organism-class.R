context("src_organism-class")

hg38light <- system.file(
    package="Organism.dplyr", "extdata", "light.hg38.knownGene.sqlite"
)

mm10light <- system.file(
    package="Organism.dplyr", "extdata", "light.mm10.ensGene.sqlite"
)

test_that("src_organism_constructor", {
    expect_error(src_organism())
    expect_error(src_organism("FOO"))

    if (interactive()) {
        message("creating expensive 'src_organism'")
        src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
    } else {
        src <- src_organism(dbpath=hg38light)
    }
    expect_equal(is(src, "src_organism"), TRUE)
    expect_equal(length(src), 3)
    expect_equal(is(src$con, "SQLiteConnection"), TRUE)
    expect_equal(class(src$path), "character")
    expect_equal(class(src$schema), "character")
    
    ## prevent overwrite
    expect_error(src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene", src$path))
    
    ## test for expected tables
    tbls <- src_tbls(src)
    expect_equal(length(tbls[startsWith(tbls, "id")]), 7)
    expect_equal(length(tbls[startsWith(tbls, "ranges")]), 4)
    
    ## test for tables not empty
    for (table in tbls) 
        expect_true(dim(head(tbl(src, table)) %>% collect())[1] > 0)
})

test_that("src_ucsc_constructor", {
    expect_error(src_ucsc())
    expect_error(src_ucsc("FOO"))

    if (interactive()) {
        message("creating expensive 'src_organism'")
        src <- src_ucsc("human")
    } else {
        src <- src_organism(dbpath=hg38light)
    }
    expect_equal(is(src, "src_organism"), TRUE)

    ## prevent overwrite
    expect_error(src_ucsc("human", genome = "hg37"))
})

test_that("mouse", {
    ## test ensGene no filter
    suppressPackageStartupMessages({
        library(TxDb.Mmusculus.UCSC.mm10.ensGene)
    })
    txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
    src <- src_organism(dbpath=mm10light)
    expect_true(is(src, "src_organism"))
    
    tx_src <- transcriptsBy(src)
    tx_txdb <- transcriptsBy(txdb)[names(tx_src)]
    expect_true(all(names(tx_src) %in% names(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    
    ## test ensGene with filter
    tx_src <- unlist(exonsBy(src, filter=list(TxIdFilter("2237"))))
    tx_txdb <- unlist(exonsBy(txdb)["2237"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})
