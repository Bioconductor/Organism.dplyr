context("src_organism-class")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

test_that("src_organism_constructor", {
    expect_error(src_organism())
    expect_error(src_organism("FOO"))

    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
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

    src <- src_ucsc("human")
    expect_equal(is(src, "src_organism"), TRUE)

    ## prevent overwrite
    expect_error(src_ucsc("human", genome = "hg37"))
})

test_that("mouse", {
    ## test knownGene
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    src <- src_ucsc("mouse", id = "knownGene")
    
    tx_src <- transcriptsBy(src)
    tx_txdb <- transcriptsBy(txdb)
    
    expect_equal(is(src, "src_organism"), TRUE)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_equal(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)), TRUE)
    
    ## test ensGene no filter
    library(TxDb.Mmusculus.UCSC.mm10.ensGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
    src <- src_ucsc("mouse", id = "ensGene")
    
    tx_src <- transcriptsBy(src)
    tx_txdb <- transcriptsBy(txdb)
    expect_equal(is(src, "src_organism"), TRUE)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_equal(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)), TRUE)
    
    ## test ensGene with filter
    tx_src <- unlist(exonsBy(src, filter=list(Tx_idFilter("3"))))
    tx_txdb <- unlist(exonsBy(txdb)["3"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_equal(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)), TRUE)
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})
