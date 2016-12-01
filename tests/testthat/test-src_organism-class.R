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

    ## test for expected tables?
    ## test for expected variable in 'id' table?
})

test_that("src_ucsc_constructor", {
    expect_error(src_ucsc())
    expect_error(src_ucsc("FOO"))

    src <- src_ucsc("human")
    expect_equal(is(src, "src_organism"), TRUE)

    ## prevent overwrite
    expect_error(src_ucsc("human", genome = "hg37"))

    ## Should the tables in the object created with src_ucsc("human") be
    ## the same as those in an object created with
    ## src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")?
    ## If yes, you could add a few tests to confirm that.
})
