library(TxDb.Hsapiens.UCSC.hg38.knownGene)

test_that("extractors return same result as txdb", {
    src <- src_organism(dbpath =  "C:/Users/YU19864/Documents/test.sqlite")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    expect_equal(length(transcripts(src)), length(transcripts(txdb)))
    expect_equal(length(exons(src)), length(exons(txdb)))
    expect_equal(length(cds(src)), length(cds(txdb)))
    expect_equal(length(promoters(src)), length(promoters(txdb)))
    expect_equal(length(transcriptsBy(src)), length(transcriptsBy(txdb)))
    expect_equal(length(exonsBy(src)), length(exonsBy(txdb)))
    expect_equal(length(cdsBy(src)), length(cdsBy(txdb)))
    expect_equal(length(intronsByTranscript(src)), 
                 length(intronsByTranscript(txdb)))
    expect_equal(length(fiveUTRsByTranscript(src)), 
                 length(fiveUTRsByTranscript(txdb)))
    expect_equal(length(threeUTRsByTranscript(src)), 
                 length(threeUTRsByTranscript(txdb)))
})