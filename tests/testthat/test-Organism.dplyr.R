library(Organism.dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

test_that("input txdb and dbpath consistent", {
    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
    expect_error(src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene", src$path))
    expect_error(src_ucsc("human", genome = "hg37"))
})

test_that("extractors return same result as txdb", {
    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
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

test_that("extractors return same result as txdb with filter", {
    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    expect_equal(length(transcripts(src, filter=list(entrez = "5728"))), 
                 length(transcripts(txdb, filter=list(gene_id = "5728"))))
    expect_equal(length(exons(src, filter=list(entrez = "5728"))), 
                 length(exons(txdb, filter=list(gene_id = "5728"))))
    expect_equal(length(cds(src, filter=list(entrez = "5728"))), 
                 length(cds(txdb, filter=list(gene_id = "5728"))))
    expect_equal(length(promoters(src, filter=list(entrez = "5728"))), 
                 length(promoters(txdb, filter=list(gene_id = "5728"))))
    expect_equal(length(unlist(transcriptsBy(src, filter=list(entrez = "5728")))), 
                 length(unlist(transcriptsBy(txdb)["5728"])))
})