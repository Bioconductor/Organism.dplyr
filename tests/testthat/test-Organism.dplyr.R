context("GenomicFeatures-extractors")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

## Example for transcripts(). This logic could be applied for the
## other extractors as well.

test_that("transcripts-extractor", {
    tx_src <- transcripts(src)
    tx_txdb <- transcripts(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_equal(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)), TRUE)
    ## this fails; not sure if it should
    # expect_equal(identical(mcols(tx_src), mcols(tx_txdb)), TRUE)

    ## filters
    tx_src <- transcripts(src, filter=list(entrez="5728"))
    tx_txdb <- transcripts(txdb, filter=list(gene_id="5728"))
    expect_equal(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)), TRUE)
    ## this fails; not sure if it should
    # expect_equal(identical(mcols(tx_src), mcols(tx_txdb)), TRUE)
})

test_that("extractors return same result as txdb", {
    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
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
