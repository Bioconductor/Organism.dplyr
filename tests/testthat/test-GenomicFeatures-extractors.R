context("GenomicFeatures-extractors")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

test_that("transcripts-extractor", {
    tx_src <- transcripts(src)
    tx_txdb <- transcripts(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))

    ## filters
    tx_src <- transcripts(src, filter=list(entrez=c("5728", "672")))
    tx_txdb <- transcripts(txdb, filter=list(gene_id=c("5728", "672")))
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$tx_id, tx_txdb$tx_id)
})

test_that("exons-extractor", {
    tx_src <- exons(src)
    tx_txdb <- exons(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- exons(src, filter=list(entrez=c("5728", "672")))
    tx_txdb <- exons(txdb, filter=list(gene_id=c("5728", "672")))
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})

test_that("cds-extractor", {
    tx_src <- cds(src)
    tx_txdb <- cds(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- cds(src, filter=list(entrez=c("5728", "672")))
    tx_txdb <- cds(txdb, filter=list(gene_id=c("5728", "672")))
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$cds_id, tx_txdb$cds_id)
})

test_that("promoters-extractor", {
    tx_src <- promoters(src)
    tx_txdb <- promoters(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- promoters(src, filter=list(entrez="5728"))
    tx_txdb <- promoters(txdb, filter=list(gene_id="5728"))
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$tx_id, tx_txdb$tx_id)
})

test_that("transcriptsBy-extractor", {
    tx_src <- transcriptsBy(src)
    tx_txdb <- transcriptsBy(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(transcriptsBy(src, filter=list(entrez="5728")))
    tx_txdb <- unlist(transcriptsBy(txdb)["5728"])
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$tx_id, tx_txdb$tx_id)
})

test_that("exonsBy-extractor", {
    tx_src <- exonsBy(src)
    tx_txdb <- exonsBy(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(exonsBy(src, filter=list(tx_id="87017")))
    tx_txdb <- unlist(exonsBy(txdb)["87017"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})

test_that("cdsBy-extractor", {
    tx_src <- cdsBy(src)
    tx_txdb <- cdsBy(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(cdsBy(src, filter=list(tx_id="87017")))
    tx_txdb <- unlist(cdsBy(txdb)["87017"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$cds_id, tx_txdb$cds_id)
})


test_that("intronsByTranscript-extractor", {
    tx_src <- intronsByTranscript(src)
    tx_txdb <- intronsByTranscript(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(intronsByTranscript(src, filter=list(tx_id="87017")))
    tx_txdb <- unlist(intronsByTranscript(txdb)["87017"])
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
})

test_that("fiveUTRsByTranscript-extractor", {
    tx_src <- fiveUTRsByTranscript(src)
    tx_txdb <- fiveUTRsByTranscript(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(fiveUTRsByTranscript(src, filter=list(tx_id="87011")))
    tx_txdb <- unlist(fiveUTRsByTranscript(txdb)["87011"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})

test_that("threeUTRsByTranscript-extractor", {
    tx_src <- threeUTRsByTranscript(src)
    tx_txdb <- threeUTRsByTranscript(txdb)
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_true(identical(mcols(tx_src), mcols(tx_txdb)))
    
    ## filters
    tx_src <- unlist(threeUTRsByTranscript(src, filter=list(tx_id="87017")))
    tx_txdb <- unlist(threeUTRsByTranscript(txdb)["87017"])
    expect_equal(length(unlist(tx_src)), length(unlist(tx_txdb)))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$exon_id, tx_txdb$exon_id)
})
