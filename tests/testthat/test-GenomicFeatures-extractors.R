context("GenomicFeatures-extractors")

suppressPackageStartupMessages({
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

hg38light <- hg38light()
src <- src_organism(dbpath=hg38light)

.test_extractor <- function(src, txdb, fun, subset) {
    suppressWarnings({
        src <- fun(src)
        txdb <- fun(txdb)
    })
    expect_true(all(mcols(src)[[subset]] %in% mcols(txdb)[[subset]]))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))

    txdb <- txdb[match(mcols(src)[[subset]], mcols(txdb)[[subset]], 0)]
    o_src <- GenomicRanges::order(granges(src))
    o_txdb <- GenomicRanges::order(granges(txdb))

    ## FIXME: metadata differs, so simpler expect_identical(tx_src,
    ## tx_txdb) fails
    expect_identical(granges(src)[o_src], granges(txdb)[o_txdb])
    expect_identical(mcols(src)[o_src,], mcols(txdb)[o_txdb,])
}

.test_extractor_egfilter <- function(src, txdb, fun, subset) {
    egid <- c("10", "100")
    src0 <- src
    txdb <- fun(txdb, filter=list(gene_id=egid))

    ## AnnotationFilterList(*Filter)
    src <- fun(src0, filter=AnnotationFilterList(EntrezFilter(egid)))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(mcols(src)[[subset]], mcols(txdb)[[subset]]))

    ## list(*Filter)
    src <- fun(src0, filter=list(EntrezFilter(egid)))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(mcols(src)[[subset]], mcols(txdb)[[subset]]))

    ## #Filter
    src <- fun(src0, filter=EntrezFilter(egid))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(mcols(src)[[subset]], mcols(txdb)[[subset]]))
}

.test_extractorBy <- function(src, txdb, funBy) {
    src <- funBy(src)
    expect_true(all(names(src) %in% names(txdb)))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))

    txdb <- txdb[names(src)]
    o_src <- GenomicRanges::order(granges(unlist(src)))
    o_txdb <- GenomicRanges::order(granges(unlist(txdb)))

    expect_identical(granges(unlist(src))[o_src], granges(unlist(txdb)[o_txdb]))
    expect_identical(mcols(unlist(src))[o_src, ], mcols(unlist(txdb))[o_txdb, ])
}

.test_extractorBy_txfilter <- function(src, txdb, funBy) {
    txid <- c(15880L, 15881L)
    src0 <- src
    txdb <- txdb[txid]

    ## AnnotationFilterList(*Filter)
    src <- funBy(src0, filter=AnnotationFilterList(TxIdFilter(txid)))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## list(*Filter)
    src <- funBy(src0, filter=list(TxIdFilter(txid)))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## *Filter
    src <- funBy(src0, filter=TxIdFilter(txid))
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))
}

test_that("validate-filter", {
    an1 <- AnnotationFilterList(SymbolFilter("ADA"), SeqNameFilter('NFkB'))
    an2 <- AnnotationFilterList(SymbolFilter("ADA"), TxEndFilter(1000000, '<'))
    expect_false(.check_filters(src, an1))
    expect_true(.check_filters(src, an2))
})

test_that("transcripts-extractor", {
    .test_extractor(src, txdb, transcripts, "tx_name")
    .test_extractor_egfilter(src, txdb, transcripts, "tx_id")
})

test_that("exons-extractor", {
    .test_extractor(src, txdb, exons, "exon_id")
    .test_extractor_egfilter(src, txdb, exons, "exon_id")
})

test_that("cds-extractor", {
    .test_extractor(src, txdb, cds, "cds_id")
    .test_extractor_egfilter(src, txdb, cds, "cds_id")
})

test_that("promoters-extractor", {
    .test_extractor(src, txdb, promoters, "tx_id")
    .test_extractor_egfilter(src, txdb, promoters, "tx_id")
})

test_that("transcriptsBy-extractor", {
    txdb <- suppressWarnings(transcriptsBy(txdb))
    .test_extractorBy(src, txdb, transcriptsBy)
    ## .test_extractorBy_txfilter(src, txdb, transcriptsBy)
    
    ## filters
    ## FIXME TxIdFilter does not work correctly here
    egid <- c("10", "100")
    tx_src <- unlist(transcriptsBy(src, filter=AnnotationFilterList(EntrezFilter(egid))))
    tx_txdb <- unlist(txdb[egid])
    expect_equal(length(tx_src), length(tx_txdb))
    expect_true(all.equal(seqinfo(tx_src), seqinfo(tx_txdb)))
    expect_equal(tx_src$tx_id, tx_txdb$tx_id)
})

test_that("exonsBy-extractor", {
    txdb <- suppressWarnings(exonsBy(txdb))
    .test_extractorBy(src, txdb, exonsBy)
    .test_extractorBy_txfilter(src, txdb, exonsBy)
})

test_that("cdsBy-extractor", {
    txdb <- suppressWarnings(cdsBy(txdb))
    .test_extractorBy(src, txdb, cdsBy)
    .test_extractorBy_txfilter(src, txdb, cdsBy)
})

test_that("intronsByTranscript-extractor", {
    txdb <- suppressWarnings(intronsByTranscript(txdb))
    .test_extractorBy(src, txdb, intronsByTranscript)
    .test_extractorBy_txfilter(src, txdb, intronsByTranscript)
})    

test_that("fiveUTRsByTranscript-extractor", {
    txdb <- suppressWarnings(fiveUTRsByTranscript(txdb))
    .test_extractorBy(src, txdb, fiveUTRsByTranscript)
    .test_extractorBy_txfilter(src, txdb, fiveUTRsByTranscript)
})

test_that("threeUTRsByTranscript-extractor", {
    txdb <- suppressWarnings(threeUTRsByTranscript(txdb))
    .test_extractorBy(src, txdb, threeUTRsByTranscript)
    .test_extractorBy_txfilter(src, txdb, threeUTRsByTranscript)
})    

#test_that(".logicOp_subset works", {
#    logicalOp_subset <- Organism.dplyr:::.logicOp_subset
#
#    expect_equal(logicalOp_subset(character(), FALSE), character())
#    expect_equal(logicalOp_subset(character(), TRUE), character())
#
#    expect_equal(logicalOp_subset("&", c(TRUE, TRUE)), "&")
#    expect_equal(logicalOp_subset("&", c(TRUE, FALSE)), character())
#    expect_equal(logicalOp_subset("&", c(FALSE, TRUE)), character())
#    expect_equal(logicalOp_subset("&", c(FALSE, FALSE)), character())
#
#    expect_equal(logicalOp_subset("|", c(TRUE, TRUE)), "|")
#    expect_equal(logicalOp_subset("|", c(TRUE, FALSE)), character())
#    expect_equal(logicalOp_subset("|", c(FALSE, TRUE)), character())
#    expect_equal(logicalOp_subset("|", c(FALSE, FALSE)), character())
#	
#    expect_equal(logicalOp_subset(c("&", "&"), c(TRUE, TRUE, TRUE)), c("&", "&"))
#    expect_equal(logicalOp_subset(c("&", "&"), c(FALSE, TRUE, TRUE)), "&")
#    expect_equal(logicalOp_subset(c("&", "&"), c(TRUE, TRUE, FALSE)), "&")
#    expect_equal(logicalOp_subset(c("&", "&"), c(TRUE, FALSE, TRUE)), "&")
#    expect_equal(logicalOp_subset(c("&", "&"), c(FALSE, TRUE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("&", "&"), c(FALSE, FALSE, TRUE)), character())
#    expect_equal(logicalOp_subset(c("&", "&"), c(TRUE, FALSE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("&", "&"), c(FALSE, FALSE, FALSE)), character())
#
#    expect_equal(logicalOp_subset(c("|", "|"), c(TRUE, TRUE, TRUE)), c("|", "|"))
#    expect_equal(logicalOp_subset(c("|", "|"), c(FALSE, TRUE, TRUE)), "|")
#    expect_equal(logicalOp_subset(c("|", "|"), c(TRUE, TRUE, FALSE)), "|")
#    expect_equal(logicalOp_subset(c("|", "|"), c(TRUE, FALSE, TRUE)), "|")
#    expect_equal(logicalOp_subset(c("|", "|"), c(FALSE, TRUE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("|", "|"), c(FALSE, FALSE, TRUE)), character())
#    expect_equal(logicalOp_subset(c("|", "|"), c(TRUE, FALSE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("|", "|"), c(FALSE, FALSE, FALSE)), character())
#
#    expect_equal(logicalOp_subset(c("&", "|"), c(TRUE, TRUE, TRUE)), c("&", "|"))
#    expect_equal(logicalOp_subset(c("&", "|"), c(FALSE, TRUE, TRUE)), "|")
#    expect_equal(logicalOp_subset(c("|", "&"), c(TRUE, TRUE, FALSE)), "|")
#    expect_equal(logicalOp_subset(c("|", "&"), c(TRUE, FALSE, TRUE)), "&")
#    expect_equal(logicalOp_subset(c("&", "|"), c(FALSE, TRUE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("&", "|"), c(FALSE, FALSE, TRUE)), character())
#    expect_equal(logicalOp_subset(c("|", "&"), c(TRUE, FALSE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("|", "&"), c(FALSE, FALSE, FALSE)), character())
#
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(TRUE, TRUE, TRUE, TRUE)), c("&", "|", "&"))
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(FALSE, TRUE, TRUE, TRUE)), c("|", "&"))
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(TRUE, TRUE, FALSE, TRUE)), c("&", "&"))
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(TRUE, FALSE, TRUE, FALSE)), "&")
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(FALSE, TRUE, FALSE, FALSE)), character())
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(FALSE, FALSE, TRUE, TRUE)), "&")
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(TRUE, FALSE, FALSE, TRUE)), "&")
#    expect_equal(logicalOp_subset(c("&", "|", "&"), c(FALSE, FALSE, FALSE, FALSE)), character())
#})
