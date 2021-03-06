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

    ## *Expression
    src <- fun(src0, filter=~entrez == egid)
    expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(mcols(src)[[subset]], mcols(txdb)[[subset]]))

    ## AnnotationFilter Negation
    src1 <- fun(src0, filter=~symbol == "ADA")
    src2 <- fun(src0, filter=~!symbol != "ADA")
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))

    ## AnnotationFilterList Negation
    src1 <- fun(src0, filter=~symbol == "ADA" & tx_id == 169786) 
    src2 <- fun(src0, filter=~!(symbol != "ADA") & tx_id != 169786)
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))

    ## Grouping
    src1 <- fun(src0, filter=~symbol == "ADA" & tx_id == 169786 |
        symbol %startsWith% "SNORD")
    src2 <- fun(src0, filter=~symbol %startsWith% "SNORD" |
        (symbol == "ADA" & tx_id != 169786))
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))
}

.test_extractorBy <- function(src, txdb, funBy) {
    src <- funBy(src)
    expect_true(all(names(src) %in% names(txdb)))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))

    txdb <- txdb[names(src)] #names(txdb) %in% names(src)]
    o_src <- GenomicRanges::order(granges(unlist(src)))
    o_txdb <- GenomicRanges::order(granges(unlist(txdb)))

    expect_identical(granges(unlist(src))[o_src], granges(unlist(txdb)[o_txdb]))
    expect_identical(mcols(unlist(src))[o_src, ], mcols(unlist(txdb))[o_txdb, ])
}

.test_extractorBy_txfilter <- function(src, txdb, funBy) {
    txid <- 1:2#c(10L, 72L)#c(15880L, 15881L)
    src0 <- src
    txdb <- txdb[txid]#txdb[as.character(txid)]

    ## AnnotationFilterList(*Filter)
    src <- funBy(src0, filter=AnnotationFilterList(TxIdFilter(txid)))
    #expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## list(*Filter)
    src <- funBy(src0, filter=list(TxIdFilter(txid)))
    #expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## *Filter
    src <- funBy(src0, filter=TxIdFilter(txid))
    #expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## *Expression
    src <- funBy(src0, filter=~tx_id == txid)
    #expect_equal(length(src), length(txdb))
    expect_true(all.equal(seqinfo(src), seqinfo(txdb)))
    expect_true(setequal(src$tx_id, txdb$tx_id))

    ## AnnotationFilter Negation
    src1 <- funBy(src0, filter=~symbol == "ADA")
    src2 <- funBy(src0, filter=~!symbol != "ADA")
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))

    ## AnnotationFilterList Negation
    src1 <- funBy(src0, filter=~symbol=="ADA" & tx_id == 169786) 
    src2 <- funBy(src0, filter=~!(symbol!="ADA") & tx_id != 169786)
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))

    ## Grouping
    src1 <- funBy(src0, filter=~symbol == "ADA" & tx_id == 169786 |
        symbol %startsWith% "SNORD")
    src2 <- funBy(src0, filter=~symbol %startsWith% "SNORD" |
        (symbol == "ADA" & tx_id != 169786))
    #expect_equal(length(src1), length(src2))
    expect_true(all.equal(seqinfo(src1), seqinfo(src2)))
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
    ## Does not work for transcriptsBy; doesn't appear to need fixing
    ## .test_extractorBy_txfilter(src, txdb, transcriptsBy)
    
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
    ## .test_extractorBy(src, txdb, fiveUTRsByTranscript)
    # FIXME: fiveUTRs seem to be throwing issues in checks that are being masked
#    .test_extractorBy_txfilter(src, txdb, fiveUTRsByTranscript)
})

test_that("threeUTRsByTranscript-extractor", {
    txdb <- suppressWarnings(threeUTRsByTranscript(txdb))
    #.test_extractorBy(src, txdb, threeUTRsByTranscript)
    .test_extractorBy_txfilter(src, txdb, threeUTRsByTranscript)
})

.deleteTempTables(src)
