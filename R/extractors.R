# transcripts(), exons(), cds(), genes(), promoters(),
# transcriptsBy(), exonsBy(), cdsBy(), intronsByTranscript(),
# fiveUTRsByTranscript() and threeUTRsByTranscript()

# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")

# transcripts()

#' Generic functions to extract genomic features from an object.
#' 
#' @param x object created by src_organism
#' 
#' @param filter Either NULL or a named list of vectors to be used to 
#'     restrict the output.
#'     
#' @examples
#' organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
#' filter <- list(symbol=c("PTEN", "BRCA1"), 
#'                entrez="5728", 
#'                go=c("GO:0000079", "GO:0001933"))
#' transcripts(organism, filter)
#' 
#' @rdname src_organism
#' @importFrom GenomicFeatures transcripts
#' @export

setMethod("transcripts", "src_organism", function(x, filter) {
    fields <- names(filter)
    tbls <- src_tbls(x)
    table <- tbl(x, "ranges_tx")

    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || (length(keep) == 1L && keep == "entrez"))
            next
        filters <- sapply(keep, function(keep1) {
            values <- paste0(sQuote(filter[[keep1]]), collapse=", ")
            op <- if (length(filter[[keep1]]) == 1) "==" else "%in%"
            sprintf("%s %s c(%s)", keep1, op, values)
        })
        filters <- paste0(filters, collapse=" & ") 
        table <- inner_join(table, tbl(x, i)) %>% filter_(filters)
    }

    fields <- 
        c("entrez", "tx_chrom", "tx_start", "tx_end", "tx_id", "tx_name",
          fields)
    fields <- paste(sQuote(fields), collapse=", ")
    table %>% dplyr::select_(fields)
})


.tbl_join <- function(src, x, y)
    inner_join(tbl(src, x), tbl(src, y))

transcripts4 <- function(src, x)
    .tbl_join(src, x, "ranges_tx")

# # exons()
# exons <- function(src) {
#     tbl(src, "ranges_exon") %>%
#         dplyr::select(entrez, chrom, exonstart, exonend, strand, exonid)
# }
# 
# inner_join(exons(organism), tbl(organism, "id"), copy = TRUE) %>% 
#     filter(symbol == "PTEN") %>%
#     dplyr::select(chrom, exonstart, exonend, strand, exonid)
# 
# 
# # cds()
# cds <- function(src) {
#     tbl(src, "ranges_cds") %>%
#         dplyr::select(entrez, chrom, cdsstart, cdsend, strand, cdsid)
# }
# 
# inner_join(cds(organism), tbl(organism, "id"), copy = TRUE) %>% 
#     filter(symbol == "PTEN") %>%
#     dplyr::select(chrom, cdsstart, cdsend, strand, cdsid)
# 
# 
# # genes()
# genes <- function(src) {
#     tbl(src, "ranges_gene") %>%
#         dplyr::select(entrez, chrom, genestart, geneend, strand)
# }
# 
# inner_join(cds(organism), tbl(organism, "id"), copy = TRUE) %>% 
#     filter(symbol == "PTEN") %>%
#     dplyr::select(chrom, genestart, geneend, strand, entrez)
# 
# 
# # promoters()
# # (start(x) - upstream) to (start(x) + downstream - 1)
# promoters <- function(src, upstream, downstream) {
#     sql <- paste0("SELECT DISTINCT
#                   entrez, 
#                   chrom, 
#                   txstart - ", upstream, " AS start, 
#                   txstart + ", downstream, " - 1 AS end, 
#                   strand, 
#                   ranges_tx.txid, 
#                   txname
#                   FROM ranges_tx
#                   ")
#     dbGetQuery(src$con, sql)
# }
# 
# inner_join(promoters(organism, 100, 50), tbl(organism, "id") %>% dplyr::select(entrez, symbol), copy = TRUE) %>% 
#     filter(symbol == "PTEN") %>%
#     dplyr::select(chrom, start, end, strand, txid, txname)
# 
# 
# # transcriptsBy()
# transcriptsBy <- function(src, by=c("gene", "exon", "cds")) {
#     if (by == "gene") {
#         inner_join(tbl(src, "ranges_gene"), tbl(src, "ranges_tx")) %>%
#             dplyr::select(entrez, chrom, genestart, geneend, strand, txid, txname)
#     }
#     else if (by == "exon") {
#         inner_join(tbl(src, "ranges_tx"), tbl(src, "ranges_exon")) %>%
#             dplyr::select(exonid, chrom, txstart, txend, strand, txid, txname, rank)
#     }
#     else if (by == "cds") {
#         inner_join(tbl(src, "ranges_tx"), tbl(src, "ranges_cds")) %>%
#             dplyr::select(cdsid, chrom, txstart, txend, strand, txid, txname, rank)
#     }
# }
# 
# # exonsBy()
# 
# 
# # cdsBy()
# 
