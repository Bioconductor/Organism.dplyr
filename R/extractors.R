# transcripts(), exons(), cds(), genes(), promoters(), transcriptsBy(), exonsBy(), 
# cdsBy(), intronsByTranscript(), fiveUTRsByTranscript() and threeUTRsByTranscript()

# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")

# transcripts()
#' @param src object created by src_organism
#' 
#' @examples
#' organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
#' inner_join(transcripts1(organism), tbl(organism, "id"), copy = TRUE) %>% 
#'      filter(symbol == "PTEN") %>%
#'      dplyr::select(entrez, chrom, txstart, txend, txid, txname)
#' 
#' @rdname src_organism
#' @export
transcripts1 <- function(src) {
    tbl(src, "ranges_tx")  %>%
        dplyr::select(entrez, chrom, txstart, txend, txid, txname)
}


#' @examples
#' transcripts2(organism) %>%
#'      filter(symbol == "PTEN") %>%
#'      dplyr::select(entrez, chrom, txstart, txend, txid, txname)
#' 
#' @rdname src_organism
#' @export
transcripts2 <- function(src) {
    inner_join(tbl(src, "ranges_tx"), tbl(src, "id"))
}

 

   
    


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
