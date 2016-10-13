.tbl_filter <- function(keep1, table, filter) {
    values <- paste0("'", filter[[keep1]], "'", collapse=", ")
    op <- if (length(filter[[keep1]]) == 1) "==" else "%in%"
    sprintf("%s %s c(%s)", keep1, op, values)
}

.tbl_join <- function(x, table, filter) {
    if (is.null(filter))
        return(table)

    if ("entrez" %in% names(filter)) {
        filters <- .tbl_filter("entrez", table, filter)
        table <- table %>% filter_(filters)
        filter <- filter[names(filter) != "entrez"]
    }

    fields <- names(filter)
    tbls <- src_tbls(x)

    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || length(keep) == 0)
            next
        filters <- sapply(keep, .tbl_filter, table, filter)
        filters <- paste0(filters, collapse=" & ")
        table <- inner_join(table, tbl(x, i)) %>% filter_(filters)
    }

    table
}

#' Generic functions to extract genomic features from an object.
#'
#' @param x A src_organism object
#'
#' @param filter Either NULL or a named list of vectors to be used to
#'     restrict the output.
#'
#' @examples
#' organism <- src_organism("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
#' filters <- list(symbol=c("PTEN", "BRCA1"),
#'                entrez="5728",
#'                go=c("GO:0000079", "GO:0001933"))
#' transcripts(organism, filters)
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures transcripts
#' @export

setMethod("transcripts", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_tx")
    table <- .tbl_join(x, table, filter)
    fields <- unique(
        c("tx_id", "tx_chrom", "tx_start", "tx_end", "tx_name",
          names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})

#' @rdname src_organism
#' @importFrom GenomicFeatures exons
#' @export

setMethod("exons", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_exon")
    .tbl_join(x, table, filter) %>%
        dplyr::select(entrez, exon_chrom, exon_start, exon_end,
                      exon_strand, tx_id, exon_id, exon_name, exon_rank)
})


#' @rdname src_organism
#' @importFrom GenomicFeatures cds
#' @export

setMethod("cds", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_cds")
    .tbl_join(x, table, filter) %>%
        dplyr::select(entrez, cds_chrom, cds_start, cds_end,
                      cds_strand, tx_id, cds_id, cds_name, exon_rank)
})


#' @rdname src_organism
#' @importFrom GenomicFeatures genes
#' @export

setMethod("genes", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_gene")
    .tbl_join(x, table, filter) %>%
        dplyr::select(entrez, tx_chrom, gene_start, gene_end,
                      tx_strand)
})


#' @param upstream For \code{promoters()}: An integer(1) value indicating the
#' number of bases upstream from the transcription start site.
#'
#' @param downstream For \code{promoters()}: An integer(1) value indicating
#' the number of bases downstream from the transcription start site.
#'
#' @examples
#' promoters(organism, upstream=100, downstream=50,
#'           filter = list(symbol="BRCA1"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures promoters
#' @importFrom dplyr mutate
#' @export

setMethod("promoters", "src_organism",
function(x, upstream, downstream, filter = NULL) {
    if(missing(upstream))
        upstream = 0
    if(missing(downstream))
        downstream = 2200

    tmp <- transcripts(organism, filter = filter) %>%
           mutate(start = if(tx_strand == "+" | tx_strand == "*")
                              tx_start - upstream
                          else if (tx_strand == "-")
                              tx_end - downstream + 1,
                  end = if(tx_strand == "+" | tx_strand == "*")
                            tx_start + downstream - 1
                        else if (tx_strand == "-")
                            tx_end + upstream)

    tmp %>%
    dplyr::select(tx_id, entrez, tx_chrom, start, end, tx_strand, tx_name)
})


#' @param by One of "gene", "exon", "cds" or "tx". Determines the grouping.
#'
#' @examples
#' transcriptsBy(organism, by = "gene", filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures transcriptsBy
#' @importFrom dplyr inner_join
#' @export

setMethod("transcriptsBy", "src_organism",
function(x, by = c("gene", "exon", "cds"), filter = NULL) {
    by <- match.arg(by)
    tx <- transcripts(x, filter = filter)

    if (by == "gene") {
        inner_join(tx, tbl(x, "ranges_gene")) %>%
            dplyr::select(entrez, tx_chrom, tx_start, tx_end,
                          tx_strand, tx_id, tx_name)
    }
    else if (by == "exon") {
        inner_join(tx, tbl(x, "ranges_exon"), by = "tx_id") %>%
            dplyr::select(exon_id, tx_chrom, tx_start, tx_end,
                          tx_strand, tx_id, tx_name, exon_rank)
    }
    else if (by == "cds") {
        inner_join(tx, tbl(x, "ranges_cds"), by = "tx_id") %>%
            dplyr::select(cds_id, tx_chrom, tx_start, tx_end,
                          tx_strand, tx_id, tx_name, exon_rank)
    }
})


#' @examples
#' exonsBy(organism, by = "gene", filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures exonsBy
#' @export

setMethod("exonsBy", "src_organism",
function(x, by = c("gene", "tx"), filter = NULL) {
    by <- match.arg(by)
    exons <- exons(x, filter = filter)

    if (by == "gene") {
        inner_join(exons, tbl(x, "ranges_gene")) %>%
            dplyr::select(entrez, exon_chrom, exon_start, exon_end,
                          exon_strand, exon_id, exon_name)
    }
    else if (by == "tx") {
        inner_join(exons, tbl(x, "ranges_tx"), by = "tx_id") %>%
            dplyr::select(tx_id, exon_chrom, exon_start, exon_end,
                          exon_strand, exon_id, exon_name, exon_rank)
    }
})


#' @examples
#' cdsBy(organism, by = "gene", filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures cdsBy
#' @export

setMethod("cdsBy", "src_organism",
function(x, by = c("gene", "tx"), filter = NULL) {
    by <- match.arg(by)
    cds <- cds(x, filter = filter)

    if (by == "gene") {
        inner_join(cds, tbl(x, "ranges_gene")) %>%
            dplyr::select(entrez, cds_chrom, cds_start, cds_end,
                          cds_strand, cds_id, cds_name)
    }
    else if (by == "tx") {
        inner_join(cds, tbl(x, "ranges_tx"), by = "tx_id") %>%
            dplyr::select(tx_id, cds_chrom, cds_start, cds_end,
                          cds_strand, cds_id, cds_name, exon_rank)
    }
})
