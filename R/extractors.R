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
#' organism <- src_organism("org.Hs.eg.db", 
#'                          "TxDb.Hsapiens.UCSC.hg38.knownGene")
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
        c("tx_id", "entrez", "tx_chrom", "tx_start", "tx_end", "tx_strand", 
          "tx_name", names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})

#' @rdname src_organism
#' @importFrom GenomicFeatures exons
#' @export

setMethod("exons", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_exon")
    table <- .tbl_join(x, table, filter)
    fields <- unique(
        c("exon_chrom", "exon_start", "exon_end", "exon_strand",
          "entrez", "tx_id", "exon_id", "exon_name", "exon_rank",
          names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})


#' @rdname src_organism
#' @importFrom GenomicFeatures cds
#' @export

setMethod("cds", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_cds")
    table <- .tbl_join(x, table, filter)
    fields <- unique(
        c("cds_id", "tx_id", "entrez", "cds_chrom", "cds_start", "cds_end", 
          "cds_strand", "cds_name", "exon_rank", names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})


#' @rdname src_organism
#' @importFrom GenomicFeatures genes
#' @export

setMethod("genes", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_gene")
    table <- .tbl_join(x, table, filter)
    fields <- unique(
        c("entrez", "tx_chrom", "gene_start", "gene_end", "tx_strand", 
          names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})


#' @param upstream For \code{promoters()}: An integer(1) value indicating 
#' the number of bases upstream from the transcription start site.
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
    if (!isSingleNumber(upstream)) 
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream)) 
        upstream <- as.numeric(upstream)
    if (!isSingleNumber(downstream)) 
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream)) 
        downstream <- as.numeric(downstream)
    if (upstream < 0 | downstream < 0) 
        stop("'upstream' and 'downstream' must be integers >= 0")
    # if (any(strand(x) == "*")) 
    #     warning("'*' ranges were treated as '+'")
    if(missing(upstream))
        upstream = 0
    if(missing(downstream))
        downstream = 2200
    
    table <- transcripts(organism, filter = filter) %>%
        mutate(start = ifelse(tx_strand == "-", 
                              tx_end - downstream + 1, 
                              tx_start - upstream),
            end = ifelse(tx_strand == "-", 
                         tx_end + upstream, 
                         tx_start + downstream - 1))
    
    fields <- unique(
        c("tx_id", "entrez", "tx_chrom", "start", "end", "tx_strand", 
          "tx_name", names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
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
        table <- inner_join(tx, tbl(x, "ranges_gene")) 
        fields <- unique(
            c("entrez", "tx_id", "tx_chrom", "tx_start", "tx_end",  
              "tx_strand", "tx_name", names(filter)))
    }
    else if (by == "exon") {
        table <- inner_join(tx, tbl(x, "ranges_exon"), by = "tx_id") 
        fields <- unique(
            c("exon_id", "tx_id", "tx_chrom", "tx_start", "tx_end",  
              "tx_strand", "tx_name", "exon_rank", names(filter)))
    }
    else if (by == "cds") {
        table <- inner_join(tx, tbl(x, "ranges_cds"), by = "tx_id") 
        fields <- unique(
            c("cds_id", "tx_id", "tx_chrom", "tx_start", "tx_end",  
              "tx_strand", "tx_name", "exon_rank", names(filter)))
    }
    do.call(select_, c(list(table), as.list(fields)))
})


#' @examples
#' exonsBy(organism, by = "gene", filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures exonsBy
#' @export

setMethod("exonsBy", "src_organism",
function(x, by = c("tx", "gene"), filter = NULL) {
    by <- match.arg(by)
    exons <- exons(x, filter = filter)

    if (by == "tx") {
        table <- inner_join(exons, tbl(x, "ranges_tx"), by = "tx_id") 
        fields <- unique(
            c("tx_id", "exon_id", "exon_chrom", "exon_start", "exon_end",  
              "exon_strand", "exon_name", "exon_rank", names(filter)))
    }
    else if (by == "gene") {
        table <- inner_join(exons, tbl(x, "ranges_gene"))
        fields <- unique(
            c("entrez", "exon_id", "exon_chrom", "exon_start", "exon_end",  
              "exon_strand", "exon_name", names(filter)))
    }
    do.call(select_, c(list(table), as.list(fields)))
})


#' @examples
#' cdsBy(organism, by = "gene", filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures cdsBy
#' @export

setMethod("cdsBy", "src_organism",
function(x, by = c("tx", "gene"), filter = NULL) {
    by <- match.arg(by)
    cds <- cds(x, filter = filter)

    if (by == "tx") {
        table <- inner_join(cds, tbl(x, "ranges_tx"), by = "tx_id") 
        fields <- unique(
            c("tx_id", "cds_id", "cds_chrom", "cds_start", "cds_end",  
              "cds_strand", "cds_name", "exon_rank", names(filter)))
    }
    else if (by == "gene") {
        table <- inner_join(cds, tbl(x, "ranges_gene"))
        fields <- unique(
            c("entrez", "cds_id", "cds_chrom", "cds_start", "cds_end", 
              "cds_strand", "cds_name", names(filter)))
    }
    do.call(select_, c(list(table), as.list(fields)))
})


setMethod("intronsByTranscript", "src_organism",
function(x, filter=NULL) {
    tx <- transcripts(x, filter=filter)
    exn <- exonsBy(x, filter=filter)

    tx_gr <- tx %>% select(tx_id, tx_chrom, tx_start, tx_end, tx_strand) %>%
        collect(n=Inf) %>% as("GRanges")
    exn_grl <- exn %>%
        select(tx_id, exon_id, exon_chrom, exon_start, exon_end,
               exon_strand) %>%
        collect(n=Inf) %>% as("GRanges") %>% split(.$tx_id)
    
    tx_gr<- tx_gr[match(names(exn_grl), mcols(tx_gr)[, "tx_id"])]
    ans <- unlist(psetdiff(tx_gr, exn_grl))

    mcols(ans)[, "tx_id"] <- names(ans)
    unname(ans) %>% as.data.frame %>% tbl_df %>%
        select(tx_id, intron_chrom=seqnames, intron_start=start,
               intron_end=end, intron_strand=strand)
})
