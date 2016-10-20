.tbl_filter <- function(keep1, table, filter) {
    values <- paste0("'", filter[[keep1]], "'", collapse=", ")
    op <- if (length(filter[[keep1]]) == 1) "==" else "%in%"
    sprintf("%s %s c(%s)", keep1, op, values)
}

.tbl_join <- function(x, table, filter, schema) {
    if (is.null(filter))
        return(table)

    if (schema %in% names(filter)) {
        filters <- .tbl_filter(schema, table, filter)
        table <- table %>% filter_(filters)
        filter <- filter[names(filter) != schema]
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
    schema <- attr(x, "schema")
    table <- .tbl_join(x, table, filter, schema) %>% arrange(tx_id)
    fields <- unique(
        c("tx_chrom", "tx_start", "tx_end", "tx_strand", 
          schema, "tx_id", "tx_name", names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})

#' @rdname src_organism
#' @importFrom GenomicFeatures exons
#' @export

setMethod("exons", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_exon")
    schema <- attr(x, "schema")
    table <- .tbl_join(x, table, filter, schema) %>% arrange(exon_id)
    fields <- unique(
        c("exon_chrom", "exon_start", "exon_end", "exon_strand",
          schema, "tx_id", "exon_id", "exon_name", "exon_rank",
          names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})


#' @rdname src_organism
#' @importFrom GenomicFeatures cds
#' @export

setMethod("cds", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_cds")
    schema <- attr(x, "schema")
    table <- .tbl_join(x, table, filter, schema) %>% arrange(cds_id)
    fields <- unique(
        c("cds_chrom", "cds_start", "cds_end", "cds_strand", 
          schema, "tx_id", "cds_id", "cds_name", "exon_rank", names(filter)))
    do.call(select_, c(list(table), as.list(fields)))
})


#' @rdname src_organism
#' @importFrom GenomicFeatures genes
#' @export

setMethod("genes", "src_organism", function(x, filter = NULL) {
    table <- tbl(x, "ranges_gene")
    schema <- attr(x, "schema")
    table <- .tbl_join(x, table, filter, schema) %>% arrange_(schema)
    fields <- unique(
        c("tx_chrom", "gene_start", "gene_end", "tx_strand", 
          schema, names(filter)))
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
#' @importFrom S4Vectors isSingleNumber
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
                         tx_start + downstream - 1)) %>%
        arrange(tx_id)
    
    schema <- attr(x, "schema")
    fields <- unique(
        c("tx_chrom", "start", "end", "tx_strand", 
          schema, "tx_id", "tx_name", names(filter)))
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
    schema <- attr(x, "schema")

    if (by == "gene") {
        table <- inner_join(tx, tbl(x, "ranges_gene")) 
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",  
              schema, "tx_id", "tx_name", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange_(schema)
    }
    else if (by == "exon") {
        table <- inner_join(tx, tbl(x, "ranges_exon"), by = "tx_id") 
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",  
              "exon_id", "tx_id", "tx_name", "exon_rank", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange(exon_id)
    }
    else if (by == "cds") {
        table <- inner_join(tx, tbl(x, "ranges_cds"), by = "tx_id") 
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",  
              "cds_id", "tx_id", "tx_name", "exon_rank", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange(cds_id)
    }
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
    schema <- attr(x, "schema")

    if (by == "tx") {
        table <- inner_join(exons, tbl(x, "ranges_tx"), by = "tx_id") 
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",  
              "tx_id", "exon_id", "exon_name", "exon_rank", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange(tx_id)
    }
    else if (by == "gene") {
        table <- inner_join(exons, tbl(x, "ranges_gene")) 
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",  
              schema, "exon_id", "exon_name", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange_(schema)
    }
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
    schema <- attr(x, "schema")

    if (by == "tx") {
        table <- inner_join(cds, tbl(x, "ranges_tx"), by = "tx_id") 
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand",  
              "tx_id", "cds_id", "cds_name", "exon_rank", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange(tx_id) 
    }
    else if (by == "gene") {
        table <- inner_join(cds, tbl(x, "ranges_gene")) 
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand", 
              schema, "cds_id", "cds_name", names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>% arrange_(schema)
    }
})


#' @examples
#' intronsByTranscript(organism, filter = list(symbol="PTEN"))
#'
#' @rdname src_organism
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom GenomicRanges split mcols mcols<-
#' @importFrom IRanges psetdiff
#' @export
setMethod("intronsByTranscript", "src_organism",
function(x, filter=NULL) {
    tx <- transcripts(x, filter=filter)
    exn <- exonsBy(x, filter=filter)

    tx_gr <- tx %>% 
        dplyr::select(tx_id, tx_chrom, tx_start, tx_end, tx_strand) %>%
        collect(n=Inf) %>% as("GRanges")
    exn_grl <- exn %>% collect(n=Inf) %>% 
        dplyr::select(tx_id, exon_id, exon_chrom, exon_start, exon_end,
               exon_strand) %>%
        as("GRanges") %>% split(.$tx_id)
    
    tx_gr<- tx_gr[match(names(exn_grl), mcols(tx_gr)[, "tx_id"])]
    ans <- unlist(psetdiff(tx_gr, exn_grl))

    mcols(ans)[, "tx_id"] <- names(ans)
    unname(ans) %>% as.data.frame %>% tbl_df %>%
        dplyr::select(tx_id, intron_chrom=seqnames, intron_start=start,
               intron_end=end, intron_strand=strand)
})


.getSplicings <- function(x, filter=NULL) {
    exon <- exons(x, filter=filter)
    cds <- cds(x, filter=filter)
    
    inner_join(exon, cds, by = c("tx_id", "exon_rank")) %>%
        dplyr::select(tx_id, exon_rank, exon_id, exon_name, exon_chrom, 
                      exon_strand, exon_start, exon_end, cds_id, cds_start,  
                      cds_end) %>% 
        collect(n=Inf) %>%
        arrange(tx_id, exon_rank)
}

#' @examples
#' fiveUTRsByTranscript(organism, filter = list(symbol="PTEN"))
#' 
#' @rdname src_organism
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @export
setMethod("fiveUTRsByTranscript", "src_organism",function(x, filter=NULL) {
    splicings <- .getSplicings(x, filter)
    
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- GenomicFeatures:::.exons_with_5utr(splicings$tx_id, exons_with_cds)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, idx)
    
    class(splicings) <- c("tbl_df", class(splicings))
    
    splicings %>% 
        mutate(start = 
                   ifelse(!is.na(cds_id) & exon_strand == "-",
                        cds_end + 1L, 
                        exon_start), 
               end = 
                   ifelse(!is.na(cds_id) & exon_strand == "+", 
                        cds_start - 1L, 
                        exon_end)) %>% 
        filter(start <= end) %>%
        dplyr::select(exon_chrom, start, end, exon_strand, 
                      tx_id, exon_id, exon_name, exon_rank)
})


#' @examples
#' threeUTRsByTranscript(organism, filter = list(symbol="PTEN"))
#' 
#' @rdname src_organism
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @export
setMethod("threeUTRsByTranscript", "src_organism",function(x, filter=NULL) {
    splicings <- .getSplicings(x, filter)
    
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- GenomicFeatures:::.exons_with_3utr(splicings$tx_id, exons_with_cds)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, idx)
    
    class(splicings) <- c("tbl_df", class(splicings))
    
    splicings %>% 
        mutate(start = 
                   ifelse(!is.na(cds_id) & exon_strand == "+",
                          cds_end + 1L, 
                          exon_start), 
               end = 
                   ifelse(!is.na(cds_id) & exon_strand == "-", 
                          cds_start - 1L, 
                          exon_end)) %>% 
        filter(start <= end) %>%
        dplyr::select(exon_chrom, start, end, exon_strand, 
                      tx_id, exon_id, exon_name, exon_rank)
})
