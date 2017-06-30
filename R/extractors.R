#' @importFrom AnnotationFilter AnnotationFilter AnnotationFilterList field

.tbl_join <- function(x, table, tbls, filter) {
    if (is.null(filter))
        return(table)

	stopifnot(.check_filters(filter))

    fields <- .fields(filter)
    if ("granges" %in% fields)
        filter <- .filter_subset(filter, !(fields %in% "granges"))

    for (i in filter)
        stopifnot(is(i, "AnnotationFilter"))

    fields <- .fields(filter)

    if (!all(fields %in% columns(x)))
        stop(paste0("'", .filter_subset(filter, !(fields %in% columns(x))), "'",
                    collapse = ","), " filter not available")

    ## filter by fields from main table
    keep <- unique(fields[fields %in% colnames(table)])
    if (length(keep) != 0) {
        filters <- .filter(filter, keep)
        table <- table %>% filter_(filters)
        filter <- .filter_subset(filter, setdiff(fields, keep))
    }

    ## filter by fields from other tables
    fields <- .filter_names(filter)
    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || length(keep) == 0)
            next

        filters <- .filter(filter, keep)

        table <- inner_join(table, tbl(x, i)) %>% filter_(filters)
        fields <- setdiff(fields, keep)
    }
    table
}

.check_filters <- function(filter) {
	filters <- vapply(lapply(value(filter), function(x) class(x)[1]), function(x) x, character(1))
	!any(!(filters %in% supportedFilters()))
}

.keep <- function(filter, fields, fields_remove) {
    drop <- setdiff(fields_remove, .filter_names(filter))
    keep <- fields[!(fields %in% drop)]
}

.filter_subset <- function(filter, fields_subset) {
	if (is.null(filter))
		return(NULL)
	if (is.character(fields_subset))
		fields_subset <- .fields(filter) %in% fields_subset
	if (is.numeric(fields_subset))
		fields_subset <- seq_len(length(filter)) %in% fields_subset
	res <- value(filter)
	names(res) <- .fields(filter)
	ops <- .logicOp_subset(logicOp(filter), fields_subset)
	do.call(AnnotationFilterList, c(res[fields_subset], list(logicOp=ops)))
}

.logicOp_subset <- function(op, fields_subset) {
	keepOp <- rep(TRUE, length(op))
	for (i in seq_len(length(fields_subset))) {
		if (!fields_subset[i]) {
			first = i-1
			second = i
			if (first == 0 | second == length(fields_subset)) {
				if (first == 0)
					keepOp[second] <- FALSE
				else
					keepOp[first] <- FALSE
			} else {
				if (((op[first] == '&') & (op[second] == '|')) | 
						((op[first] == '|') & (op[second] == '&'))) {
					if (op[first] == '&')
						keepOp[second] <- FALSE
					else
						keepOp[first] <- FALSE
				} else {
					keepOp[second] <- FALSE
				}
			}
		}	
	}
	if (!any(fields_subset))
		return(character())
	if (table(fields_subset)["TRUE"] <= 1)
		character()
	else
		op[keepOp]
}

.filter_list <- function(filter) {
    if (!is.null(filter) & is(filter, "AnnotationFilter"))
        return(AnnotationFilterList(filter))
	if (!is.null(filter) & is(filter, "list") &
			!is(filter, "AnnotationFilterList"))
		return(do.call(AnnotationFilterList, filter))
	filter
}

.filter_names <- function(filter) {
    setdiff(.fields(filter), "granges")
}

.filter <- function(filter, keep) {
	fields <- .fields(filter)
	ops <- logicOp(filter)
    filter <- lapply(.filter_subset(filter, fields %in% keep), .convertFilter)
    #paste0(filter, collapse=" & ")
	paste(c(sprintf("%s %s ", head(filter, -1), ops), tail(filter, 1)), collapse="")
}

.return_tbl <- function(table, filter) {
    filter <- .filter_list(filter)
    if ("granges" %in% .fields(filter))
        message("filter by 'granges' only supported by methods returning ",
                "GRanges or GRangesList")
    table
}

.updateSeqinfo <- function(x, gr) {
    seqinfo <- .getSeqinfo(x)
    seqlevels(gr) <- seqlevels(seqinfo)
    seqinfo(gr) <- seqinfo
    gr
}

#' @importFrom IRanges subsetByOverlaps
.toGRanges <- function(x, table, filter) {
    if (!is.null(filter)) {
        filter <- .filter_list(filter)
        fields <- .fields(filter)
        condition <- endsWith(fields, "start") | endsWith(fields, "end")
        if (any(condition) && !fields[condition] %in% table)
            stop("use GRanges as filter instead of ", fields[condition])
    }

    gr <- table %>% collect(n=Inf) %>% as("GRanges")
    if ("granges" %in% .fields(filter)) {
        gr <- subsetByOverlaps(gr, .value(.filter_subset(filter, "granges")[[1]]))
    }
    .updateSeqinfo(x, gr)
}

.transcripts <- function(x, filter = NULL) {
    table <- tbl(x, "ranges_tx")
    tbls <- setdiff(src_tbls(x), "ranges_tx")
    table <- .tbl_join(x, table, tbls, filter)
}

.transcripts_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    table <- .transcripts(x, filter)
    fields <- unique(c(
        "tx_chrom", "tx_start", "tx_end", "tx_strand",
        "tx_id", "tx_name", .filter_names(filter)))
    do.call(select_, c(list(table), as.list(fields))) %>%
        arrange_(~ tx_id)
}


#' Extract genomic features from src_organism objects
#'
#' Generic functions to extract genomic features from an object. This
#' page documents the methods for \code{\link{src_organism}} objects
#' only.
#'
#' These are the main functions for extracting transcript information
#' from a \code{\link{src_organism}} object, inherited from
#' \code{\link[GenomicFeatures]{transcripts}} in
#' \code{GenomicFeatures} package. Two versions of results are
#' provided: \code{\link[tibble]{tibble}} (\code{transcripts_tbl()})
#' and \code{\link{GRanges}} or \code{\link{GRangesList}}
#' (\code{transcripts()}).
#'
#' @param x A \code{src_organism} object
#'
#' @param filter Either NULL or a named list of vectors to be used to
#'     restrict the output. Filter can also be a \code{\link{GRanges}} object
#'     using "GRangesFilter" (see examples).
#'
#' @return functions with \code{_tbl} return a \code{\link[tibble]{tibble}}
#'     object, other methods return a \code{\link{GRanges}} or
#'     \code{\link{GRangesList}} object.
#'
#' @seealso \code{\link{src_organism}} for creating a \code{src_organism}
#'     object.
#'
#' @author Yubo Cheng.
#'
#' @rdname extractors
#'
#' @examples
#' \dontrun{src <- src_ucsc("human")}
#' src <- src_organism(dbpath=hg38light())
#'
#' ## transcript coordinates with filter in tibble format
#' filters <- list(SymbolFilter(c("A1BG", "CDH2")))
#' transcripts_tbl(src, filters)
#'
#' transcripts_tbl(src, list(SymbolFilter("SNORD", "startsWith")))
#' transcripts_tbl(src, list(GoFilter("GO:0005615")))
#' transcripts_tbl(src, filter=list(
#'      SymbolFilter("SNORD", "startsWith"),
#'      TxStartFilter(25070000, "<")
#' ))
#' @export
transcripts_tbl <- function(x, filter = NULL) {
    .return_tbl(.transcripts_tbl(x, filter), filter)
}


#' @importFrom GenomicRanges GRanges
#'
#' @examples
#' ## transcript coordinates with filter in granges format
#' filters <- list(GRangesFilter(GenomicRanges::GRanges("chr15:1-25070000")))
#' transcripts(src, filters)
#'
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomeInfoDb seqinfo<- seqlevels seqlevels<-
#' @rdname extractors
#' @exportMethod transcripts
setMethod("transcripts", "src_organism",
    function(x, filter = NULL)
{
    .toGRanges(x, .transcripts_tbl(x, filter), filter)
})


.exons <- function(x, filter = NULL) {
    table <- tbl(x, "ranges_exon")
    tbls <- setdiff(src_tbls(x), "ranges_exon")
    table <- .tbl_join(x, table, tbls, filter)
}

.exons_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    fields_remove <- c(x$schema, "tx_id", "exon_name", "exon_rank")
    fields <- unique(c(
        "exon_chrom", "exon_start", "exon_end", "exon_strand",
        x$schema, "tx_id", "exon_id", "exon_name", "exon_rank",
        .filter_names(filter)))
    keep <- .keep(filter, fields, fields_remove)

    table <- .exons(x, filter)
    do.call(select_, c(list(table), as.list(keep))) %>% arrange_(~ exon_id)
}

#' @rdname extractors
#' @export
exons_tbl <- function(x, filter = NULL) {
    .return_tbl(.exons_tbl(x, filter), filter)
}


#' @rdname extractors
#' @importFrom GenomicFeatures exons
#' @exportMethod exons
setMethod("exons", "src_organism",
    function(x, filter = NULL)
{
    .toGRanges(x, .exons_tbl(x, filter), filter)
})


.cds <- function(x, filter = NULL) {
    table <- tbl(x, "ranges_cds")
    tbls <- setdiff(src_tbls(x), "ranges_cds")
    table <- .tbl_join(x, table, tbls, filter)
}

.cds_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    fields_remove <- c(x$schema, "tx_id", "cds_name", "exon_rank")
    fields <- unique(c(
        "cds_chrom", "cds_start", "cds_end", "cds_strand",
        x$schema, "tx_id", "cds_id", "cds_name", "exon_rank",
        .filter_names(filter)))
    keep <- .keep(filter, fields, fields_remove)

    table <- .cds(x, filter)
    do.call(select_, c(list(table), as.list(keep))) %>% arrange_(~ cds_id)
}

#' @rdname extractors
#' @export
cds_tbl <- function(x, filter = NULL) {
    .return_tbl(.cds_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures cds
#' @exportMethod cds
setMethod("cds", "src_organism",
    function(x, filter = NULL)
{
    .toGRanges(x, .cds_tbl(x, filter), filter)
})

.genes_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    table <- tbl(x, "ranges_gene")
    # genes() from GenomicFeatures
    # table <- dbGetQuery(
    #     x$con,
    #     "SELECT * FROM ranges_gene
    #      WHERE entrez NOT IN
    #          (SELECT entrez FROM ranges_gene
    #           GROUP BY entrez
    #           HAVING count(*) > 1)")
    tbls <- setdiff(src_tbls(x), "ranges_gene")
    table <- .tbl_join(x, table, tbls, filter)
    fields <- unique(c(
        "gene_chrom", "gene_start", "gene_end", "gene_strand",
        x$schema, .filter_names(filter)))
#' @rdname filter
#' @export
GRangesFilter <- function(value) {
    new("GRangesFilter",
        field="granges",
        value=value)
}

setValidity("GRangesFilter", function(object) {
    value <- .value(object)
    txt <- character()
    if (!is(value, "GRanges"))
        txt <- c(txt, "'value' must be 'GRanges' object")
    if (length(txt)) txt else TRUE
})

#' @rdname filter
#' @exportMethod show
setMethod("show", "GRangesFilter",
    function(object)
{
    cat("class:", class(object),
        "\nvalue:\n")
    print(.value(object))
})
    do.call(select_, c(list(table), as.list(fields))) %>% arrange_(x$schema)
}

#' @rdname extractors
#' @export
genes_tbl <- function(x, filter = NULL) {
    .return_tbl(.genes_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures genes
#' @exportMethod genes
setMethod("genes", "src_organism",
    function(x, filter = NULL)
{
    .toGRanges(x, .genes_tbl(x, filter), filter)
})

.promoters_tbl <- function(x, upstream, downstream, filter = NULL) {
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
    if(missing(upstream))
        upstream = 0
    if(missing(downstream))
        downstream = 2200

    filter <- .filter_list(filter)
    table <- .transcripts_tbl(x, filter = filter) %>%
        mutate_(start = ~ ifelse(tx_strand == "-",
                              tx_end - downstream + 1,
                              tx_start - upstream),
               end = ~ ifelse(tx_strand == "-",
                            tx_end + upstream,
                            tx_start + downstream - 1)) %>% collect(n=Inf)

    table <- rename_(table, chrom = ~ tx_chrom)
    table <- rename_(table, strand = ~ tx_strand)

    fields <- unique(
        c("chrom", "start", "end", "strand", "tx_id", "tx_name",
          .filter_names(filter)))
    do.call(select_, c(list(table), as.list(fields))) %>% arrange_(~ tx_id)
}

#' @param upstream For \code{promoters()}: An integer(1) value indicating
#' the number of bases upstream from the transcription start site.
#'
#' @param downstream For \code{promoters()}: An integer(1) value indicating
#' the number of bases downstream from the transcription start site.
#'
#' @rdname extractors
#' @export
promoters_tbl <- function(x, upstream, downstream, filter = NULL) {
    .return_tbl(.promoters_tbl(x, upstream, downstream, filter), filter)
}

#' @examples
#' ## promoters
#' promoters(src, upstream=100, downstream=50,
#'           filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures promoters
#' @importFrom S4Vectors isSingleNumber
#' @exportMethod promoters
setMethod("promoters", "src_organism",
    function(x, upstream, downstream, filter = NULL)
{
    .toGRanges(x, .promoters_tbl(x, upstream, downstream, filter), filter)
})


.toGRangesList <- function(x, type, by, filter = NULL) {
    fun <- switch (type, transcripts = .transcriptsBy_tbl,
                   exons = .exonsBy_tbl,
                   cds = .cdsBy_tbl,
                   fiveUTRs = .fiveUTRsByTranscript_tbl,
                   threeUTRs = .threeUTRsByTranscript_tbl)

    if (type %in% c("transcripts", "exons", "cds"))
        table <- .toGRanges(x, fun(x, by, filter), filter)
    else
        table <- .toGRanges(x, fun(x, filter), filter)

    f <- switch(by, gene=x$schema, exon="exon_id", cds="cds_id", tx="tx_id")
    grp <- mcols(table)[[f]]
    mcols(table)[[f]] <- NULL
    split(table, grp)
}


.transcriptsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    tx <- .transcripts(x, filter = filter)

    if (by == "gene") {
        table <- inner_join(tx, tbl(x, "ranges_gene"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              x$schema, "tx_id", "tx_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema)
    }
    else if (by == "exon") {
        table <- left_join(tx, tbl(x, "ranges_exon"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              "exon_id", "tx_id", "tx_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ exon_id)
    }
    else if (by == "cds") {
        table <- left_join(tx, tbl(x, "ranges_cds"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              "cds_id", "tx_id", "tx_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ cds_id)
    }
}

#' @rdname extractors
#' @export
transcriptsBy_tbl <- function(x, by = c("gene", "exon", "cds"), filter = NULL)
{
    by <- match.arg(by)
    .return_tbl(.transcriptsBy_tbl(x, by, filter), filter)
}

#' @param by One of "gene", "exon", "cds" or "tx". Determines the grouping.
#'
#' @examples
#' ## transcriptsBy
#' transcriptsBy(src, by = "exon", filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures transcriptsBy
#' @exportMethod transcriptsBy
setMethod("transcriptsBy", "src_organism",
    function(x, by = c("gene", "exon", "cds"), filter = NULL)
{
    by <- match.arg(by)
    .toGRangesList(x, "transcripts", by, filter)
})

.exonsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    exons <- .exons(x, filter = filter)

    if (by == "tx") {
        table <- left_join(exons, tbl(x, "ranges_tx"))
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",
              "tx_id", "exon_id", "exon_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ tx_id)
    }
    else if (by == "gene") {
        table <- inner_join(exons, tbl(x, "ranges_gene"))
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",
              x$schema, "exon_id", "exon_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema)
    }
}

#' @rdname extractors
#' @export
exonsBy_tbl <- function(x, by = c("tx", "gene"), filter = NULL) {
    by <- match.arg(by)
    .return_tbl(.exonsBy_tbl(x, by, filter), filter)
}

#' @examples
#' ## exonsBy
#' exonsBy(src, filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures exonsBy
#' @exportMethod exonsBy
setMethod("exonsBy", "src_organism",
    function(x, by = c("tx", "gene"), filter = NULL)
{
    by <- match.arg(by)
    .toGRangesList(x, "exons", by, filter)
})

.cdsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    cds <- .cds(x, filter = filter)

    if (by == "tx") {
        table <- left_join(cds, tbl(x, "ranges_tx"))
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand",
              "tx_id", "cds_id", "cds_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ tx_id)
    }
    else if (by == "gene") {
        table <- inner_join(cds, tbl(x, "ranges_gene"))
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand",
              x$schema, "cds_id", "cds_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema)
    }
}

#' @rdname extractors
#' @export
cdsBy_tbl <- function(x, by = c("tx", "gene"), filter = NULL) {
    by <- match.arg(by)
    .return_tbl(.cdsBy_tbl(x, by, filter), filter)
}


#' @rdname extractors
#' @importFrom GenomicFeatures cdsBy
#' @exportMethod cdsBy
setMethod("cdsBy", "src_organism",
    function(x, by = c("tx", "gene"), filter = NULL)
{
    by <- match.arg(by)
    .toGRangesList(x, "cds", by, filter)
})


#' @rdname extractors
#' @export
intronsByTranscript_tbl <-
    function(x, filter = NULL)
{
    ans <- unlist(intronsByTranscript(x, filter))
    mcols(ans)[, "tx_id"] <- names(ans)
    unname(unlist(ans)) %>% as.data.frame %>% tbl_df %>%
        dplyr::select_(.dots = c('tx_id',
                          intron_chrom = 'seqnames',
                          intron_start = 'start',
                          intron_end = 'end',
                          intron_strand = 'strand'))
}

#' @examples
#' ## intronsByTranscript
#' intronsByTranscript(src, filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom GenomicRanges split mcols mcols<-
#' @importFrom IRanges psetdiff
#' @exportMethod intronsByTranscript
setMethod("intronsByTranscript", "src_organism",
    function(x, filter=NULL)
{
    filter <- .filter_list(filter)
    tx <- .transcripts(x, filter=filter)
    exn <- .exonsBy_tbl(x, by="tx", filter=filter)

    table <- tx %>% collect(n=Inf)

    fields <- unique(
        c("tx_id", "tx_chrom", "tx_start", "tx_end", "tx_strand",
          .filter_names(filter)))
    tx_gr <- do.call(select_, c(list(table), as.list(fields))) %>%
        as("GRanges")

    tx_gr <- .updateSeqinfo(x, tx_gr)

    table <- exn %>% collect(n=Inf)

    fields <- unique(
        c("exon_chrom", "exon_start", "exon_end", "exon_strand", "tx_id",
          "exon_id", .filter_names(filter)))
    exn_grl <- do.call(select_, c(list(table), as.list(fields))) %>%
        as("GRanges")
    exn_grl <- split(exn_grl, exn_grl$tx_id)

    tx_gr<- tx_gr[match(names(exn_grl), mcols(tx_gr)[, "tx_id"])]
    ans <- psetdiff(tx_gr, exn_grl)
    ans
})


.getSplicings <- function(x, filter=NULL) {
    exon <- .exons(x, filter=filter)
    cds <- .cds(x, filter=filter)

    exon_txid <- exon %>% dplyr::select_(~ tx_id) %>% collect(n = Inf)
    exon_txid <- exon_txid[["tx_id"]]
    cds_txid <- cds %>% dplyr::select_(~ tx_id) %>% collect(n = Inf)
    cds_txid <- cds_txid[["tx_id"]]
    exclude <- setdiff(exon_txid, cds_txid)

    table <- left_join(exon, cds,
                       by = c("tx_id", "exon_rank", .filter_names(filter)))

    fields <- unique(
        c("tx_id", "exon_rank", "exon_id", "exon_name", "exon_chrom",
          "exon_strand", "exon_start", "exon_end", "cds_id", "cds_start",
          "cds_end", .filter_names(filter)))
    splicings <- do.call(select_, c(list(table), as.list(fields))) %>%
        collect(n=Inf)

    if(length(exclude) != 0)
        splicings <- splicings %>% filter_(~ !tx_id %in% exclude)
    splicings
}

.UTRsByTranscript <- function(x, filter = NULL, strand1, strand2) {
    splicings <- .getSplicings(x, filter)

    exons_with_cds <- which(!is.na(splicings$cds_id))
    ifelse (strand1 == "-",
            idx <- GenomicFeatures:::.exons_with_5utr
                    (splicings$tx_id, exons_with_cds),
            idx <- GenomicFeatures:::.exons_with_3utr
                    (splicings$tx_id, exons_with_cds))

    splicings <- S4Vectors:::extract_data_frame_rows(splicings, idx)

    table <- splicings %>%
        mutate_(start = ~
                   ifelse(!is.na(cds_id) & exon_strand == strand1,
                          cds_end + 1L,
                          exon_start),
               end = ~
                   ifelse(!is.na(cds_id) & exon_strand == strand2,
                          cds_start - 1L,
                          exon_end)) %>%
        filter_(~ start <= end) %>% collect(n=Inf) %>% tbl_df

    table <- rename_(table, chrom = ~ exon_chrom)
    table <- rename_(table, strand = ~ exon_strand)

    fields <- unique(
        c("chrom", "start", "end", "strand", "tx_id", "exon_id", "exon_name",
          "exon_rank", .filter_names(filter)))
    do.call(select_, c(list(table), as.list(fields))) %>%
        arrange_(~ tx_id, ~ exon_rank)
}

.fiveUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    .UTRsByTranscript(x, filter, "-", "+")
}

#' @rdname extractors
#' @export
fiveUTRsByTranscript_tbl <- function(x, filter = NULL) {
    .return_tbl(.fiveUTRsByTranscript_tbl(x, filter), filter)
}

#' @examples
#' ## fiveUTRsByTranscript
#' fiveUTRsByTranscript(src, filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @exportMethod fiveUTRsByTranscript
setMethod("fiveUTRsByTranscript", "src_organism",
    function(x, filter=NULL)
{
    .toGRangesList(x, "fiveUTRs", "tx", filter)
})

.threeUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    .UTRsByTranscript(x, filter, "+", "-")
}

#' @rdname extractors
#' @export
threeUTRsByTranscript_tbl <- function(x, filter = NULL) {
    .return_tbl(.threeUTRsByTranscript_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @exportMethod threeUTRsByTranscript
setMethod("threeUTRsByTranscript", "src_organism",
    function(x, filter=NULL)
{
    .toGRangesList(x, "threeUTRs", "tx", filter)
})
