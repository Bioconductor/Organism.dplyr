#' @importFrom AnnotationFilter AnnotationFilter AnnotationFilterList field
#'      value convertFilter
#' @importFrom AnnotationDbi columns
#' @importFrom dplyr %>% as_tibble inner_join full_join filter_
.tbl_join <- function(x, filter, main_table, table_names) {
    if (is.null(filter))
        return(main_table)

    if (!.check_filters(x, filter))
        stop("Some filters are not avaiable.")

    res <- lapply(filter, function(i) {
        if (is(i, "AnnotationFilterList"))
            .tbl_join(x, i, main_table, table_names)
        else {
            f_field <- field(i)
                if (f_field == "granges") {
                    gr <- main_table %>% collect(n=Inf) %>% as("GRanges")
                    if (!is.null(i)) {
                        gr <- subsetByOverlaps(
                            gr,
                            value(i),
                            type=condition(i)
                        )
                    }
                    #.updateSeqinfo(x, gr)
                    gr <- as_tibble(gr)
                    gr <- gr[, setdiff(colnames(gr), "width")]
                    colnames(gr)[1:4] <- colnames(main_table)[c(1, 3, 4, 2)]
                    name <- .getNewTableName()
                    dbWriteTable(x$db, name, gr)
                    tbl(x$db, name)
                }
                else {
                    dplyr_filter <- convertFilter(i)
                    selected <-
                        vapply(table_names, function(j) f_field %in% j, logical(1))
                    selected_table <- names(selected)[selected]
                    main_table <- .iterTable(x, main_table)
                    val <- tbl(x$db, selected_table[1]) %>% filter_(dplyr_filter)
                    .getAllTableValues(x, val, table_names, selected)
                }
            }
        }
    )
    ops <- logicOp(filter)
    table <- .iterTable(x, main_table)
    main_table <- .innerJoin(main_table, res[[1]])
    res <- res[-1]
    for (i in seq_along(res)) {
        main_table <- .iterTable(x, main_table)
        if(ops[[i]] == '&')
            main_table <- .innerJoin(main_table, res[[i]])
        else
            main_table <- .fullJoin(main_table, res[[i]])
    }
    main_table
}

.check_filters <- function(x, filter) {
    fields <- c(as.character(.supportedFilters()[,2]))
    filters <- vapply(value(filter), function(i) {
        is(i, "AnnotationFilterList") || field(i) %in% fields
        },
        logical(1)
    )
    all(filters)
}

.getAllTableValues <- function(x, table, table_names, selected) {
    table_names <- table_names[!selected]
    for (i in names(table_names)) {
        tmp <- tbl(x$db, i)
        table <- .iterTable(x, table)
        table <- left_join(table, tmp)
    }
    table
}

.filter_list <- function(filter) {
    if (!is.null(filter) && is(filter, "AnnotationFilter"))
        return(AnnotationFilterList(filter))
    if (!is.null(filter) && is(filter, "list") &
        !is(filter, "AnnotationFilterList"))
        return(do.call(AnnotationFilterList, filter))
    filter
}

.filter_names <- function(filter) {
    setdiff(.fields(filter), "granges")
}

.tableNames <- function(x, filter, main_ranges) {
    all_select_values <- unique(.fields(filter))
    tbls <- setdiff(src_tbls(x), main_ranges)
    names(tbls) <- tbls
    colnames <- lapply(tbls, function(tbl) colnames(tbl(x, tbl, .load_tbl_only=TRUE)))
    mapped <- Map(
        function(nms, value) nms[nms %in% value],
        colnames,
        MoreArgs=list(c(x$schema, all_select_values))
    )
    mapped <- mapped[lengths(mapped) > 1]
    main_map <- list(colnames(tbl(x, main_ranges, .load_tbl_only=TRUE)))
    names(main_map) <- main_ranges
    c(main_map, mapped)
}

.return_tbl <- function(table, filter) {
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
    ## Filter out any rows that contain NA in chrom, start, end, or strand
    table <- table %>% filter_at(vars(c(1, 2, 3, 4)), all_vars(!is.na(.)))

    gr <- table %>% collect(n=Inf) %>% as("GRanges")
#    .updateSeqinfo(x, gr)
}

.checkCompatibleStartEnds <- function(type, filter) {
    fields <- .fields(filter)
    starts <- fields[grep("start", fields)]
    ends <- fields[grep("end", fields)]
    starts <- grepl(type, starts)
    ends <- grepl(type, ends)
    all(starts) && all(ends)
}

.checkDepth <- function(table) {
    max <- 85
    output <- capture.output(table %>% show_query(), type="message")
    output <- tail(output, -1)
    max < length(output)
}

.iterTable <- function(x, table) {
    if (.checkDepth(table)) {
        table <- table %>% collect(n=Inf)
        nam <- .getNewTableName()
        dbWriteTable(x$db, nam, table)
        res <- tbl(x$db, nam)
        #class(res) <- c("tbl_organism", class(res))
        res
    }
    else
        table
}

.getNewTableName <- local({
    id <- 0L
    function() {
        id <<- id + 1
        paste0("table", id)
    }
})

.xscripts <- function(x, main_ranges, filter = NULL) {
    table_names <- .tableNames(x, filter, main_ranges)
    add_tables <- setdiff(names(table_names), dbListTables(x$db))
    for (i in add_tables) {
        table <- tbl(x, i, .load_tbl_only=TRUE) %>% collect(n=Inf)
        dbWriteTable(x$db, i, table)
    }
    #db <- src_dbi(db)
    table <- tbl(x$db, main_ranges)
    table <- .tbl_join(x, filter, table, table_names)
    table
}

.transcripts_tbl <- function(x, filter = NULL) {
    table <- .xscripts(x, "ranges_tx", filter)
    fields <- unique(c(
        "tx_chrom", "tx_start", "tx_end", "tx_strand",
        "tx_id", "tx_name", .filter_names(filter)))
    table <- .iterTable(x, table)
    do.call(dplyr::select, c(list(table), as.list(fields))) %>%
        arrange_(~ tx_id) %>% distinct()
}

.parseFilterInput <- function(filter) {
    if (is.language(filter))
        filter <- AnnotationFilter(filter)
    filter <- .filter_list(filter)
    if (is(filter, "AnnotationFilterList"))
        filter <- distributeNegation(filter)
    filter
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
    filter <- .parseFilterInput(filter)
    .transcripts_tbl(x, filter)
}


#' @importFrom GenomicRanges GRanges
#'
#' @param granges A \code{GRangesFilter} object to subset the resulting
#'     \code{GRanges} object by.  Is \code{NULL} if missing.
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
    function(x, filter = NULL) {
        filter <- .parseFilterInput(filter)
        if(!.checkCompatibleStartEnds("tx", filter))
            stop("Only TxStartFilter and TxEndFilters can be used as start and 
                end filters for transcripts method")
        .toGRanges(x, .transcripts_tbl(x, filter), filter)
})

.exons_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    table <- .xscripts(x, "ranges_exon", filter)
    fields <- unique(c(
        "exon_chrom", "exon_start", "exon_end", "exon_strand",
        x$schema, "tx_id", "exon_id", "exon_name", "exon_rank",
        .filter_names(filter)))
    fields_remove <- c(x$schema, "tx_id", "exon_name", "exon_rank")
    fields <- setdiff(fields, fields_remove)
    table <- .iterTable(x, table)
    do.call(select_, c(list(table), as.list(fields))) %>%
        arrange_(~ exon_id) %>% distinct()
}

#' @rdname extractors
#' @export
exons_tbl <- function(x, filter = NULL) {
    filter <- .parseFilterInput(filter)
    .return_tbl(.exons_tbl(x, filter), filter)
}


#' @rdname extractors
#' @importFrom GenomicFeatures exons
#' @exportMethod exons
setMethod("exons", "src_organism", function(x, filter = NULL, granges = NULL) {
    filter <- .parseFilterInput(filter)
    if(!.checkCompatibleStartEnds("exon", filter))
        stop("Only ExonStartFilter and ExonEndFilters can be used as start and 
            end filters for exons method")
    .toGRanges(x, .exons_tbl(x, filter), filter, granges)
})

.cds_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    table <- .xscripts(x, "ranges_cds", filter)
    fields <- unique(c(
        "cds_chrom", "cds_start", "cds_end", "cds_strand",
        x$schema, "tx_id", "cds_id", "cds_name", "exon_rank",
        .filter_names(filter)))
    fields_remove <- c(x$schema, "tx_id", "cds_name", "exon_rank")
    fields <- setdiff(fields, fields_remove)
    table <- .iterTable(x, table)
    do.call(select_, c(list(table), as.list(fields))) %>% arrange_(~ cds_id) %>%
        distinct()
}

#' @rdname extractors
#' @export
cds_tbl <- function(x, filter = NULL) {
    filter <- .parseFilterInput(filter)
    .return_tbl(.cds_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures cds
#' @exportMethod cds
setMethod("cds", "src_organism", function(x, filter = NULL, granges = NULL) {
    filter <- .parseFilterInput(filter)
    if(!.checkCompatibleStartEnds("cds", filter))
        stop("Only CdsStartFilter and CdsEndFilters can be used as start and 
            end filters for cds method")
    .toGRanges(x, .cds_tbl(x, filter), filter, granges)
})

.genes_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    table <- .xscripts(x, "ranges_gene", filter)
    fields <- unique(c(
        "gene_chrom", "gene_start", "gene_end", "gene_strand",
        x$schema, .filter_names(filter)))
    table <- .iterTable(x, table)
    do.call(select_, c(list(table), as.list(fields))) %>% arrange_(x$schema) %>%
        distinct()
}

#' @rdname extractors
#' @export
genes_tbl <- function(x, filter = NULL) {
    filter <- .parseFilterInput(filter)
    .return_tbl(.genes_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures genes
#' @exportMethod genes
setMethod("genes", "src_organism", function(x, filter = NULL, granges = NULL) {
    filter <- .parseFilterInput(filter)
    if(!.checkCompatibleStartEnds("gene", filter))
        stop("Only GeneStartFilter and GeneEndFilters can be used as start and 
            end filters for genes method")
    .toGRanges(x, .genes_tbl(x, filter), filter, granges)
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
    table <- .iterTable(x, table)
    do.call(select_, c(list(table), as.list(fields))) %>% arrange_(~ tx_id) %>%
        distinct()
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
    filter <- .parseFilterInput(filter)
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
    function(x, upstream, downstream, filter = NULL, granges = NULL) {
        filter <- .parseFilterInput(filter)
        .toGRanges(x, .promoters_tbl(x, upstream, downstream, filter), filter,
            granges)
})


.toGRangesList <- function(x, type, by, filter = NULL, granges = NULL) {
    fun <- switch (type, transcripts = .transcriptsBy_tbl,
                   exons = .exonsBy_tbl,
                   cds = .cdsBy_tbl,
                   fiveUTRs = .fiveUTRsByTranscript_tbl,
                   threeUTRs = .threeUTRsByTranscript_tbl)

    if (type %in% c("transcripts", "exons", "cds"))
        table <- .toGRanges(x, fun(x, by, filter), filter, granges)
    else
        table <- .toGRanges(x, fun(x, filter), filter, granges)

    f <- switch(by, gene=x$schema, exon="exon_id", cds="cds_id", tx="tx_id")
    grp <- mcols(table)[[f]]
    mcols(table)[[f]] <- NULL
    split(table, grp)
}


.transcriptsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    tx <- .xscripts(x, "ranges_tx", filter = filter)

    if (by == "gene") {
        table <- inner_join(tx, tbl(x, "ranges_gene"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              x$schema, "tx_id", "tx_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema) %>% distinct()
    }
    else if (by == "exon") {
        table <- left_join(tx, tbl(x, "ranges_exon"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              "exon_id", "tx_id", "tx_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ exon_id) %>% distinct()
    }
    else if (by == "cds") {
        table <- left_join(tx, tbl(x, "ranges_cds"))
        fields <- unique(
            c("tx_chrom", "tx_start", "tx_end", "tx_strand",
              "cds_id", "tx_id", "tx_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ cds_id) %>% distinct()
    }
}

#' @rdname extractors
#' @export
transcriptsBy_tbl <- function(x, by = c("gene", "exon", "cds"), filter = NULL) {
    filter <- .parseFilterInput(filter)
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
    function(x, by = c("gene", "exon", "cds"), filter = NULL, granges = NULL) {
        filter <- .parseFilterInput(filter)
        by <- match.arg(by)
        .toGRangesList(x, "transcripts", by, filter, granges)
})

.exonsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    exons <- .xscripts(x, "ranges_exon", filter = filter)

    if (by == "tx") {
        table <- left_join(exons, tbl(x, "ranges_tx"))
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",
              "tx_id", "exon_id", "exon_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ tx_id) %>% distinct()
    }
    else if (by == "gene") {
        table <- inner_join(exons, tbl(x, "ranges_gene"))
        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand",
              x$schema, "exon_id", "exon_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema) %>% distinct()
    }
}

#' @rdname extractors
#' @export
exonsBy_tbl <- function(x, by = c("tx", "gene"), filter = NULL) {
    filter <- .parseFilterInput(filter)
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
    function(x, by = c("tx", "gene"), filter = NULL, granges = NULL) {
        filter <- .parseFilterInput(filter)
        by <- match.arg(by)
        .toGRangesList(x, "exons", by, filter, granges)
})

.cdsBy_tbl <- function(x, by, filter = NULL) {
    filter <- .filter_list(filter)
    cds <- .xscripts(x, "ranges_cds", filter = filter)

    if (by == "tx") {
        table <- left_join(cds, tbl(x, "ranges_tx"))
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand",
              "tx_id", "cds_id", "cds_name", "exon_rank",
              .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(~ tx_id) %>% distinct()
    }
    else if (by == "gene") {
        table <- inner_join(cds, tbl(x, "ranges_gene"))
        fields <- unique(
            c("cds_chrom", "cds_start", "cds_end", "cds_strand",
              x$schema, "cds_id", "cds_name", .filter_names(filter)))
        do.call(select_, c(list(table), as.list(fields))) %>%
            arrange_(x$schema) %>% distinct()
    }
}

#' @rdname extractors
#' @export
cdsBy_tbl <- function(x, by = c("tx", "gene"), filter = NULL) {
    filter <- .parseFilterInput(filter)
    by <- match.arg(by)
    .return_tbl(.cdsBy_tbl(x, by, filter), filter)
}


#' @rdname extractors
#' @importFrom GenomicFeatures cdsBy
#' @exportMethod cdsBy
setMethod("cdsBy", "src_organism",
    function(x, by = c("tx", "gene"), filter = NULL, granges = NULL) {
        filter <- .parseFilterInput(filter)
        by <- match.arg(by)
        .toGRangesList(x, "cds", by, filter, granges)
})


#' @rdname extractors
#' @export
intronsByTranscript_tbl <-
    function(x, filter = NULL) {
        filter <- .parseFilterInput(filter)
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
    function(x, filter=NULL) {
        filter <- .parseFilterInput(filter)
        filter <- .filter_list(filter)
        tx <- .xscripts(x, "ranges_tx", filter=filter)
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
    exon <- .xscripts(x, "ranges_exon", filter=filter)
    cds <- .xscripts(x, "ranges_cds", filter=filter)

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
        arrange_(~ tx_id, ~ exon_rank) %>% distinct()
}

.fiveUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    .UTRsByTranscript(x, filter, "-", "+")
}

#' @rdname extractors
#' @export
fiveUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .parseFilterInput(filter)
    .return_tbl(.fiveUTRsByTranscript_tbl(x, filter), filter)
}

#' @examples
#' ## fiveUTRsByTranscript
#' fiveUTRsByTranscript(src, filter = list(SymbolFilter("ADA")))
#'
#' @rdname extractors
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @exportMethod fiveUTRsByTranscript
setMethod("fiveUTRsByTranscript", "src_organism", function(x, filter=NULL) {
    filter <- .parseFilterInput(filter)
    .toGRangesList(x, "fiveUTRs", "tx", filter)
})

.threeUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .filter_list(filter)
    .UTRsByTranscript(x, filter, "+", "-")
}

#' @rdname extractors
#' @export
threeUTRsByTranscript_tbl <- function(x, filter = NULL) {
    filter <- .parseFilterInput(filter)
    .return_tbl(.threeUTRsByTranscript_tbl(x, filter), filter)
}

#' @rdname extractors
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @exportMethod threeUTRsByTranscript
setMethod("threeUTRsByTranscript", "src_organism", function(x, filter=NULL) {
    filter <- .parseFilterInput(filter)
    .toGRangesList(x, "threeUTRs", "tx", filter)
})
