########################################################
########################################################
## Function to place filters and associated functions  #
##                                                     #
########################################################

#' @importFrom AnnotationFilter AnnotationFilter AnnotationFilterList field
#'      value convertFilter logicOp
#' @importFrom AnnotationDbi columns
#' @importFrom DBI dbListTables dbRemoveTable dbWriteTable
#' @importFrom dplyr "%>%" as_tibble arrange_at
#' @importFrom rlang .data parse_expr
.tbl_join <- function(x, filter, main_table, table_names) {
    if (is.null(filter))
        return(main_table)

    if (!.check_filters(x, filter))
        stop("Some filters are not avaiable.")

    ## Catch AnnotatioFilterList consisting of column filters.
    ## These filters do not need to be evaluated.
    if (.is_columns_afl(filter))
        return(.getAllTableValues(x, main_table, table_names, FALSE))

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
                    expr <- parse_expr(dplyr_filter)
                    val <-
                        tbl(x$db, selected_table[1]) %>%
                        filter(expr)
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

.is_columns_afl <- function(filter) {
    while (is(value(filter)[[1]], "AnnotationFilterList"))
        filter <- value(filter)[[1]]
    truth_vec <- vapply(filter, function(i) {
        if (!is(i, "AnnotationFilter"))
            FALSE
        else{
            if(condition(i) == "!=" && length(value(i)) == 0L)
                TRUE
            else
                FALSE
        }
    }, logical(1))
    all(truth_vec)
}

.getAllTableValues <- function(x, table, table_names, selected) {
    table_names <- table_names[!selected]
    for (i in names(table_names)) {
        tmp <- tbl(x$db, i)
        table <- .iterTable(x, table)
        table <- .leftJoin(table, tmp)
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

.updateSeqinfo <- function(x, gr) {
    seqinfo <- .getSeqinfo(x)
    seqlevels(gr) <- seqlevels(seqinfo)
    seqinfo(gr) <- seqinfo
    gr
}

.checkCompatibleStartEnds <- function(type, filter) {
    fields <- .fields(filter)
    starts <- fields[grep("start", fields)]
    ends <- fields[grep("end", fields)]
    starts <- grepl(type, starts)
    ends <- grepl(type, ends)
    all(starts) && all(ends)
}

########################################################
########################################################
## xscript* functions:                                 #
##  -- Perform different actions based on "main_ranges"#
##       for different extractor types                 #
##       (e.g.) main_ranges = "ranges_cds" for cds)    #
########################################################

.xscripts <- function(x, main_ranges, filter = NULL) {
    .removeTables()
    table_names <- .tableNames(x, filter, main_ranges)
    add_tables <- setdiff(names(table_names), dbListTables(x$db))
    for (i in add_tables) {
        table <- tbl(x, i, .load_tbl_only=TRUE) %>% collect(n=Inf)
        dbWriteTable(x$db, i, table)
    }
    table <- tbl(x$db, main_ranges)
    table <- .tbl_join(x, filter, table, table_names)
    table
}

.xscripts_tbl <- function(x, main_ranges, filter = NULL) {
    table <- .xscripts(x, main_ranges, filter)
    fields <- unique(unlist(c(.TBL_RETURNS[main_ranges],.filter_names(filter))))
    table <- .iterTable(x, table)
    table <- .cleanOutput(x, table)
    arrange_value <- .TBL_ARRANGE[main_ranges]
    if (is.null(arrange_value))
        arrange_value <- x$schema
    arrange_value <- unlist(arrange_value[lengths(arrange_value) != 0L])
    table %>%
        dplyr::select(fields) %>%
        distinct() %>%
        collect(n = Inf) %>%
        arrange_at(vars(arrange_value))
}

.xscriptsBy_tbl <- function(x, main_ranges, by, filter = NULL) {
    filter <- .filter_list(filter)
    table <- .xscripts(x, main_ranges, filter = filter)
    table <- .iterTable(x, table)
    if (main_ranges == "fiveUTRs") {
        filter <- .filter_list(filter)
        return(.UTRsByTranscript(x, filter, "-", "+"))
    }
    else if (main_ranges == "threeUTRs") {
        filter <- .filter_list(filter)
        return(.UTRsByTranscript(x, filter, "+", "-"))
    }
    else
        .byFun(x, table, main_ranges, by, filter = filter)
}

.byFun <- function(x, table, main_ranges, by, filter = NULL) {    
    .SPECIAL_FIELDS <- list(
        tx = c("tx_id", "exon_rank"),
        exon = c("exon_id", "exon_rank"),
        cds = c("cds_id", "exon_rank"),
        gene = c(x$schema)
    )
    other_main <- switch(by, tx = "ranges_tx",
                         exon = "ranges_exon",
                         cds = "ranges_cds",
                         gene = "ranges_gene")
    name <- strsplit(main_ranges, "_")[[1]][2]
    gr_names <- c("chrom", "start", "end", "strand", "id", "name")
    gr_names <- paste0(name, "_", gr_names)

    joinFun <- .leftJoin
    if (by == "gene")
        joinFun <- .innerJoin
    table <- joinFun(table, tbl(x, other_main))
    fields <- unique(
        unlist(c(gr_names, .SPECIAL_FIELDS[by], .filter_names(filter))))
    table %>%
        dplyr::select(fields) %>%
        arrange_at(vars(.SPECIAL_FIELDS[by][[1]])) %>%
        distinct()
}

########################################################
########################################################
## .toGRanges* functions:                              #
##  -- returns GRanges or GRangesList to extractors    #
########################################################

#' @importFrom dplyr across vars
#' @importFrom IRanges subsetByOverlaps
.toGRanges <- function(x, table, filter) {
    ## Filter out any rows that contain NA in chrom, start, end, or strand
    gr <-
        table %>%
        collect(n=Inf) %>%
        filter(across(1:4, ~ !is.na(.x))) %>%
        as("GRanges")
    .updateSeqinfo(x, gr)
}

.toGRangesList <- function(x, type, by, filter = NULL) {
    main_ranges <- switch (type, transcripts = "ranges_tx",
                   exons = "ranges_exon",
                   cds = "ranges_cds",
                   fiveUTRs = "fiveUTRs",
                   threeUTRs = "threeUTRs")

    if (type == "fiveUTRs") {
        table <- .UTRsByTranscript(x, filter, "-", "+")
        table <- .toGRanges(x, table, filter)
    }
    else if (type == "threeUTRs") {
        table <- .UTRsByTranscript(x, filter, "+", "-")
        table <- .toGRanges(x, table, filter)
    } else {
        table <- .xscriptsBy_tbl(x, main_ranges, by, filter)
        table <- .iterTable(x, table, doit=TRUE)
        table <- .toGRanges(x, table, filter)
    }

    f <- switch(by, gene=x$schema, exon="exon_id", cds="cds_id", tx="tx_id")
    grp <- mcols(table)[[f]]
    mcols(table)[[f]] <- NULL
    split(table, grp)
}

########################################################
########################################################
## Helpful lists for specific table values:            #
########################################################

.TBL_RETURNS <- list(
    ranges_tx = c("tx_chrom", "tx_start", "tx_end", "tx_strand", "tx_id", "tx_name"),
    ranges_exon = c("exon_chrom", "exon_start", "exon_end", "exon_strand", "exon_id"),
    ranges_cds = c("cds_chrom", "cds_start", "cds_end", "cds_strand", "cds_id"),
    ranges_gene = c("gene_chrom", "gene_start", "gene_end", "gene_strand"),
    promoters = c("chrom", "start", "end", "strand", "tx_id", "tx_name") 
)

.TBL_ARRANGE <- list(
    ranges_tx = "tx_id",
    ranges_exon = "exon_id",
    ranges_cds = "cds_id",
    ranges_gene = character()
)

########################################################
########################################################
## Promoter helper functions                           #
########################################################

#' @importFrom S4Vectors isSingleNumber
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
    table <-
        .xscripts_tbl(x, "ranges_tx", filter = filter) %>%
        mutate(
            start = ifelse(
                .data$tx_strand == "-",
                .data$tx_end - downstream + 1,
                .data$tx_start - upstream),
            end = ifelse(
                .data$tx_strand == "-",
                .data$tx_end + upstream,
                .data$tx_start + downstream - 1)
        ) %>%
        collect(n=Inf)

    table <-
        table %>%
        rename(chrom = .data$tx_chrom, strand = .data$tx_strand)

    fields <- unique(
        c("chrom", "start", "end", "strand", "tx_id", "tx_name",
          .filter_names(filter)))
    table <- .cleanOutput(x, table)
    table %>%
        dplyr::select(fields) %>%
        distinct() %>%
        arrange(.data$tx_id)
}

########################################################
########################################################
## UTR extractor helper functions:                     #
########################################################

.getSplicings <- function(x, filter=NULL) {
    exon <- .xscripts(x, "ranges_exon", filter=filter)
    cds <- .xscripts(x, "ranges_cds", filter=filter)

    exon_txid <- exon %>% dplyr::select("tx_id") %>% collect(n = Inf)
    exon_txid <- exon_txid[["tx_id"]]
    cds_txid <- cds %>% dplyr::select("tx_id") %>% collect(n = Inf)
    cds_txid <- cds_txid[["tx_id"]]
    exclude <- setdiff(exon_txid, cds_txid)

    table <- .leftJoin(exon, cds)
                       # by = c("tx_id", "exon_rank", .filter_names(filter)))

    fields <- unique(
        c("tx_id", "exon_rank", "exon_id", "exon_name", "exon_chrom",
          "exon_strand", "exon_start", "exon_end", "cds_id", "cds_start",
          "cds_end", .filter_names(filter)))
    splicings <- table %>% dplyr::select(fields) %>% collect(n=Inf)

    if(length(exclude) != 0)
        splicings <- splicings %>% filter(!.data$tx_id %in% exclude)
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
        mutate(start =
                   ifelse(!is.na(.data$cds_id) & .data$exon_strand == strand1,
                          .data$cds_end + 1L,
                          .data$exon_start),
               end =
                   ifelse(!is.na(.data$cds_id) & .data$exon_strand == strand2,
                          .data$cds_start - 1L,
                          .data$exon_end)) %>%
        filter(.data$start <= .data$end) %>% collect(n=Inf) %>% tibble::as_tibble()

    table <- rename(table, chrom = .data$exon_chrom)
    table <- rename(table, strand = .data$exon_strand)

    fields <- unique(
        c("chrom", "start", "end", "strand", "tx_id", "exon_id", "exon_name",
          "exon_rank", .filter_names(filter)))
    table %>% dplyr::select(fields) %>%
        arrange(.data$tx_id, .data$exon_rank) %>% distinct()
}
