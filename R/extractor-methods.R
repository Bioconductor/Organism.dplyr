#' @name Genomic-Extractors
#'
#' @title Extract genomic features from src_organism objects
#'
#' @aliases cds cds_tbl cdsBy cdsBy_tbl exons exons_tbl exonsBy
#'      exonsBy_tbl fiveUTRsByTranscript fiveUTRsByTranscript_tbl genes
#'      genes_tbl intronsByTranscript intronsByTranscript_tbl
#'      promoters promoters_tbl threeUTRsByTranscript threeUTRsByTranscript_tbl
#'      transcripts transcripts_tbl transcriptsBy transcriptsBy_tbl 
#'
#' @description
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
#' @usage
#'
#'  cds(x, ...)
#'  exons(x, ...)
#'  genes(x, ...)
#'  transcripts(x, ...)
#'  cds_tbl(x, filter=NULL, columns=NULL)
#'  exons_tbl(x, filter=NULL, columns=NULL)
#'  genes_tbl(x, filter=NULL, columns=NULL)
#'  transcripts_tbl(x, filter=NULL, columns=NULL)
#'  cdsBy(x, by=c("tx", "gene"), ...)
#'  exonsBy(x, by=c("tx", "gene"), ...)
#'  transcriptsBy(x, by=c("gene", "exon", "cds"), ...)
#'  cdsBy_tbl(x, by=c("tx", "gene"), filter=NULL, columns=NULL)
#'  exonsBy_tbl(x, by=c("tx", "gene"), filter=NULL, columns=NULL)
#'  transcriptsBy_tbl(x, by=c("gene", "exon", "cds"), filter=NULL, columns=NULL)
#'  promoters_tbl(x, upstream, downstream, filter=NULL, columns=NULL)
#'  intronsByTranscript_tbl(x, filter=NULL, columns=NULL)
#'  fiveUTRsByTranscript(x, ...)
#'  fiveUTRsByTranscript_tbl(x, filter=NULL, columns=NULL)
#'  threeUTRsByTranscript(x, ...)
#'  threeUTRsByTranscript_tbl(x, filter=NULL, columns=NULL)
#'
#' @param x A \code{src_organism} object
#'
#' @param filter Either NULL, \code{AnnotationFilter}, or
#'      \code{AnnotationFilterList} to be used to restrict the output. Filters
#'      consists of \code{AnnotationFilter}s and can be a \code{\link{GRanges}}
#'      object using "GRangesFilter" (see examples).
#'
#' @param columns A character vector indicating columns to be included in output
#'  \code{GRanges} object or \code{tbl}.
#'
#' @param upstream For \code{promoters()}: An integer(1) value indicating
#'  the number of bases upstream from the transcription start site.
#'
#' @param downstream For \code{promoters()}: An integer(1) value indicating
#'  the number of bases downstream from the transcription start site.
#'
#' @param by One of "gene", "exon", "cds" or "tx". Determines the grouping.
#'
#' @param ... Additional arguments to S4methods.  In this case, the same as
#'  \code{filter}.
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
#' @rdname Genomic-Extractors
#'
#' @examples
#' \dontrun{src <- src_ucsc("human")}
#' src <- src_organism(dbpath=hg38light())
#'
#' ## transcript coordinates with filter in tibble format
#' filters <- AnnotationFilter(~symbol == c("A1BG", "CDH2"))
#' transcripts_tbl(src, filters)
#'
#' transcripts_tbl(src, AnnotationFilter(~symbol %startsWith% "SNORD"))
#' transcripts_tbl(src, AnnotationFilter(~go == "GO:0005615"))
#' transcripts_tbl(src, filter=AnnotationFilter(
#'      ~symbol %startsWith% "SNORD" & tx_start < 25070000))
#'
#' ## transcript coordinates with filter in granges format
#' filters <- GRangesFilter(GenomicRanges::GRanges("chr15:1-25070000"))
#' transcripts(src, filters)
#'
#' ## promoters
#' promoters(src, upstream=100, downstream=50,
#'           filter = SymbolFilter("ADA"))
#'
#' ## transcriptsBy
#' transcriptsBy(src, by = "exon", filter = SymbolFilter("ADA"))
#'
#' ## exonsBy
#' exonsBy(src, filter = SymbolFilter("ADA"))
#'
#' ## intronsByTranscript
#' intronsByTranscript(src, filter = SymbolFilter("ADA"))
#'
#' @export cds cds_tbl cdsBy cdsBy_tbl exons exons_tbl exonsBy
#' @export exonsBy_tbl fiveUTRsByTranscript fiveUTRsByTranscript_tbl genes
#' @export genes_tbl intronsByTranscript intronsByTranscript_tbl
#' @export promoters promoters_tbl threeUTRsByTranscript
#' @export threeUTRsByTranscript_tbl transcripts transcripts_tbl transcriptsBy
#' @export transcriptsBy_tbl 
#'
#' @importFrom GenomeInfoDb seqinfo<- seqlevels seqlevels<-
#' @importFrom GenomicFeatures transcripts exons cds genes intronsByTranscript
#'      fiveUTRsByTranscript transcriptsBy exonsBy cdsBy promoters
#'      threeUTRsByTranscript
NULL

.tblFunctions <- function(where=topenv()) {
    tbl_names <- c('transcripts', 'exons', 'cds', 'genes')
    ranges_names <- c('ranges_tx', 'ranges_exon', 'ranges_cds', 'ranges_gene')
    by_names <- list(
        transcripts = c("gene", "exon", "cds"),
        exons = c("tx", "gene"),
        cds = c("tx", "gene"),
        genes = NULL
    )

    for (i in seq_along(tbl_names)) {
        eval(substitute({
            f <- function(x, filter = NULL, columns=NULL) {
                filter <- .parseFilterInput(filter, columns)
                table <- .xscripts_tbl(x, RANGES, filter)
                class(table) <- c("tbl_organism", class(table))
                table
            }
            name <- paste0(TBL, "_tbl")
            assign(name, f, WHERE)
            setMethod(TBL, "src_organism", where=WHERE,
                function(x, filter=NULL, columns=NULL) {
                    filter <- .parseFilterInput(filter, columns)
                    .toGRanges(x, .xscripts_tbl(x, RANGES, filter), filter)
                }
            )
            if (TBL != "genes") {
                f <- function() {
                    filter <- .parseFilterInput(filter, columns)
                    by <- match.arg(by)
                    table <- .xscriptsBy_tbl(x, RANGES, by, filter)
                    class(table) <- c("tbl_organism", class(table))
                    table
                }
                formals(f) <- alist(x = , by = BY, filter = NULL, columns = NULL)
                name <- paste0(TBL, "By_tbl")
                assign(name, f, WHERE)
                f <- function() {
                    filter <- .parseFilterInput(filter, columns)
                    by <- match.arg(by)
                    .toGRangesList(x, type = TBL, by = by, filter = filter)
                }
                formals(f) <- alist(x = , by = BY, filter = NULL, columns = NULL)
                setMethod(BY_METHOD, "src_organism", where=WHERE, f)
            }
        }, list(TBL = tbl_names[i],
                RANGES = ranges_names[i],
                WHERE = where,
                BY = by_names[[i]],
                BY_METHOD = paste0(tbl_names[i], "By"))))
    }
}

.tblFunctionsBy_UTR <- function(where=topenv()) {
    tbl_names <- c("fiveUTRs", "threeUTRs")
    strand <- list(
        fiveUTRs = c("-", "+"),
        threeUTRs = c("+", "-")
    )

    for (i in seq_along(tbl_names)) {
        eval(substitute({
            f <- function(x, filter = NULL, columns = NULL) {
                filter <- .parseFilterInput(filter, columns)
                filter <- .filter_list(filter)
                .UTRsByTranscript(x, filter, STRAND1, STRAND2)
            }
            name <- paste0(TBL, "ByTranscript_tbl")
            assign(name, f, WHERE)
            setMethod(METHOD_NAME, "src_organism", where=WHERE,
                function(x, filter = NULL, columns = NULL) {
                    filter <- .parseFilterInput(filter, columns)
                    .toGRangesList(x, TBL, "tx", filter)
                }
            )
        }, list(TBL = tbl_names[i],
                STRAND1 = strand[[i]][1],
                STRAND2 = strand[[i]][2],
                WHERE = where,
                METHOD_NAME = paste0(tbl_names[i], "ByTranscript"))))
    }
}

promoters_tbl <- function(x, upstream, downstream, filter = NULL, columns = NULL) {
    filter <- .parseFilterInput(filter, columns)
    .promoters_tbl(x, upstream, downstream, filter)
}

#' @importFrom S4Vectors mcols
#' @rdname Genomic-Extractors
setMethod("promoters", "src_organism",
    function(x, upstream, downstream, filter = NULL, columns = NULL) {
        filter <- .parseFilterInput(filter, columns)
        gr <- .toGRanges(x, .promoters_tbl(x, upstream, downstream, filter), filter)
        names(gr) <- mcols(gr)$tx_name
        gr
})

#' @importFrom GenomicRanges GRanges mcols<-
intronsByTranscript_tbl <-
    function(x, filter = NULL, columns = NULL) {
        filter <- .parseFilterInput(filter, columns)
        ans <- unlist(intronsByTranscript(x, filter))
        mcols(ans)[, "tx_id"] <- names(ans)
        unname(unlist(ans)) %>% as.data.frame %>% tibble::as_tibble %>%
            dplyr::select_(.dots = c('tx_id',
                              intron_chrom = 'seqnames',
                              intron_start = 'start',
                              intron_end = 'end',
                              intron_strand = 'strand'))
}

#' @rdname Genomic-Extractors
#'
#' @importFrom GenomicRanges mcols split
#' @importFrom IRanges psetdiff
setMethod("intronsByTranscript", "src_organism",
    function(x, filter=NULL, columns = NULL) {
        filter <- .parseFilterInput(filter, columns)
        filter <- .filter_list(filter)
        tx <- .xscripts(x, "ranges_tx", filter=filter)
        exn <- .xscriptsBy_tbl(x, "ranges_exon", by="tx", filter=filter)

        table <- tx %>% collect(n=Inf)

        fields <- unique(
            c("tx_id", "tx_chrom", "tx_start", "tx_end", "tx_strand",
              .filter_names(filter)))
        tx_gr <- do.call(select_, c(list(table), as.list(fields))) %>%
            filter_at(vars(c(1,2,3,4)), all_vars(!is.na(.))) %>% as("GRanges")

        tx_gr <- .updateSeqinfo(x, tx_gr)

        table <- exn %>% collect(n=Inf)

        fields <- unique(
            c("exon_chrom", "exon_start", "exon_end", "exon_strand", "tx_id",
             "exon_id", .filter_names(filter)))
        exn_grl <- do.call(select_, c(list(table), as.list(fields))) %>%
            filter_at(vars(c(1,2,3,4)), all_vars(!is.na(.))) %>% as("GRanges")
        exn_grl <- split(exn_grl, exn_grl$tx_id)

        tx_gr<- tx_gr[match(names(exn_grl), mcols(tx_gr)[, "tx_id"])]
        ans <- psetdiff(tx_gr, exn_grl)
        ans
})

#' @importFrom AnnotationFilter condition distributeNegation value
.parseFilterInput <- function(filter, columns) {
    if (is.language(filter))
        filter <- AnnotationFilter(filter)
    filter <- .filter_list(filter)
    if (is(filter, "AnnotationFilterList"))
        filter <- distributeNegation(filter)
    if (!is.null(columns) && length(columns) > 0L) {
        null_filters <- lapply(columns, function(i) {
            template <- ~filter != NULL
            template[[2]][[2]] <- as.name(i)
            AnnotationFilter(template)
        })
        null_filters <- do.call("AnnotationFilterList", null_filters)
        filter <- AnnotationFilterList(filter, null_filters)
    }
    filter
}

.tblFunctions(topenv())
.tblFunctionsBy_UTR(topenv())
