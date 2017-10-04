#' @importFrom dplyr inner_join union
## RSQLite safe joins
.fullJoin <- function(table1, table2) {
    ## just a union() works since all tables should have same columns
    union(table1, table2)
}

.innerJoin <- function(table1, table2) {
    ## find columns that contain only NAs to prevent a bad join on them
    trouble_makers <- c('tx_type', 'exon_name', 'cds_name', 'refseq',
                        'ensembltrans', 'enzyme', 'ensemblprot', 'prosite')
    col1 <- colnames(table1)
    col2 <- colnames(table2)
    trouble_makers <-
        trouble_makers[trouble_makers %in% col1 & trouble_makers %in% col2]
    if(length(trouble_makers) > 0) {
        nam <- col1[!col1 %in% trouble_makers]
        table1 <- do.call(dplyr::select, c(list(table1), as.list(nam)))
    }
    inner_join(table1, table2)
}
