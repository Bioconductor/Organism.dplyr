#################################################################
#################################################################
## Helper functions to manage temporary sqlite database in an   #
##   src_organism object.                                       #
#################################################################

#' @importFrom utils capture.output
#' @importFrom dplyr show_query
.checkDepth <- function(table) {
    max <- 65
    output <- capture.output(table %>% show_query(), type="message")
    output <- tail(output, -1)
    max < length(output)
}

.iterTable <- function(x, table, doit=FALSE) {
    if (doit | .checkDepth(table)) {
        table <- table %>% collect(n=Inf)
        nam <- .getNewTableName()
        dbWriteTable(x$db, nam, table)
        res <- tbl(x$db, nam)
        # class(res) <- c("tbl_organism", class(res))
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

.getNewOutputName <- local({
    id <- 0L
    function() {
        id <<- id + 1
        paste0("output", id)
    }
})

.cleanOutput <- function(x, table) {
    table <- table %>% collect(n=Inf)
    name <- .getNewOutputName()
    dbWriteTable(x$db, name, table)
    .deleteTempTables(x)
    e <- new.env()
    attr(e, "src") <- x
    attr(e, "name") <- name
    attr(table, "finalizer") <- e
    reg.finalizer(attr(table, "finalizer"), function(object) {
        dbRemoveTable(attr(object, "src")$db, attr(object, "name"))
    })
    table
}

.deleteTempTables <- function(x) {
    tables <- dbListTables(x$db)
    tables <- tables[grep('table', tables)]
    for (i in tables)
        dbRemoveTable(x$db, i)
}
