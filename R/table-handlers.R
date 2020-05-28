#################################################################
#################################################################
## Helper functions to manage temporary sqlite database in an   #
##   src_organism object.                                       #
#################################################################

#' @importFrom utils capture.output
#' @importFrom dbplyr sql_render
#' @importFrom dplyr show_query
.checkDepth <- function(table) {
    max <- 65
    output <- dbplyr::sql_render(table)
    #output <- tail(output, -1)
    output <- strsplit(output, "\n")
    max < length(output)
}

.iterTable <- function(x, table, doit=FALSE) {
    if (doit | .checkDepth(table)) {
        table <- table %>% collect(n=Inf)
        nam <- .getNewTableName()
        dbWriteTable(x$con, nam, table)
        res <- tbl(x$con, nam)
        # class(res) <- c("tbl_organism", class(res))
        res
    }
    else
        table
}

.getNewTableName <- local({
    pid <- Sys.getpid()
    id <- 0L
    function() {
        id <<- id + 1
        paste0("table", id, "_", pid)
    }
})

.getNewOutputName <- local({
    pid <- Sys.getpid()
    id <- 0L
    function() {
        id <<- id + 1
        paste0("output", id, "_", pid)
    }
})

.tablesForRemoval <- new.env(parent = emptyenv())

.removeTables <- function() {
    for (tbl in names(.tablesForRemoval)) {
        dbRemoveTable(.tablesForRemoval[[tbl]], tbl)
        rm(envir = .tablesForRemoval, list = tbl)
    }
}

.cleanOutput <- function(x, table) {
    table <- table %>% collect(n=Inf)
    name <- .getNewOutputName()
    dbWriteTable(x$db, name, table)
    table <- tbl(x$db, name)
    .deleteTempTables(x)
    e <- new.env()
    e[["src"]] <- x
    e[["name"]] <- name
    attr(table, "finalizer") <- e
    # onexit <- ifelse(x$dbpath != "", TRUE, FALSE)
    reg.finalizer(attr(table, "finalizer"), function(object) {
        .tablesForRemoval[[object[["name"]] ]] <- object[["src"]]$db
    })
    table
}

.deleteTempTables <- function(x) {
    tables <- dbListTables(x$con)
    tables <- tables[grep('table', tables)]
    for (i in tables)
        dbRemoveTable(x$con, i)
}
