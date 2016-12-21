## transcripts_tbl(src, filter=SymbolFilter("PTEN"))
## transcripts_tbl(src, filter=SymbolFilter("BRCA", "startsWith"))
## transcripts_tbl(src, filter=list(
##     SymbolFilter("BRCA", "startsWith"),
##     EntrezFilter(123)))

#' @export
setClass("BasicFilter",
         representation(
             "VIRTUAL",
             field="character",
             condition="character",
             value="character",
             .valueIsCharacter="logical"
         ),
         prototype=list(
             condition= "==",
             value=character(),
             .valueIsCharacter=TRUE
         )
)

.condition <- function(x) x@condition
.value <- function(x) x@value

setValidity("BasicFilter", function(object) {
    txt <- character()
    if (length(.condition(object)) != 1L)
        txt <- c(txt, "'condition' must be length 1")
    if (!.condition(object) %in% .OPS)
        txt <- c(txt,
                 sprintf("'condition' must be one of %s",
                         paste("'", .OPS, "'", collapse=", ")))
    if (length(txt)) txt else TRUE
})

## FIXME: show,BasicFilter-method

.filterFactory <- function(field, isCharacter) {
    class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
    .valueIsCharacter <- isCharacter
    function(value, condition="==") {
        new(class,
            field=field,
            condition=condition,
            value=as.character(value),
            .valueIsCharacter=.valueIsCharacter)
    }
}

.OPS <- c("==", "!=", "startsWith")
.CHAR_FIELDS <- c("symbol", "entrez", "ensembl")
.INT_FIELDS <- character()

#' @export SymbolFilter EntrezFilter EnsemblFilter
.onLoad <- function(libname, pkgname) {
    for (field in .CHAR_FIELDS) {
        class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
        setClass(class, contains="BasicFilter")
        assign(class, .filterFactory(field, TRUE),
               envir=topenv(parent.frame()))
    }
}

.convertFilter <- function(filter) {
    field <- class(filter)[1]
    field <- substr(field, 1, nchar(field) - 6)
    value <- .value(filter)
    condition <- .condition(filter)
    
    op <- switch(
        condition,
        "==" = if (length(value) == 1) "==" else "%in%", 
        "!=" = if (length(value) == 1) "!=" else "%in%",
        "startsWith" = "%like%"
    )
    
    if (condition != "startsWith") 
        value <- paste0("'", value, "'", collapse=", ")
    
    if ((condition == "!=") && op == "%in%")
        sprintf("!%s %s c(%s)", field, op, value)
    else if (condition == "startsWith")
        sprintf("%s %s %s", field, op, paste0("'", value, "%", "'"))
    else 
        sprintf("%s %s c(%s)", field, op, value)
}


.tbl_filter <- function(keep1, value) {
    values <- paste0("'", value, "'", collapse=", ")
    op <- if (length(value) == 1) "==" else "%in%"
    sprintf("%s %s c(%s)", keep1, op, values)
}

.filter_names <- function(filter) {
    setdiff(names(filter), "granges")
}

.tbl_join <- function(x, table, tbls, filter) {
    if (is.null(filter))
        return(table)
    
    if (is(filter, "BasicFilter")) {
        field <- class(filter)[1]
        fields <- substr(field, 1, nchar(field) - 6)
    } else
        fields <- names(filter)
    
    ## filter by fields from main table
    fields1 <- fields[fields %in% colnames(table)]
    if (length(fields1) != 0) {
        filters <- sapply(fields1, .tbl_filter, filter[[fields1]])
        filters <- paste0(filters, collapse=" & ")
        table <- table %>% filter_(filters)
        filter <- filter[setdiff(fields, fields1)]
    }
    
    ## filter by fields from other tables
    # fields <- .filter_names(filter)
    for (i in tbls) {
        keep <- fields[fields %in% colnames(tbl(x, i))]
        if (is.null(keep) || length(keep) == 0)
            next
        
        if (is(filter, "BasicFilter")) {
            filters <- .convertFilter(filter)
        } else {
            filters <- sapply(keep, .tbl_filter, filter[[keep]])
            filters <- paste0(filters, collapse=" & ")
        }
        
        table <- inner_join(table, tbl(x, i)) %>% filter_(filters)
        fields <- setdiff(fields, keep)
    }
    
    table
}





