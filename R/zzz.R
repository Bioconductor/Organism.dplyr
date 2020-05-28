#' @importFrom AnnotationFilter supportedFilters
.onLoad <-
    function(...)
{
    .filter_init()
#    .tbl_Functions()
}

.onUnload <-
    function(...)
{
    .removeTables()
}
