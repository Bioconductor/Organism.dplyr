## create filter functions
.onLoad <- function(libname, pkgname) {
    for (field in c(.CHAR_FIELDS, .INT_FIELDS))  {
        class <- paste0(sub("([a-z])", "\\U\\1", field, perl=TRUE), "Filter")
        setClass(class, contains="BasicFilter")
        if (field %in% .CHAR_FIELDS)
            assign(class, .filterFactory(field, TRUE),
               envir=topenv(parent.frame()))
        else
            assign(class, .filterFactory(field, FALSE),
               envir=topenv(parent.frame()))

    }
}
