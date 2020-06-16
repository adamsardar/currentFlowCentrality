# An excuse to have data.table imported smoothly. This is so as to prevent problems with testing and the `:=` operator
#' @import data.table
dummyFunc <- function(...){
  x <- data.table()
  invisible(I(...))
}