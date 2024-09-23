constructCollectionUtils <- function() {
  
  public <- list()
  
  public$lapply <- function(x, func, reporter = "generic lapply", verbose = TRUE) {
    count <- 0
    x %>% lapply(function(currEntry) {
      count <<- count + 1
      if (verbose) {
        print(paste0(reporter, ": processing item ", count, " of ", length(x)))
      }
      func(currEntry)
    })
  }
  
  public$lapplyWithName <- function(x, func, reporter = "generic lapplyWithName", verbose = TRUE) {
    # func should take in 2 arguments - function(currEntryName, currItem)
    entries <- x %>% names()
    names(entries) <- entries
    entries %>% public$lapply(function(currEntryName) {
      func(currEntryName, x[[currEntryName]])
    }, reporter, verbose)
  }
  
  public$sapply <- function(x, func, reporter = "generic sapply") {
    count <- 0
    x %>% sapply(function(currEntry) {
      count <<- count + 1
      print(paste0(reporter, ": processing item ", count, " of ", length(x)))
      func(currEntry)
    })
  }
  
  public$foreach <- function(x, func) {
    for (i in x) {
      func(i)
    }
  }
  
  return(public)
}
