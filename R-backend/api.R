library(plumber)
library(jsonlite)
library(dplyr)
library(GEOquery)
library(R.utils)
library(tcltk)
library(oligo)
library(affy)
library(affyio)
library(Biobase)

# Load custom functions
source("GEO_lib_functions.R")
source("graph_functions.R")

#* Download GEO datasets with logs
#* @param gse_ids Comma-separated string of GSE IDs, e.g., GSE123,GSE456
#* @post /datasets
function(req, res) {
  body <- req$postBody
  if (is.null(body) || nchar(body) == 0) {
    res$status <- 400
    return(list(error = "No GSE IDs provided"))
  }
  
  parsed <- tryCatch(fromJSON(body), error = function(e) {
    res$status <- 400
    return(list(error = "Invalid JSON"))
  })
  
  if (is.null(parsed$gse_ids)) {
    res$status <- 400
    return(list(error = "Missing 'gse_ids' field"))
  }
  
  gse_ids <- trimws(unlist(strsplit(parsed$gse_ids, ",")))
  
  result <- tryCatch({
    download_GEO(gse_ids, type = c("series_matrix", "CEL"))
  }, error = function(e) {
    res$status <- 500
    return(list(error = e$message))
  })
  
  return(list(
    success = TRUE,
    datasets = names(result$datasets),
    logs = result$logs
  ))
}
