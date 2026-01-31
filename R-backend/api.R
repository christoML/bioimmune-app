library(plumber)
library(jsonlite)
library(dplyr)
library(GEOquery)
library(R.utils)
library(oligo)
library(affy)
library(affyio)
library(Biobase)

# Load custom functions
source("GEO_lib_functions.R")
source("graph_functions.R")

# -----------------------------
# Endpoint: Download GEO datasets
# -----------------------------
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

# -----------------------------
# Endpoint: Get Metadata (group_vector)
# -----------------------------
#* Get group vector metadata from downloaded datasets
#* @post /get_geo_groups
function(req, res) {
  tryCatch({
    # Call the function from GEO_lib_functions.R
    # This function should handle reading datasets and saving R_objects/group_vector.rds
    result <- get_geo_groups()  
    
    # Return only the summary to frontend
    list(success = TRUE, logs = result$logs)
    
  }, error = function(e) {
    res$status <- 400
    list(success = FALSE, error = e$message)
  })
}

# -----------------------------
# Endpoint: Gene Matrix Preprocessing (RMA)
# -----------------------------
#* Run gene_mat_preprocess on .CEL or series_matrix files
#* @param type:string Either "CEL" or "series_matrix"
#* @post /gene_mat_preprocess
function(req, res, type = "CEL") {
  tryCatch({
    # Ensure type is either CEL or series_matrix
    type <- match.arg(type, c("CEL", "series_matrix"))

    # Call the existing function
    result <- gene_mat_preprocess(input_type = type)

    # Return logs to frontend
    list(success = TRUE, logs = result$logs)

  }, error = function(e) {
    res$status <- 400
    list(success = FALSE, error = e$message)
  })
}
