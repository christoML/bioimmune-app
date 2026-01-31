library(plumber)
library(jsonlite)

# Load plumber API
r <- plumb("api.R")

# ------------------------
# CORS filter
# ------------------------
r$filter("cors", function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "http://localhost:5173")
  res$setHeader("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type")

  if (req$REQUEST_METHOD == "OPTIONS") {
    return(list())
  }

  forward()
})

# Run server
r$run(host = "0.0.0.0", port = 8000)
