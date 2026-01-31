#=========================================================
# FUNCTION: download_GEO (multi-type version)
#=========================================================
download_GEO <- function(gse_ids, type = c("series_matrix", "CEL")) {
  type <- match.arg(type, several.ok = TRUE)
  
  # Always download inside 'datasets' folder relative to this script
  download_dir <- file.path(dirname(normalizePath("GEO_lib_functions.R")), "datasets")
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  options(timeout = 600)  # increase timeout for large files
  
  logs <- c()          # collect logs
  gse_objects <- list()  # store datasets
  
  # Logger collects messages
  logger <- function(msg) {
    logs <<- c(logs, msg)
    cat(msg, "\n")  # optional: still print to R console
  }
  
  logger("📦 Starting GEO download...")
  
  for (gse_id in gse_ids) {
    logger(paste0("\n🔽 Processing ", gse_id, "..."))
    gse_objects[[gse_id]] <- list()
    
    for (t in type) {
      tryCatch({
        if (t == "series_matrix") {
          gse_obj <- getGEO(gse_id, destdir = download_dir)
          gse_objects[[gse_id]][[t]] <- gse_obj
          
          # Unzip .gz files
          gz_files <- list.files(download_dir, pattern = "\\.gz$", full.names = TRUE)
          if (length(gz_files) > 0) for (f in gz_files) gunzip(f, remove = TRUE, overwrite = TRUE)
          
          logger(paste0("✅ ", gse_id, " (series_matrix) processed successfully."))
          
        } else if (t == "CEL") {
          gse_dir <- file.path(download_dir, paste0(gse_id, "-CEL"))
          dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)
          getGEOSuppFiles(gse_id, baseDir = gse_dir, makeDirectory = FALSE)
          
          # Extract TAR files
          tar_files <- list.files(gse_dir, pattern = "\\.tar$", full.names = TRUE)
          if (length(tar_files) > 0) for (tarf in tar_files) { untar(tarf, exdir = gse_dir); file.remove(tarf) }
          
          # Unzip .gz files
          gz_files <- list.files(gse_dir, pattern = "\\.gz$", full.names = TRUE)
          if (length(gz_files) > 0) for (f in gz_files) gunzip(f, remove = TRUE, overwrite = TRUE)
          
          gse_objects[[gse_id]][[t]] <- gse_dir
          logger(paste0("✅ ", gse_id, " (CEL) processed successfully."))
        }
      }, error = function(e) {
        logger(paste0("❌ Failed for ", gse_id, " (", t, "): ", e$message))
        gse_objects[[gse_id]][[t]] <- NULL
      })
    }
  }
  
  logger("\n📁 All downloads complete.")
  logger(paste0("✅ Successfully processed datasets: ", paste(names(gse_objects), collapse = ", ")))
  
  return(list(datasets = gse_objects, logs = logs))
}


#=========================================================
# FUNCTION: get_geo_groups()
#=========================================================
get_geo_groups <- function(pattern = NULL) {
  
  # ---- GUI: select folder once ----
  data_dir <- tclvalue(
    tkchooseDirectory(title = "Select folder containing series_matrix.txt files")
  )
  if (data_dir == "") stop("No folder selected. Aborting.")
  
  sm_files <- list.files(
    data_dir,
    pattern = "series_matrix\\.txt$",
    full.names = TRUE
  )
  
  if (length(sm_files) == 0) {
    stop("No series_matrix.txt files found in selected folder.")
  }
  
  all_groups <- list()
  
  for (sm_file in sm_files) {
    gse_id <- sub("_series_matrix.txt$", "", basename(sm_file))
    cat("📄 Reading:", gse_id, "\n")
    
    gse <- getGEO(filename = sm_file, GSEMatrix = TRUE)
    eset <- if (inherits(gse, "ExpressionSet")) gse else gse[[1]]
    pheno <- pData(eset)
    
    # ---- Identify metadata column ----
    if ("characteristics_ch1.2" %in% colnames(pheno)) {
      col_use <- "characteristics_ch1.2"
    } else if ("characteristics_ch1" %in% colnames(pheno)) {
      col_use <- "characteristics_ch1"
    } else if (!is.null(pattern)) {
      col_use <- grep(pattern, colnames(pheno), value = TRUE)[1]
    } else {
      stop("Cannot find a suitable metadata column. Please specify 'pattern'.")
    }
    
    meta_col <- pheno[[col_use]]
    
    group <- ifelse(
      grepl("control", meta_col, ignore.case = TRUE), "NC",
      ifelse(grepl("rheumatoid", meta_col, ignore.case = TRUE), "RA", NA)
    )
    
    keep <- !is.na(group)
    group <- factor(group[keep], levels = c("NC", "RA"))
    names(group) <- rownames(pheno)[keep]
    
    cat("Group summary for", gse_id, ":\n")
    print(table(group))
    
    all_groups[[length(all_groups) + 1]] <- group
  }
  
  group_vector <- do.call(c, all_groups)
  
  return(group_vector)
}

#================FUNCTION===============#
#====preprocess_and_combat_GUI==========#
preprocess_and_combat_GUI <- function(group_vector = NULL, input_paths = NULL) {
  #-----------------------------
  # Step 0: Select output directory
  #-----------------------------
  out_dir <- tclvalue(tkchooseDirectory(title = "Select folder to save gene-level XLSX files"))
  if (is.null(out_dir) || length(out_dir) == 0 || out_dir == "") stop("No folder selected. Aborting.")
  cat("📂 Output directory:", out_dir, "\n")
  
  exprs_list <- list()
  annotation_list <- list()
  
  #-----------------------------
  # Step 1: Select CEL folders if not provided
  #-----------------------------
  cel_folders <- input_paths
  if (is.null(cel_folders)) {
    parent_dir <- tclvalue(tkchooseDirectory(title = "Select folder containing GSE -CEL folders"))
    if (is.null(parent_dir) || length(parent_dir) == 0 || parent_dir == "") stop("No folder selected. Aborting.")
    
    cel_folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
    cel_folders <- cel_folders[grepl("-CEL$", cel_folders)]
  }
  
  if (length(cel_folders) == 0) stop("No '-CEL' folders found.")
  cat("Selected CEL folders:\n")
  print(cel_folders)
  
  #-----------------------------
  # Step 2: Read and RMA normalize each CEL dataset
  #-----------------------------
  for (cel_path in cel_folders) {
    gse_id <- basename(cel_path)
    cat("\n📦 Processing CEL dataset:", gse_id, "\n")
    
    cel_files <- list.celfiles(cel_path, full.names = TRUE)
    if (length(cel_files) == 0) {
      warning("No CEL files found for ", gse_id)
      next
    }
    
    # Check platform
    cat("Checking CEL file chip types...\n")
    chip_types <- sapply(cel_files, function(f) read.celfile.header(f)$cdfName)
    platform <- unique(chip_types)
    print(table(chip_types))
    
    if (!any(chip_types %in% c("HG-U133A", "HG-U133_Plus_2"))) {
      cat("Skipping", gse_id, "because it is neither HG-U133A (GPL96) nor HG-U133_Plus_2 (GPL570).\n")
      next
    }
    
    if (any(platform == "HG-U133A")) {
      cel_use <- names(chip_types[chip_types == "HG-U133A"])
      annotation_db <- hgu133a.db
      cat("Found", length(cel_use), "HG-U133A CEL files.\n")
    } else if (any(platform == "HG-U133_Plus_2")) {
      cel_use <- names(chip_types[chip_types == "HG-U133_Plus_2"])
      annotation_db <- hgu133plus2.db
      cat("Found", length(cel_use), "HG-U133_Plus_2 CEL files.\n")
    }
    
    # Read CELs and RMA normalize
    raw <- ReadAffy(filenames = cel_use)
    norm_data <- affy::rma(raw)
    exprs_mat <- exprs(norm_data)
    
    exprs_list[[gse_id]] <- exprs_mat
    annotation_list[[gse_id]] <- annotation_db
  }
  
  #-----------------------------
  # Step 3: Merge probe-level matrices
  #-----------------------------
  common_probes <- Reduce(intersect, lapply(exprs_list, rownames))
  cat("Number of common probes across datasets:", length(common_probes), "\n")
  
  exprs_common <- do.call(cbind, lapply(exprs_list, function(x) x[common_probes, , drop=FALSE]))
  
  # Remove any prefixes or suffixes, keep just GSM ID
  colnames(exprs_common) <- sub(".*(GSM\\d+).*", "\\1", colnames(exprs_common))
  # Filter only samples present in group_vector
  keep_samples <- names(group_vector)[names(group_vector) %in% colnames(exprs_common)]
  exprs_common <- exprs_common[, keep_samples]
  group_vector <- group_vector[keep_samples]
  
  #-----------------------------
  # Step 4: Create batch vector
  #-----------------------------
  # Create a batch vector based on which GSE each GSM came from
  batch <- sapply(colnames(exprs_common), function(gsm_id) {
    # Search through exprs_list for the GSE that contains this GSM
    gse <- names(exprs_list)[sapply(exprs_list, function(mat) {
      # Clean the column names of mat the same way
      clean_mat_cols <- sub(".*(GSM\\d+).*", "\\1", colnames(mat))
      gsm_id %in% clean_mat_cols
    })]
    # Return the GSE name (should be length 1)
    gse
  })
  batch <- factor(batch)
  cat("Batch composition:\n")
  print(table(batch))
  
  
  #-----------------------------
  # Step 5: ComBat batch correction at probe level
  #-----------------------------
  mod <- model.matrix(~ group_vector)
  combat_probe <- ComBat(dat = as.matrix(exprs_common),
                         batch = batch,
                         mod = mod,
                         par.prior = TRUE,
                         prior.plots = FALSE)
  cat("✅ Probe-level ComBat complete.\n")
  
  #-----------------------------
  # Step 6: Map probes → genes AFTER ComBat
  #-----------------------------
  annotation_combined <- c()
  for (gse in names(exprs_list)) {
    probes <- rownames(exprs_list[[gse]])
    annot_db <- annotation_list[[gse]]
    annotation_combined <- c(annotation_combined, mapIds(annot_db,
                                                         keys = probes,
                                                         column = "SYMBOL",
                                                         keytype = "PROBEID",
                                                         multiVals = "first"))
  }
  
  annotation_combined <- annotation_combined[rownames(combat_probe)]
  valid <- !is.na(annotation_combined)
  combat_probe <- combat_probe[valid, , drop=FALSE]
  gene_symbols <- annotation_combined[valid]
  
  combat_gene <- rowsum(combat_probe, group = gene_symbols)
  
  #-----------------------------
  # Step 7: Save final gene-level matrix
  #-----------------------------
  out_file <- file.path(out_dir, "ComBat_final_gene_matrix.xlsx")
  write.xlsx(combat_gene, out_file, rowNames = TRUE)
  cat("📁 Final gene-level matrix saved to:", out_file, "\n")
  
  #-----------------------------
  # Step 8: Optional QC plotting
  #-----------------------------
  par(mfrow = c(1,2))
  boxplot(exprs_common, main="Before ComBat", outline=FALSE, col="lightcoral")
  boxplot(combat_probe, main="After ComBat", outline=FALSE, col="lightseagreen")
  par(mfrow = c(1,1))
  
  return(combat_gene)
}


#=========================================================
# FUNCTION: gene_mat_preprocess()
#=========================================================
gene_mat_preprocess <- function(input_type = c("CEL", "series_matrix"), input_paths = NULL) {
  input_type <- match.arg(input_type)
  
  #------------------ STEP 0 — Select output directory ------------------
  out_dir <- tclvalue(tkchooseDirectory(title = "Select folder to save gene-level XLSX files"))
  if (out_dir == "") stop("No folder selected. Aborting.")
  
  cat("📂 Output directory:", out_dir, "\n")
  
  exprs_list <- list()
  
  #------------------ STEP 1 — Process CEL folders ------------------
  if (input_type == "CEL") {
    cel_folders <- input_paths
    if (is.null(cel_folders)) {
      cel_folders <- list.dirs(
        tclvalue(tkchooseDirectory(title = "Select folder containing GSE -CEL folders")),
        full.names = TRUE,
        recursive = FALSE
      )
      cel_folders <- cel_folders[grepl("-CEL$", cel_folders)]
    }
    if (length(cel_folders) == 0) stop("No '-CEL' folders found.")
    
    exprs_list <- list()
    
    for (cel_path in cel_folders) {
      gse_id <- basename(cel_path)
      cat("\n📦 Processing CEL dataset:", gse_id, "\n")
      
      cel_files <- list.celfiles(cel_path, full.names = TRUE)
      if (length(cel_files) == 0) {
        warning("No CEL files found for ", gse_id)
        next
      }
      
      # Check platform
      cat("Checking CEL file chip types...\n")
      chip_types <- sapply(cel_files, function(f) read.celfile.header(f)$cdfName)
      platform <- unique(chip_types)
      print(table(chip_types))
      
      # Keep only supported platforms
      if (!any(chip_types %in% c("HG-U133A", "HG-U133_Plus_2"))) {
        cat("Skipping", gse_id,
            "because it is neither HG-U133A (GPL96) nor HG-U133_Plus_2 (GPL570).\n")
        next
      }
      
      if (any(platform == "HG-U133A")) {
        cel_use <- names(chip_types[chip_types == "HG-U133A"])
        annotation_db <- hgu133a.db
        cat("Found", length(cel_use), "HG-U133A CEL files.\n")
      } else if (any(platform == "HG-U133_Plus_2")) {
        cel_use <- names(chip_types[chip_types == "HG-U133_Plus_2"])
        annotation_db <- hgu133plus2.db
        cat("Found", length(cel_use), "HG-U133_Plus_2 CEL files.\n")
      }
      
      # Read CELs (raw, no normalization)
      raw_data <- ReadAffy(filenames = cel_use)
      exprs_data <- exprs(raw_data)
      norm_data <- affy::rma(raw_data)
      exprs_data <- exprs(norm_data)
      cat("RMA complete. Matrix dimensions:", dim(exprs_data), "\n")
      
     #------------------ Probe → Gene Mapping ------------------
    cat("Mapping probes to gene symbols...\n")
    probe_ids <- rownames(exprs_data)
    gene_symbols <- mapIds(
      annotation_db,
      keys = probe_ids,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
    
    valid <- !is.na(gene_symbols)
    exprs_data <- exprs_data[valid, , drop = FALSE]
    gene_symbols <- gene_symbols[valid]
    
    # Collapse duplicates (average by gene symbol)
    exprs_gene <- rowsum(exprs_data, group = gene_symbols)
    
    cat("Final gene-level matrix dimensions:", dim(exprs_gene), "\n")
    
    #------------------ Clean Column Names ------------------
    cat("Cleaning sample (column) names...\n")
    coln <- colnames(exprs_gene)
    coln <- sub("(GSM\\d+).*", "\\1", coln)
    coln <- sub("^X", "", coln)
    colnames(exprs_gene) <- coln
      
      # Save XLSX
      out_file <- file.path(out_dir, paste0(gse_id, "_gene_matrix.xlsx"))
      write.xlsx(exprs_gene, out_file, rowNames = TRUE)
      cat("✅ Saved:", out_file, "\n")
      
      exprs_list[[gse_id]] <- exprs_gene
    }
  }
  
  
  #------------------ STEP 2 — Process Series Matrix Files ------------------
  if (input_type == "series_matrix") {
    # If no files provided, ask user to select folder
    if (is.null(input_paths)) {
      series_dir <- tclvalue(tkchooseDirectory(title = "Select folder containing series_matrix files"))
      if (series_dir == "") stop("No folder selected. Aborting.")
      input_paths <- list.files(series_dir, pattern = "series_matrix.*\\.txt$", full.names = TRUE)
    }
    
    if (length(input_paths) == 0) stop("No series_matrix files found.")
    
    exprs_list <- list()
    
    for (f in input_paths) {
      gse_id <- tools::file_path_sans_ext(basename(f))
      cat("\n📦 Processing series_matrix dataset:", gse_id, "\n")
      
      # Load the local series_matrix file directly
      expr <- read.delim(f, comment.char = "!", stringsAsFactors = FALSE, row.names = 1)
      
      # Determine platform and annotation
      annot_db <- hgu133a.db
      
      # Map probes → genes
      probe_ids <- rownames(expr)
      gene_symbols <- mapIds(
        annot_db,
        keys = probe_ids,
        column = "SYMBOL",
        keytype = "PROBEID",
        multiVals = "first"
      )
      
      valid <- !is.na(gene_symbols)
      expr <- expr[valid, , drop = FALSE]
      gene_symbols <- gene_symbols[valid]
      
      # Collapse duplicates by averaging
      exprs_gene <- rowsum(expr, group = gene_symbols)
      
      # Clean GSM column names
      coln <- colnames(exprs_gene)
      coln <- sub("(GSM\\d+).*", "\\1", coln)
      coln <- sub("^X", "", coln)
      colnames(exprs_gene) <- coln
      
      # Save XLSX
      out_file <- file.path(out_dir, paste0(gse_id, "_gene_matrix.xlsx"))
      write.xlsx(exprs_gene, out_file, rowNames = TRUE)
      cat("✅ Saved:", out_file, "\n")
      
      exprs_list[[gse_id]] <- exprs_gene
    }

    
    cat("\n📁 All series_matrix datasets processed.\n")
  }
  
  cat("\n📁 All datasets processed. Gene-level XLSX files are ready.\n")
  
  
  #------------plots for RMA-------------------#
  # Set up side-by-side plots
  par(mfrow = c(1, 2), mar = c(10, 4, 4, 2)) # rotate x labels
  
  # Boxplot of raw CEL intensities
  boxplot(raw_data,
          outline = FALSE,
          col = "lightblue",
          las = 2,
          main = "Raw CEL Intensities",
          ylab = "Expression")
  
  # Boxplot of RMA-normalized data
  boxplot(exprs_data,
          outline = FALSE,
          col = "lightgreen",
          las = 2,
          main = "RMA-normalized Intensities",
          ylab = "Expression")
  
  # Reset plotting layout
  par(mfrow = c(1, 1))
  return(exprs_list)
}


#=========================================================
# FUNCTION: merge_and_combat()
#=========================================================
merge_and_combat2 <- function(group_vector,
                             use_combat = TRUE) {
  
  #------------------ STEP 1 — Select Directory ------------------
  data_dir <- tclvalue(
    tkchooseDirectory(title = "Select folder containing gene-level XLSX files")
  )
  if (data_dir == "") stop("No folder selected. Aborting.")
  
  cat("📂 Selected directory:", data_dir, "\n")
  
  #------------------ STEP 2 — Locate XLSX Files ------------------
  files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) < 2)
    stop("Need at least two gene-matrix XLSX files to merge.")
  
  cat("Found", length(files), "files:\n")
  print(basename(files))
  
  #------------------ STEP 3 — Read Matrices ------------------
  exprs_list <- lapply(files, function(f) {
    as.data.frame(read.xlsx(f, rowNames = TRUE))
  })
  
  names(exprs_list) <- basename(files)
  
  #------------------ STEP 4 — Find Common Genes ------------------
  common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
  cat("Number of common genes:", length(common_genes), "\n")
  
  exprs_list_common <- lapply(exprs_list, function(x) {
    x[common_genes, , drop = FALSE]
  })
  
  #------------------ STEP 5 — Merge ------------------
  combined_exprs <- do.call(cbind, exprs_list_common)
  
  # Clean GSM names
  colnames(combined_exprs) <- sub(".*(GSM\\d+).*", "\\1", colnames(combined_exprs))
  
  cat("Combined dataset dimensions:", dim(combined_exprs), "\n")
  
  #------------------ STEP 6 — Align group vector ------------------
  keep_samples <- intersect(colnames(combined_exprs), names(group_vector))
  
  if (length(keep_samples) < 2)
    stop("Too few samples after matching group_vector.")
  
  combined_exprs <- combined_exprs[, keep_samples, drop = FALSE]
  group_vector <- factor(group_vector[keep_samples])
  
  combined_exprs <- combined_exprs[, names(group_vector)]
  
  cat("Samples per group:\n")
  print(table(group_vector))
  
  #------------------ STEP 7 — Create batch vector ------------------
  batch <- sapply(colnames(combined_exprs), function(gsm) {
    names(exprs_list_common)[
      sapply(exprs_list_common, function(mat) gsm %in% colnames(mat))
    ][1]
  })
  batch <- factor(batch)
  
  cat("Batch composition:\n")
  print(table(batch))
  
  #==============================
  # OPTIONAL: ComBat correction
  #==============================
  if (use_combat) {
    
    cat("\n🧬 Applying ComBat batch correction...\n")
    
    mod <- model.matrix(~ group_vector)
    
    combat_exprs <- ComBat(
      dat = as.matrix(combined_exprs),
      batch = batch,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
    
    out_file <- file.path(data_dir, "Merged_ComBat_matrix.xlsx")
    write.xlsx(combat_exprs, out_file, rowNames = TRUE)
    
    cat("✅ ComBat-corrected matrix saved to:\n", out_file, "\n")
    
    # Boxplots
    par(mfrow = c(1, 2), mar = c(10, 4, 4, 2))
    boxplot(combined_exprs, outline = FALSE, las = 2,
            main = "Before ComBat", ylab = "Expression")
    boxplot(combat_exprs, outline = FALSE, las = 2,
            main = "After ComBat", ylab = "Expression")
    par(mfrow = c(1, 1))
    
    return(combat_exprs)
    
  } else {
    
    cat("\n📦 Skipping ComBat (merge-only mode)\n")
    
    out_file <- file.path(data_dir, "Merged_noComBat_matrix.xlsx")
    write.xlsx(combined_exprs, out_file, rowNames = TRUE)
    
    cat("✅ Merged matrix saved to:\n", out_file, "\n")
    
    boxplot(combined_exprs, outline = FALSE, las = 2,
            main = "Merged expression (no ComBat)",
            ylab = "Expression")
    
    return(combined_exprs)
  }
}















#=========================================================
# FUNCTION: merge_and_combat()
#=========================================================
merge_and_combat <- function(group_vector = NULL,use_combat = TRUE) {
  
  #------------------ STEP 1 — Select Directory ------------------
  data_dir <- tclvalue(
    tkchooseDirectory(title = "Select folder containing gene-level XLSX files")
  )
  if (data_dir == "") stop("No folder selected. Aborting.")
  
  cat("📂 Selected directory:", data_dir, "\n")
  
  #------------------ STEP 2 — Locate XLSX Files ------------------
  files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) < 2)
    stop("Need at least two gene-matrix XLSX files for ComBat.")
  
  cat("Found", length(files), "files:\n")
  print(basename(files))
  
  #------------------ STEP 3 — Read and Clean Matrices ------------------
  exprs_list <- lapply(files, function(f) {
    df <- as.data.frame(read.xlsx(f, rowNames = TRUE))
    return(df)
  })
  
  # Name the list by file base (without suffix)
  names(exprs_list) <- sub("_RMA_matrix\\.xlsx$", "", basename(files))
  
  #------------------ STEP 4 — Find Common Genes ------------------
  common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
  cat("Number of common genes:", length(common_genes), "\n")
  
  exprs_list_common <- lapply(exprs_list, function(x) x[common_genes, , drop = FALSE])
  
  #------------------ STEP 5 — merge ------------------
  combined_exprs <- do.call(cbind, exprs_list_common)
  # Remove any prefix before GSM
  colnames(combined_exprs) <- sub(".*(GSM\\d+)$", "\\1", colnames(combined_exprs))
  
  
  cat("Combined dataset dimensions:", dim(combined_exprs), "\n")
  
  #------------------ STEP 6 — Filter Columns to Match Group Vector ------------------
  # Keep only samples present in group_vector (remove any NA samples)
  keep_samples <- names(group_vector)[names(group_vector) %in% colnames(combined_exprs)]
  combined_exprs <- combined_exprs[, keep_samples]
  group_vector <- group_vector[keep_samples]
  
  # Make sure order matches
  combined_exprs <- combined_exprs[, names(group_vector)]
  
  #------------------ STEP 7 — Create Batch Vector ------------------
  batch <- sapply(colnames(combined_exprs), function(gsm) {
    # Find which dataset this GSM belongs to
    gse <- names(exprs_list_common)[sapply(exprs_list_common, function(mat) gsm %in% sub(".*(GSM\\d+)$", "\\1", colnames(mat)))]
    gse
  })
  batch <- factor(batch)
  
  cat("Batch composition:\n")
  print(table(batch))
  
  #------------------ STEP 8 — Apply ComBat Correction ------------------
  cat("\nApplying ComBat batch correction with biological group preserved...\n")
  mod <- model.matrix(~ group_vector)
  combat_exprs <- ComBat(
    dat = as.matrix(combined_exprs),
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  cat("✅ Batch correction complete.\n")
  
  #------------------ STEP 9 — Save and Return ------------------
  out_file <- file.path(data_dir, "ComBat_final_RMA_matrix.xlsx")
  write.xlsx(combat_exprs, out_file, rowNames = TRUE)
  cat("📁 Corrected expression matrix saved to:", out_file, "\n")
  
  # Set plotting area for side-by-side
  par(mfrow = c(1, 2), mar = c(10, 4, 4, 2))  # rotate x labels for readability
  
  # Boxplot before ComBat
  boxplot(combined_exprs,
          outline = FALSE,
          col = "lightcoral",
          las = 2,
          main = "Before ComBat",
          ylab = "Expression")
  
  # Boxplot after ComBat
  boxplot(combat_exprs,
          outline = FALSE,
          col = "lightseagreen",
          las = 2,
          main = "After ComBat",
          ylab = "Expression")
  
  # Reset plotting layout
  par(mfrow = c(1, 1))
  
  return(combat_exprs)
}


#=========================================================
# FUNCTION: run_limma_deg()
#=========================================================
run_limma_deg <- function(expr_mat,
                          group_vector,
                          group_case,
                          group_control,
                          logfc_cutoff = 1,
                          p_cutoff = 0.05,
                          adjust_method = "BH") {
  
  # -----------------------------
  # 1. Sanity checks
  # -----------------------------
  if (is.null(colnames(expr_mat)))
    stop("Expression matrix must have sample IDs as column names.")
  
  if (is.null(names(group_vector)))
    stop("group_vector must be a named vector (names = sample IDs).")
  
  # -----------------------------
  # 2. Align group vector to expression matrix
  # -----------------------------
  common_samples <- intersect(colnames(expr_mat), names(group_vector))
  
  if (length(common_samples) < 2)
    stop("Too few matching samples between expression matrix and group_vector.")
  
  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  group_vector <- group_vector[common_samples]
  group_vector <- factor(group_vector)
  
  cat("Samples per group:\n")
  print(table(group_vector))
  
  # -----------------------------
  # 3. Design matrix
  # -----------------------------
  design <- model.matrix(~ 0 + group_vector)
  colnames(design) <- levels(group_vector)
  
  # -----------------------------
  # 4. Limma pipeline
  # -----------------------------
  fit <- lmFit(expr_mat, design)
  contrast_formula <- paste0(group_case, " - ", group_control)
  contrast <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # -----------------------------
  # 5. Extract results
  # -----------------------------
  deg_results <- topTable(
    fit2,
    number = Inf,
    adjust.method = adjust_method,
    sort.by = "P"
  )
  
  # -----------------------------
  # 6. Identify significant DEGs
  # -----------------------------
  degs <- deg_results[
    abs(deg_results$logFC) >= logfc_cutoff &
      deg_results$adj.P.Val < p_cutoff,
  ]
  
  # -----------------------------
  # 7. Split into up and downregulated
  # -----------------------------
  up_degs <- degs[degs$logFC > 0, ]
  down_degs <- degs[degs$logFC < 0, ]
  
  cat("Total genes tested:", nrow(deg_results), "\n")
  cat("Significant DEGs:", nrow(degs), "\n")
  cat("Upregulated genes:", nrow(up_degs), "\n")
  cat("Downregulated genes:", nrow(down_degs), "\n")
  
  return(list(
    deg_results = deg_results,
    degs = degs,
    up_degs = up_degs,
    down_degs = down_degs,
    fit = fit2,
    design = design
  ))
}

#=========================================================
# FUNCTION: run_wgcna()
#=========================================================
run_wgcna <- function(
    datExpr,
    group_vector,
    trait_case = "RA",
    mad_fraction = 0.25,
    softPower = 8,
    minModuleSize = 30,
    mergeCutHeight = 0.25
) {
  
  enableWGCNAThreads()
  options(stringsAsFactors = FALSE)

  # -------------------------------
  # 2. Quality control
  # -------------------------------
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  # -------------------------------
  # 3. MAD filtering (Figure only)
  # -------------------------------
  geneMAD <- apply(datExpr, 2, mad)
  
  topGenes <- names(sort(geneMAD, decreasing = TRUE))[
    1:floor(mad_fraction * ncol(datExpr))
  ]
  
  datExprMAD <- datExpr[, topGenes]
  
  # -------------------------------
  # 4–8. Network construction
  # -------------------------------
  adjacency <- adjacency(
    datExprMAD,
    power = softPower,
    type = "signed"
  )
  
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize
  )
  
  dynamicColors <- labels2colors(dynamicMods)
  
  # -------------------------------
  # 9–10. Module eigengenes & merge
  # -------------------------------
  MEList <- moduleEigengenes(datExprMAD, colors = dynamicColors)
  MEs <- orderMEs(MEList$eigengenes)
  
  merge <- mergeCloseModules(
    datExprMAD,
    dynamicColors,
    cutHeight = mergeCutHeight,
    verbose = 3
  )
  
  mergedColors <- merge$colors
  mergedMEs <- orderMEs(merge$newMEs)
  
  # -------------------------------
  # 11–13. Module–trait
  # -------------------------------
  clinicalTraits <- data.frame(
    RA = as.numeric(group_vector == "RA"),
    NC = as.numeric(group_vector == "NC")
  )
  
  moduleTraitCor <- cor(mergedMEs, clinicalTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor,
    nSamples = nrow(datExprMAD)
  )
  
  sigModule <- rownames(moduleTraitCor)[
    which.max(abs(moduleTraitCor[, 1]))
  ]
  
  sigColor <- gsub("ME", "", sigModule)
  
  # -------------------------------
  # 14–16. GS & MM
  # -------------------------------
  GS <- as.data.frame(cor(datExprMAD, clinicalTraits, use = "p"))
  colnames(GS) <- "GS"
  
  GS_pvalue <- corPvalueStudent(
    as.matrix(GS),
    nSamples = nrow(datExprMAD)
  )
  
  MM <- as.data.frame(cor(datExprMAD, mergedMEs, use = "p"))
  MM_pvalue <- corPvalueStudent(
    as.matrix(MM),
    nSamples = nrow(datExprMAD)
  )
  
  geneInfo <- data.frame(
    Gene = colnames(datExprMAD),
    Module = mergedColors,
    GS = GS$GS,
    GS_pvalue = GS_pvalue,
    MM = abs(MM[, sigModule]),
    MM_pvalue = MM_pvalue[, sigModule]
  )
  
  geneInfo_mod <- geneInfo[geneInfo$Module == sigColor, ]
  geneInfo_mod <- geneInfo_mod[
    order(-geneInfo_mod$MM, -geneInfo_mod$GS),
  ]
  
  hubGenes <- geneInfo_mod$Gene
  
  #--------------------------------
  # Plot everything
  #--------------------------------
  #WGCNA Plotting
  dev.new()
  plot_soft_threshold(sft_data)
  dev.new()
  plot_gene_dendrogram(datExprMAD, softPower, dissTOM)
  dev.new()
  plot_module_dendrogram(geneTree, dynamicColors, mergedColors)
  dev.new()
  plot_module_trait_heatmap(moduleTraitCor, moduleTraitPvalue, mergedMEs, trait_name = c("RA", "NC"))
  dev.new()
  plot_GS_MM(geneInfo_mod, sigModule, sigColor)
  
  # -------------------------------
  # Return everything
  # -------------------------------
  return(list(
    datExprMAD = datExprMAD,
    geneTree = geneTree,
    dissTOM=dissTOM,
    dynamicColors = dynamicColors,
    mergedColors = mergedColors,
    mergedMEs = mergedMEs,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    sigModule = sigModule,
    sigColor = sigColor,
    geneInfo = geneInfo,
    geneInfo_mod = geneInfo_mod,
    hubGenes = hubGenes
  ))
}

#=========================================================
# FUNCTION: run_lasso_cv()
#=========================================================
run_lasso_cv <- function(expr_mat,group_vector,case_label = "RA",nfolds = 10,seed = 123) {
  
  set.seed(seed)
  
  # Binary outcome
  y <- as.numeric(group_vector == case_label)
  
  cv.lasso <- cv.glmnet(
    x = as.matrix(expr_mat),
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = nfolds
  )
  
  lambda_min <- cv.lasso$lambda.min
  
  lasso_model <- glmnet(
    x = as.matrix(expr_mat),
    y = y,
    family = "binomial",
    alpha = 1,
    lambda = lambda_min
  )
  
  coef_mat <- coef(lasso_model)
  
  lasso_genes <- rownames(coef_mat)[coef_mat[, 1] != 0]
  lasso_genes <- setdiff(lasso_genes, "(Intercept)")
  
  return(list(
    selected_genes = lasso_genes,
    lambda_min = lambda_min,
    cv_model = cv.lasso
  ))
}

#=========================================================
# FUNCTION: run_svm_rfe()
#=========================================================
run_svm_rfe <- function(expr_mat,group_vector,folds = 5,max_features = 30,seed = 123) {

  set.seed(seed)
  
  ctrl <- rfeControl(
    functions = caretFuncs,
    method = "cv",
    number = folds
  )
  
  sizes <- seq(1, min(max_features, ncol(expr_mat)), by = 1)
  
  svm_profile <- rfe(
    x = expr_mat,
    y = factor(group_vector),
    sizes = sizes,
    rfeControl = ctrl,
    method = "svmLinear"
  )
  
  svm_genes <- predictors(svm_profile)
  
  return(list(
    selected_genes = svm_genes,
    rfe_model = svm_profile
  ))
}

#=========================================================
# FUNCTION: run_random_forest()
#=========================================================
run_random_forest <- function(expr_mat,group_vector,ntree = 500,importance_quantile = 0.75,seed = 123) {

  set.seed(seed)
  
  rf_model <- randomForest(
    x = expr_mat,
    y = factor(group_vector),
    importance = TRUE,
    ntree = ntree
  )
  
  rf_importance <- importance(rf_model, type = 2)
  
  cutoff <- quantile(
    rf_importance[, "MeanDecreaseGini"],
    importance_quantile
  )
  
  rf_genes <- rownames(rf_importance)[
    rf_importance[, "MeanDecreaseGini"] > cutoff
  ]
  
  return(list(
    selected_genes = rf_genes,
    rf_model = rf_model,
    importance = rf_importance
  ))
}

