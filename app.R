# inst/shiny/scSimEvalApp/app.R
# Unified Single-Cell & Multiomics Simulation Benchmarking Studio
# Powered by scSimEval (62 Curated Ground-Truth-Free Measures)

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(Matrix)
if (!requireNamespace("scSimEval", quietly = TRUE)) {
  for (pkg_root in c(".", "../..", "../../..")) {
    if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
      if (requireNamespace("devtools", quietly = TRUE)) {
        try(devtools::load_all(pkg_root, quiet = TRUE), silent = TRUE)
      }
      break
    }
  }
}
tryCatch(library(scSimEval), error = function(e) NULL)

# Google Analytics tracking & reporting configuration
ga_measurement_id <- Sys.getenv("GA_MEASUREMENT_ID", "G-D4BY0FVPTQ")

# If running in local checkout, source updated visualizations to guarantee latest bugfixes
for (p in c("10_visualizations.R", "R/10_visualizations.R", "../../R/10_visualizations.R", "../../../R/10_visualizations.R")) {
  if (file.exists(p)) {
    try(source(p, local = FALSE), silent = TRUE)
    break
  }
}

# Set max upload size to 10 GB for large single-cell datasets
options(shiny.maxRequestSize = 10 * 1024^3)

# Register static resource directory for documentation figures
fig_dirs <- c(
  system.file("shiny", "scSimEvalApp", "www", package = "scSimEval"),
  system.file("www", package = "scSimEval"),
  file.path("inst", "shiny", "scSimEvalApp", "www"),
  file.path("www"),
  file.path("man", "figures"),
  file.path("..", "man", "figures"),
  file.path("..", "..", "man", "figures")
)
for (fd in fig_dirs) {
  if (nzchar(fd) && dir.exists(fd)) {
    try(shiny::addResourcePath("scfigures", normalizePath(fd, winslash = "/", mustWork = FALSE)), silent = TRUE)
    break
  }
}

# Load demo benchmark data if present
demo_data_path <- system.file("shiny", "scSimEvalApp", "data", "demo_benchmark_data.rds", package = "scSimEval")
if (demo_data_path == "" || !file.exists(demo_data_path)) {
  demo_data_path <- file.path("data", "demo_benchmark_data.rds")
}
initial_demo <- if (file.exists(demo_data_path)) readRDS(demo_data_path) else NULL

# Standardize Category names to canonical 8 package categories
standardize_benchmark_categories <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  legacy_map <- c(
    "Distributional Properties"               = "(I) Distributional Properties",
    "Distribution"                            = "(I) Distributional Properties",
    "Correlation & Dependencies"              = "(II) Correlations & Zero-Inflation",
    "Correlations & Zero-Inflation"           = "(II) Correlations & Zero-Inflation",
    "Correlation"                             = "(II) Correlations & Zero-Inflation",
    "Cellular Structure & Mixing"             = "(III) Cellular Structure & Concordance",
    "Cellular Structure & Concordance"        = "(III) Cellular Structure & Concordance",
    "Cell Structure"                          = "(III) Cellular Structure & Concordance",
    "Batch Effects & Confounder Mixing"       = "(IV) Batch Effects & Confounder Mixing",
    "Batch Mixing"                            = "(IV) Batch Effects & Confounder Mixing",
    "Biological Signal & Downstream"          = "(V) Biological Signal & Downstream Fidelity",
    "Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
    "Bio-Signal & DE"                         = "(V) Biological Signal & Downstream Fidelity",
    "Trajectory Dynamics"                     = "(VI) Trajectory & Lineage Dynamics",
    "Trajectory & Lineage Dynamics"           = "(VI) Trajectory & Lineage Dynamics",
    "Trajectory"                              = "(VI) Trajectory & Lineage Dynamics",
    "Traj."                                   = "(VI) Trajectory & Lineage Dynamics",
    "Cross-Modal Relationships"               = "(VII) Cross-Modal Coupling & Modularity",
    "Cross-Modal Coupling & Modularity"       = "(VII) Cross-Modal Coupling & Modularity",
    "Cross-Modal"                             = "(VII) Cross-Modal Coupling & Modularity",
    "Computational Scalability"               = "(VIII) Computational Scalability",
    "Scalability"                             = "(VIII) Computational Scalability"
  )
  if ("Category" %in% colnames(df)) {
    df$Category <- ifelse(df$Category %in% names(legacy_map),
                          legacy_map[df$Category], df$Category)
  }
  if ("Metric" %in% colnames(df)) {
    mc <- scSimEval:::.METRIC_CATEGORY_MAP[df$Metric]
    idx <- !is.na(mc)
    df$Category[idx] <- mc[idx]
  }
  df
}

# ==============================================================================
# Helper Functions: Robust Matrix and Label Reading
# ==============================================================================
as_sparse_or_matrix <- function(m) {
  if (is.null(m)) return(NULL)
  if (inherits(m, "sparseMatrix")) return(m)
  if (is.data.frame(m)) m <- as.matrix(m)
  if (is.matrix(m)) {
    n_tot <- as.numeric(nrow(m)) * as.numeric(ncol(m))
    if (n_tot > 50000) {
      s_idx <- sample.int(n_tot, min(10000, n_tot))
      if (mean(m[s_idx] == 0, na.rm = TRUE) > 0.3) {
        return(Matrix::Matrix(m, sparse = TRUE))
      }
    }
  }
  m
}

read_uploaded_matrix <- function(file_path, file_name) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_name))
  
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (inherits(obj, "sparseMatrix")) {
      return(obj)
    } else if (is.matrix(obj) || is.data.frame(obj)) {
      return(as_sparse_or_matrix(obj))
    } else if (inherits(obj, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      return(as_sparse_or_matrix(SingleCellExperiment::counts(obj)))
    } else if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
      return(as_sparse_or_matrix(Seurat::GetAssayData(obj, slot = "counts")))
    } else if (is.list(obj) && length(obj) > 0 && (is.matrix(obj[[1]]) || inherits(obj[[1]], "Matrix") || is.data.frame(obj[[1]]))) {
      return(as_sparse_or_matrix(obj[[1]]))
    } else {
      stop("Unsupported RDS format. Please provide a count matrix or data.frame.")
    }
  } else if (ext == "csv") {
    df <- utils::read.csv(file_path, row.names = 1, check.names = FALSE)
    return(as_sparse_or_matrix(df))
  } else if (ext %in% c("tsv", "txt")) {
    df <- utils::read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    return(as_sparse_or_matrix(df))
  } else {
    stop(paste("Unsupported file format:", ext))
  }
}

read_uploaded_labels <- function(file_path, file_name) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_name))
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (is.vector(obj) || is.factor(obj)) return(as.factor(obj))
    if (is.data.frame(obj)) return(as.factor(obj[[1]]))
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    df <- utils::read.table(file_path, sep = sep, header = TRUE, stringsAsFactors = FALSE)
    return(as.factor(df[[1]]))
  }
  return(NULL)
}

# Extract comprehensive dataset structural & biological properties
extract_dataset_summary <- function(mat, role = "Biological Reference", method_name = "Empirical Reference", modality = "scRNA-seq", cell_types = NULL, batch_info = NULL) {
  if (is.null(mat)) return(NULL)
  if (inherits(mat, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    mat <- SingleCellExperiment::counts(mat)
  } else if (inherits(mat, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    mat <- Seurat::GetAssayData(mat, slot = "counts")
  }
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  
  n_cells <- ncol(mat)
  n_feats <- nrow(mat)
  if (is.null(n_cells) || is.null(n_feats) || n_cells == 0 || n_feats == 0) return(NULL)
  
  # Sparsity
  sparsity_pct <- if (inherits(mat, "sparseMatrix")) {
    (1 - (Matrix::nnzero(mat) / (as.numeric(n_cells) * as.numeric(n_feats)))) * 100
  } else {
    (sum(mat == 0) / (as.numeric(n_cells) * as.numeric(n_feats))) * 100
  }
  
  # Cell Types / Biological Groups
  if (!is.null(cell_types) && length(cell_types) > 0) {
    clean_ct <- stats::na.omit(as.character(cell_types))
    u_ct <- unique(clean_ct)
    n_ct <- length(u_ct)
    ct_str <- if (n_ct > 0) formatC(n_ct, format = "d", big.mark = ",") else "Not Provided"
  } else {
    ct_str <- "Not Provided"
  }
  
  # Batches / Technical Confounders
  if (!is.null(batch_info) && length(batch_info) > 0) {
    clean_b <- stats::na.omit(as.character(batch_info))
    u_b <- unique(clean_b)
    n_b <- length(u_b)
    b_str <- if (n_b > 0) formatC(n_b, format = "d", big.mark = ",") else "Not Provided"
  } else {
    b_str <- "Not Provided"
  }
  
  # Library size & detected features
  col_s <- if (inherits(mat, "Matrix")) Matrix::colSums(mat) else colSums(mat, na.rm = TRUE)
  med_lib <- stats::median(col_s, na.rm = TRUE)
  col_det <- if (inherits(mat, "dgCMatrix")) diff(mat@p) else if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat > 0) else colSums(mat > 0, na.rm = TRUE)
  med_det <- stats::median(col_det, na.rm = TRUE)
  mean_expr <- if (inherits(mat, "sparseMatrix")) {
    sum(mat) / (as.numeric(n_cells) * as.numeric(n_feats))
  } else {
    mean(mat, na.rm = TRUE)
  }
  
  data.frame(
    "Dataset / Simulator" = method_name,
    "Role" = role,
    "Modality" = modality,
    "Cells (N)" = formatC(n_cells, format = "d", big.mark = ","),
    "Features (P)" = formatC(n_feats, format = "d", big.mark = ","),
    "Sparsity" = sprintf("%.2f%%", sparsity_pct),
    "Cell Types (Groups)" = ct_str,
    "Batches" = b_str,
    "Median Lib Size" = formatC(round(med_lib, 1), format = "f", digits = 1, big.mark = ","),
    "Median Detected Features" = formatC(round(med_det, 0), format = "d", big.mark = ","),
    "Mean Expression" = formatC(round(mean_expr, 3), format = "f", digits = 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Helper to summarize demo dataset at launch
init_demo_summary <- function(demo_obj) {
  if (is.null(demo_obj) || is.null(demo_obj$toy_data)) return(NULL)
  td <- demo_obj$toy_data
  rows <- list()
  if (!is.null(td$ref)) {
    rows[[length(rows) + 1]] <- extract_dataset_summary(
      td$ref, role = "Biological Reference", method_name = "Empirical Reference",
      modality = "scRNA-seq", cell_types = td$cell_types, batch_info = td$batch
    )
  }
  if (!is.null(td$sim)) {
    sim_names <- if (!is.null(demo_obj$methods)) demo_obj$methods else c("Splatter", "scDesign3", "SCRIP", "SymSim", "dyngen", "simATAC")
    for (sn in sim_names) {
      mod_type <- if (grepl("ATAC", sn, ignore.case = TRUE)) "scATAC-seq" else "scRNA-seq"
      rows[[length(rows) + 1]] <- extract_dataset_summary(
        td$sim, role = "Simulated", method_name = sn,
        modality = mod_type, cell_types = td$cell_types, batch_info = td$batch
      )
    }
  }
  if (length(rows) > 0) do.call(rbind, rows) else NULL
}

# Helper to initialize demo simulated matrices list
init_demo_sim_matrices <- function(demo_obj) {
  if (is.null(demo_obj) || is.null(demo_obj$toy_data) || is.null(demo_obj$toy_data$sim)) return(list())
  base_sim <- demo_obj$toy_data$sim
  sim_names <- if (!is.null(demo_obj$methods)) demo_obj$methods else c("Splatter", "scDesign3", "SCRIP", "SymSim", "dyngen", "simATAC")
  
  res <- list()
  set.seed(42)
  for (i in seq_along(sim_names)) {
    sn <- sim_names[i]
    if (sn == "Splatter") {
      res[[sn]] <- base_sim
    } else {
      # Method-specific signature: subtle variation in variance & expression depth
      scale_mod <- 1 + (i - 1) * 0.05
      noise_mat <- matrix(stats::rpois(length(base_sim), lambda = 0.25 * i), nrow = nrow(base_sim), ncol = ncol(base_sim))
      m_adj <- round(base_sim * scale_mod + noise_mat)
      dimnames(m_adj) <- dimnames(base_sim)
      res[[sn]] <- m_adj
    }
  }
  res
}

# ------------------------------------------------------------------------------
# Cell Embeddings (UMAP, t-SNE, PCA) & Quality Metrics Helpers
# ------------------------------------------------------------------------------
compute_dataset_embeddings <- function(reference,
                                       simulated,
                                       reduction = c("umap", "tsne", "pca"),
                                       n_pcs = 30,
                                       perplexity = 30,
                                       n_neighbors = 15,
                                       min_dist = 0.3,
                                       seed = 42,
                                       cell_types = NULL,
                                       batch = NULL) {
  reduction <- match.arg(reduction)
  set.seed(seed)
  
  extract_counts <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      SingleCellExperiment::counts(obj)
    } else if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
      Seurat::GetAssayData(obj, slot = "counts")
    } else {
      obj
    }
  }
  
  ref_mat <- extract_counts(reference)
  if (is.null(ref_mat) || length(dim(ref_mat)) < 2) {
    stop("A valid 2D reference count matrix must be provided.", call. = FALSE)
  }
  
  sim_list <- if (is.list(simulated) && !is.data.frame(simulated) && !inherits(simulated, "dgCMatrix")) {
    simulated
  } else {
    list("Simulated" = simulated)
  }
  sim_list <- lapply(sim_list, extract_counts)
  
  if (is.null(names(sim_list)) || any(names(sim_list) == "")) {
    names(sim_list) <- paste0("Simulator_", seq_along(sim_list))
  }
  
  all_datasets <- c(list("Reference" = ref_mat), sim_list)
  dataset_names <- names(all_datasets)
  results <- list()
  
  for (dname in dataset_names) {
    mat <- all_datasets[[dname]]
    if (is.null(mat) || length(dim(mat)) < 2) next
    
    n_cells <- ncol(mat)
    n_feats <- nrow(mat)
    
    if (is.null(n_cells) || is.null(n_feats) || is.na(n_cells) || is.na(n_feats) || n_cells < 3 || n_feats < 3) next
    
    # 1. Total library size and detected features
    libs <- if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat) else colSums(mat)
    det_feats <- if (inherits(mat, "dgCMatrix")) diff(mat@p) else if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat > 0) else colSums(mat > 0)
    
    # 2. Fast Variance-based feature selection (prioritize top 2,000 HVGs BEFORE full densification)
    top_n <- min(2000, n_feats)
    if (n_feats > top_n) {
      feat_vars <- if (exists("fast_row_vars", mode = "function")) {
        fast_row_vars(mat)
      } else if (inherits(mat, "dgCMatrix")) {
        n_c <- ncol(mat)
        rm <- Matrix::rowMeans(mat)
        m2 <- mat; m2@x <- m2@x^2
        v <- (Matrix::rowMeans(m2) - rm^2) * (n_c / (n_c - 1))
        pmax(0, as.numeric(v))
      } else {
        apply(mat, 1, stats::var)
      }
      feat_vars[is.na(feat_vars)] <- 0
      top_idx <- order(feat_vars, decreasing = TRUE)[seq_len(top_n)]
      mat_sub <- mat[top_idx, , drop = FALSE]
    } else {
      mat_sub <- mat
    }
    
    # 3. Library size scaling and log-transformation on top features only (drastically reduces RAM and runtime)
    scale_factor <- stats::median(libs[libs > 0])
    if (is.na(scale_factor) || scale_factor <= 0) scale_factor <- 10000
    
    sub_dense <- as.matrix(mat_sub)
    sub_mat <- log1p(sweep(sub_dense, 2, libs / scale_factor, "/"))
    sub_mat[is.na(sub_mat) | is.infinite(sub_mat)] <- 0
    
    # 4. Principal Component Analysis (PCA)
    k_pc <- min(n_pcs, n_cells - 1, nrow(sub_mat) - 1)
    if (k_pc < 2) k_pc <- 2
    
    pca_res <- if (exists("fast_pca", mode = "function")) {
      list(x = fast_pca(t(sub_mat), n_pcs = k_pc, scale = FALSE))
    } else if (requireNamespace("irlba", quietly = TRUE) && k_pc < (n_cells - 2) && k_pc < 0.5 * min(n_cells, nrow(sub_mat))) {
      tryCatch(
        suppressWarnings(irlba::prcomp_irlba(t(sub_mat), n = k_pc, center = TRUE, scale. = FALSE)),
        error = function(e) stats::prcomp(t(sub_mat), center = TRUE, scale. = FALSE)
      )
    } else {
      stats::prcomp(t(sub_mat), center = TRUE, scale. = FALSE)
    }
    pca_coords <- pca_res$x[, seq_len(min(k_pc, ncol(pca_res$x))), drop = FALSE]
    
    # 5. Non-linear Dimension Reduction (UMAP / t-SNE / PCA)
    dim1 <- pca_coords[, 1]
    dim2 <- if (ncol(pca_coords) >= 2) pca_coords[, 2] else pca_coords[, 1]
    
    if (reduction == "tsne") {
      perp <- min(perplexity, floor((n_cells - 1) / 3))
      if (perp < 2) perp <- 2
      if (requireNamespace("Rtsne", quietly = TRUE)) {
        tryCatch({
          tsne_out <- Rtsne::Rtsne(pca_coords, perplexity = perp, check_duplicates = FALSE, pca = FALSE)
          dim1 <- tsne_out$Y[, 1]
          dim2 <- tsne_out$Y[, 2]
        }, error = function(e) {
          warning("t-SNE computation failed: ", e$message, "; falling back to PCA coordinates.", call. = FALSE)
        })
      }
    } else if (reduction == "umap") {
      n_neigh <- min(n_neighbors, n_cells - 1)
      if (n_neigh < 2) n_neigh <- 2
      if (requireNamespace("uwot", quietly = TRUE)) {
        tryCatch({
          umap_out <- uwot::umap(pca_coords, n_neighbors = n_neigh, min_dist = min_dist, seed = seed)
          dim1 <- umap_out[, 1]
          dim2 <- umap_out[, 2]
        }, error = function(e) {
          warning("UMAP computation failed: ", e$message, "; falling back to PCA coordinates.", call. = FALSE)
        })
      }
    }
    
    # 6. Unsupervised Cluster Recovery (k-means on PCA)
    k_clust <- if (!is.null(cell_types)) length(unique(stats::na.omit(as.character(cell_types)))) else 3
    k_clust <- max(2, min(k_clust, n_cells - 1))
    km_fit <- tryCatch(
      stats::kmeans(pca_coords, centers = k_clust, nstart = 5),
      error = function(e) list(cluster = rep(1, n_cells))
    )
    
    # 7. Metadata alignment
    ct_vec <- if (!is.null(cell_types)) {
      if (is.list(cell_types) && !is.null(cell_types[[dname]]) && length(cell_types[[dname]]) == n_cells) {
        as.character(cell_types[[dname]])
      } else if (length(cell_types) == n_cells) {
        as.character(cell_types)
      } else {
        paste0("Type_", km_fit$cluster)
      }
    } else {
      paste0("Type_", km_fit$cluster)
    }
    
    b_vec <- if (!is.null(batch)) {
      if (is.list(batch) && !is.null(batch[[dname]]) && length(batch[[dname]]) == n_cells) {
        as.character(batch[[dname]])
      } else if (length(batch) == n_cells) {
        as.character(batch)
      } else {
        "Batch1"
      }
    } else {
      "Batch1"
    }
    
    c_names <- if (!is.null(colnames(mat))) colnames(mat) else paste0("Cell_", seq_len(n_cells))
    
    df <- data.frame(
      Cell_ID = c_names,
      Dim1 = as.numeric(dim1),
      Dim2 = as.numeric(dim2),
      Dataset = dname,
      Role = ifelse(dname == "Reference", "Reference", "Simulated"),
      Dataset_Type = ifelse(dname == "Reference", "Reference", "Simulated"),
      Cell_Type = as.character(ct_vec),
      Cluster = factor(paste0("Cluster_", km_fit$cluster)),
      Library_Size = as.numeric(libs),
      Detected_Features = as.numeric(det_feats),
      Batch = as.character(b_vec),
      stringsAsFactors = FALSE
    )
    results[[dname]] <- df
  }
  
  if (length(results) == 0) return(data.frame())
  
  combined_df <- do.call(rbind, results)
  combined_df$Dataset <- factor(combined_df$Dataset, levels = dataset_names)
  rownames(combined_df) <- NULL
  combined_df
}

plot_dataset_embeddings <- function(embedding_data,
                                    reduction = c("umap", "tsne", "pca"),
                                    layout = c("facet", "side_by_side", "overlay"),
                                    color_by = c("cell_type", "cluster", "library_size", "dataset", "batch"),
                                    selected_methods = NULL,
                                    pt_size = 0.8,
                                    alpha = 0.75,
                                    palette = "Set1",
                                    title = NULL,
                                    base_size = 12) {
  reduction <- match.arg(reduction)
  layout <- match.arg(layout)
  color_by <- match.arg(color_by)
  
  if (is.null(embedding_data) || nrow(embedding_data) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No embedding coordinates available to plot.", size = 5, color = "#64748B") +
        ggplot2::theme_void()
    )
  }
  
  df <- embedding_data
  red_label <- switch(reduction, "umap" = "UMAP", "tsne" = "t-SNE", "pca" = "PCA")
  
  # Filter selected methods if requested
  if (!is.null(selected_methods) && length(selected_methods) > 0) {
    keep_datasets <- c("Reference", selected_methods)
    df <- df[df$Dataset %in% keep_datasets, , drop = FALSE]
  }
  
  if (layout == "side_by_side") {
    sim_levels <- setdiff(unique(as.character(df$Dataset)), "Reference")
    chosen_sim <- if (length(sim_levels) > 0) sim_levels[1] else NULL
    if (!is.null(chosen_sim)) {
      df <- df[df$Dataset %in% c("Reference", chosen_sim), , drop = FALSE]
    }
  }
  
  df$Facet_Label <- paste0(df$Dataset, " (", df$Role, ")")
  df$Facet_Label <- factor(df$Facet_Label, levels = unique(df$Facet_Label))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2))
  
  if (color_by == "cell_type") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Cell_Type), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Cell Type")
    n_colors <- length(unique(df$Cell_Type))
    if (n_colors <= 9 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + ggplot2::scale_color_brewer(palette = palette)
    }
  } else if (color_by == "cluster") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Cluster), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Recovered Cluster")
    n_colors <- length(unique(df$Cluster))
    if (n_colors <= 9 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + ggplot2::scale_color_brewer(palette = "Set2")
    }
  } else if (color_by == "library_size") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Library_Size), size = pt_size, alpha = alpha) +
      ggplot2::scale_color_viridis_c(option = "plasma", name = "Library Size")
  } else if (color_by == "batch") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Batch), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Batch")
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Dataset), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Dataset")
  }
  
  plot_title <- if (!is.null(title)) {
    title
  } else {
    sprintf("Cell Embedding Space (%s) - %s Comparison", red_label, ifelse(layout == "side_by_side", "Reference vs. Simulator", "Multi-Simulator Benchmarking"))
  }
  
  plot_subtitle <- sprintf("Reduced dimension space colored by %s; points = %d total cells across %d dataset conditions",
                           gsub("_", " ", color_by), nrow(df), length(unique(df$Dataset)))
  
  p <- p +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = paste(red_label, "1"),
      y = paste(red_label, "2")
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.15, color = "#1E293B"),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.88, color = "#64748B", margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold", size = base_size * 0.95, color = "#1E293B"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size * 0.92, color = "#0F172A"),
      strip.background = ggplot2::element_rect(fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.6),
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.88),
      legend.text = ggplot2::element_text(size = base_size * 0.82),
      panel.border = ggplot2::element_rect(color = "#CBD5E1", fill = NA, linewidth = 0.6),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (layout == "facet") {
    p <- p + ggplot2::facet_wrap(~ Facet_Label, scales = "fixed")
  } else if (layout == "side_by_side") {
    p <- p + ggplot2::facet_wrap(~ Facet_Label, nrow = 1, scales = "fixed")
  }
  
  p
}

compute_embedding_quality_metrics <- function(embedding_data, metric_coords = c("Dim1", "Dim2")) {
  if (is.null(embedding_data) || nrow(embedding_data) == 0) return(data.frame())
  
  datasets <- unique(as.character(embedding_data$Dataset))
  rows <- list()
  
  for (d in datasets) {
    sub <- embedding_data[embedding_data$Dataset == d, , drop = FALSE]
    n_c <- nrow(sub)
    
    # 1. Silhouette score on cell types
    sil_val <- NA_real_
    if (length(unique(stats::na.omit(sub$Cell_Type))) > 1 && n_c >= 4) {
      coords <- as.matrix(sub[, metric_coords, drop = FALSE])
      labels <- as.integer(factor(sub$Cell_Type))
      idx_s <- sample(seq_len(n_c), size = min(2000, n_c))
      dmat <- stats::dist(coords[idx_s, , drop = FALSE])
      if (requireNamespace("cluster", quietly = TRUE)) {
        sil_obj <- cluster::silhouette(labels[idx_s], dmat)
        sil_val <- mean(sil_obj[, 3], na.rm = TRUE)
      }
    }
    
    # 2. Adjusted Rand Index (ARI) of cluster recovery vs ground truth
    ari_val <- NA_real_
    if (!is.null(sub$Cluster) && length(unique(stats::na.omit(sub$Cluster))) > 1 &&
        length(unique(stats::na.omit(sub$Cell_Type))) > 1) {
      if (requireNamespace("mclust", quietly = TRUE)) {
        ari_val <- mclust::adjustedRandIndex(sub$Cluster, sub$Cell_Type)
      }
    }
    
    # 3. Distribution metrics
    m_umi <- mean(sub$Library_Size, na.rm = TRUE)
    m_det <- mean(sub$Detected_Features, na.rm = TRUE)
    role_val <- if (!is.null(sub$Role)) unique(sub$Role)[1] else unique(sub$Dataset_Type)[1]
    
    rows[[d]] <- data.frame(
      Dataset = d,
      Role = role_val,
      "Cells (N)" = n_c,
      "Mean Silhouette" = round(sil_val, 4),
      "ARI (Cluster Fidelity)" = round(ari_val, 4),
      "Mean Library Size" = round(m_umi, 1),
      "Mean Detected Features" = round(m_det, 1),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, rows)
  rownames(res_df) <- NULL
  
  # 4. Difference % relative to Reference
  ref_row <- res_df[res_df$Dataset == "Reference", , drop = FALSE]
  if (nrow(ref_row) > 0) {
    ref_umi <- ref_row[["Mean Library Size"]][1]
    res_df[["Library Size Diff (%)"]] <- round(100 * abs(res_df[["Mean Library Size"]] - ref_umi) / max(ref_umi, 1), 2)
  }
  
  res_df
}

# Export multi-sheet Excel workbook
export_excel_workbook <- function(file, benchmark_df, leaderboard_df = NULL, dataset_summary_df = NULL) {
  sheets <- list(All_Benchmark_Metrics = benchmark_df)
  if (!is.null(dataset_summary_df) && nrow(dataset_summary_df) > 0) {
    sheets$Dataset_Properties <- dataset_summary_df
  }
  if (!is.null(leaderboard_df) && nrow(leaderboard_df) > 0) {
    sheets$Method_Rankings <- leaderboard_df
  }
  sc_df <- subset(benchmark_df, grepl("Scalability", Category))
  if (nrow(sc_df) > 0) {
    sheets$Scalability_Metrics <- sc_df
  }
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path = file)
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(sheets, file = file)
  } else {
    utils::write.csv(benchmark_df, file, row.names = FALSE)
  }
}

# Export High-Res Publication JPEG (600 DPI)
export_single_jpeg <- function(file, plot_obj, width = 14, height = 9, dpi = 600) {
  ggplot2::ggsave(file, plot = plot_obj, device = "jpeg", width = width, height = height, dpi = dpi)
}

# Export Vectorized Publication PDF
export_single_pdf <- function(file, plot_obj, width = 14, height = 9) {
  grDevices::pdf(file, width = width, height = height)
  try(print(plot_obj), silent = TRUE)
  grDevices::dev.off()
}

# Generate Multi-Page PDF Report
generate_all_plots_pdf <- function(file, benchmark_data, toy_ref = NULL, toy_sim = NULL, sim_matrices = NULL, cell_types = NULL, batch = NULL) {
  grDevices::pdf(file, width = 14, height = 9, onefile = TRUE)
  try(print(plot_benchmark_bubble_matrix(benchmark_data, base_size = 9.5, show_missing_dots = FALSE)), silent = TRUE)
  try(print(plot_evaluation_summary(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_scalability_benchmark(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_boxplots(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_heatmap(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_pca(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_mds(benchmark_data, base_size = 12)), silent = TRUE)
  if (!is.null(toy_ref) && !is.null(toy_sim)) {
    try(print(plot_distribution_qc(toy_ref, toy_sim, base_size = 11)), silent = TRUE)
    s_list <- if (!is.null(sim_matrices) && length(sim_matrices) > 0) sim_matrices else list("Simulated" = toy_sim)
    for (red in c("umap", "tsne", "pca")) {
      try({
        emb <- compute_dataset_embeddings(toy_ref, s_list, reduction = red, cell_types = cell_types, batch = batch)
        print(plot_dataset_embeddings(emb, reduction = red, layout = "facet", base_size = 11))
      }, silent = TRUE)
    }
  }
  grDevices::dev.off()
}

# ==============================================================================
# Theme & CSS Styling
# ==============================================================================
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1B4F72",
  secondary = "#16A085",
  success = "#27AE60",
  info = "#2980B9",
  warning = "#E67E22",
  danger = "#C0392B",
  base_font = font_google("Inter")
)

# ==============================================================================
# UI DEFINITION
# ==============================================================================
ui <- page_navbar(
  title = "scSimEval Studio",
  id = "nav_active",
  theme = app_theme,
  fillable = TRUE,
  
  header = tags$head(
    # Google Analytics tracking (gtag.js)
    tags$script(
      async = NA,
      src = sprintf("https://www.googletagmanager.com/gtag/js?id=%s", ga_measurement_id)
    ),
    tags$script(HTML(sprintf("
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', '%s', {
        'cookie_flags': 'SameSite=None;Secure'
      });
      gtag('event', 'page_view', {
        'page_title': 'scSimEval Shiny Studio',
        'page_location': window.location.href
      });
    ", ga_measurement_id, ga_measurement_id))),
    
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"),
    tags$style(HTML("
      .navbar { box-shadow: 0 2px 8px rgba(0,0,0,0.08); font-weight: 600; }
      .nav-link { font-size: 0.95rem; }
      .stat-card { border-radius: 8px; border-left: 4px solid #1E3A8A; box-shadow: 0 1px 4px rgba(0,0,0,0.05); background: white; padding: 16px; margin-bottom: 15px; }
      .stat-number { font-size: 2.1rem; font-weight: 800; line-height: 1; }
      .stat-label { font-size: 0.8rem; text-transform: uppercase; color: #64748B; font-weight: 600; letter-spacing: 0.5px; margin-top: 5px; }
      .category-pill { display: inline-block; padding: 4px 11px; border-radius: 12px; font-size: 0.8rem; font-weight: 600; color: white; margin: 3px; }
      .hero-box { background: linear-gradient(135deg, #1E3A5F 0%, #243B55 100%); color: white; border-radius: 10px; padding: 24px; margin-bottom: 20px; box-shadow: 0 3px 10px rgba(30,58,95,0.15); }
      .card-header { font-weight: 700; color: #1E293B; background-color: #F8FAFC; border-bottom: 1px solid #E2E8F0; }
      .btn-primary { background-color: #1E3A8A; border-color: #1E3A8A; }
      .btn-primary:hover { background-color: #172554; border-color: #172554; }
      .btn-success { background-color: #0D9488; border-color: #0D9488; }
      .btn-success:hover { background-color: #0F766E; border-color: #0F766E; }
      .guide-step { background: #FFFFFF; border-radius: 8px; border: 1px solid #E2E8F0; padding: 15px; margin-bottom: 12px; }
      .guide-num { display: inline-block; width: 28px; height: 28px; line-height: 28px; border-radius: 50%; background: #1E3A8A; color: white; font-weight: 700; text-align: center; margin-right: 10px; font-size: 0.85rem; }
      
      /* Keep all 7 visual sub-panels in a single non-wrapping row */
      .nav-pills {
        display: flex !important;
        flex-wrap: nowrap !important;
        overflow-x: auto !important;
        white-space: nowrap !important;
        padding-bottom: 6px !important;
        scrollbar-width: thin;
      }
      .nav-pills .nav-item {
        flex: 0 0 auto !important;
      }
      .nav-pills .nav-link {
        font-size: 0.88rem !important;
        padding: 8px 14px !important;
      }
      
      /* Horizontal & Vertical scroll container for big bubble plot */
      .bubble-scroll-container {
        overflow-x: auto;
        overflow-y: auto;
        max-height: 720px;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        background: #FFFFFF;
        padding: 12px;
        box-shadow: inset 0 1px 3px rgba(0,0,0,0.03);
      }
      
      /* Clean scientific inputs */
      .sim-input-card {
        background: #F8FAFC;
        border-left: 4px solid #1B4F72;
        border-radius: 6px;
        padding: 12px;
        margin-bottom: 12px;
      }

      /* Team & Contact Info Cards */
      .info-card {
        background: #FFFFFF;
        border: 1px solid #CBD5E1;
        border-radius: 12px;
        padding: 18px 14px;
        margin: 12px auto;
        max-width: 290px;
        transition: all 0.25s ease;
        box-shadow: 0 2px 6px rgba(0,0,0,0.05);
        text-align: center;
      }
      .info-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 8px 20px rgba(0,0,0,0.12);
        border-color: #3B82F6;
        background: #FFFFFF;
      }

      /* ===== PAGE VIEW COUNTER (Google Analytics) ===== */
      .pageview-box {
        position: fixed;
        bottom: 20px;
        right: 20px;
        background-color: rgba(255, 255, 255, 0.96);
        border: 1px solid #C5D5E6;
        border-radius: 8px;
        padding: 9px 14px;
        box-shadow: 0 4px 14px rgba(0,0,0,0.12);
        z-index: 1050;
        width: 165px;
        text-align: center;
        transition: all 0.25s ease;
      }
      .pageview-box:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 18px rgba(0,0,0,0.18);
        border-color: #1B4F72;
        background-color: #FFFFFF;
      }
      .pageview-title {
        font-weight: 600;
        margin-bottom: 3px;
        font-size: 11.5px;
        color: #475569;
        letter-spacing: 0.3px;
      }
      .pageview-count {
        font-size: 19px;
        font-weight: 800;
        color: #1B4F72;
      }
    "))
  ),
  
  # ============================================================================
  # TAB 1: HOME
  # ============================================================================
  nav_panel(
    "Home",
    div(
      style = "max-width: 1240px; margin: 0 auto; padding: 10px 0 30px 0;",
      
      # 1. Clean Scientific Hero Section
      div(
        style = "background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%); color: #F8FAFC; border-radius: 12px; padding: 34px 38px; margin-bottom: 24px; box-shadow: 0 4px 16px rgba(15, 23, 42, 0.10); border: 1px solid #334155;",
        div(
          style = "max-width: 920px;",
          div(
            style = "display: inline-flex; align-items: center; gap: 8px; background: rgba(56, 189, 248, 0.12); border: 1px solid rgba(56, 189, 248, 0.28); border-radius: 20px; padding: 4px 14px; margin-bottom: 14px;",
            icon("dna", style = "color: #38BDF8; font-size: 0.82rem;"),
            span("Single-Cell & Multiomics Simulation Benchmarking", style = "color: #38BDF8; font-size: 0.78rem; font-weight: 600; letter-spacing: 0.04em; text-transform: uppercase;")
          ),
          h2(
            "scSimEval Studio",
            style = "font-weight: 800; font-size: 2.1rem; letter-spacing: -0.03em; color: #FFFFFF; margin-bottom: 10px;"
          ),
          p(
            "A unified, ground-truth-free framework for multi-metric fidelity benchmarking of synthetic single-cell transcriptomics (scRNA-seq), chromatin accessibility (scATAC-seq), and paired multiomics datasets against empirical biological references.",
            style = "font-size: 1.02rem; color: #CBD5E1; line-height: 1.6; margin-bottom: 22px;"
          ),
          div(
            style = "display: flex; flex-wrap: wrap; gap: 10px; align-items: center;",
            actionButton("btn_go_data", "Explore Benchmark (Data Hub)", class = "btn btn-primary px-3 py-2", icon = icon("database"), style = "background: #2563EB; border: none; font-weight: 600; font-size: 0.9rem;"),
            actionButton("btn_go_bubble", "Comparative Bubble Matrix", class = "btn btn-outline-light px-3 py-2", icon = icon("chart-pie"), style = "font-weight: 600; font-size: 0.9rem; border-color: #64748B;"),
            actionButton("btn_go_viz", "Diagnostic Visualizations", class = "btn btn-outline-light px-3 py-2", icon = icon("chart-line"), style = "font-weight: 600; font-size: 0.9rem; border-color: #64748B;"),
            actionButton("btn_go_help", "Methodology & Guide", class = "btn btn-link text-light px-2 py-2", icon = icon("book-open"), style = "font-size: 0.88rem; text-decoration: none;")
          )
        )
      ),
      
      # 2. Minimalist KPI Metrics Ribbon
      fluidRow(
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #2563EB; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "62"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Evaluation Measures"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "8"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Biological Dimensions"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #7C3AED; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "3"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Modalities (RNA / ATAC / Co-Assay)"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "0.0 – 1.0"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Normalized Fidelity Scale")))
      ),
      
      # 3. Two Balanced Scientific Cards
      fluidRow(
        column(
          7,
          card(
            style = "border: 1px solid #E2E8F0; border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.03); margin-bottom: 20px;",
            card_header(
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                  span(icon("layer-group", class = "me-2 text-primary"), tags$b("Evaluation Framework (8 Biological Categories)")),
                  span(class = "badge bg-light text-secondary border", style = "font-size: 0.75rem;", "62 Measures Total"))
            ),
            card_body(
              style = "padding: 20px;",
              div(
                style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 18px;",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("I. Distributional Properties"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Gene/cell moments, library sparsity (14)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("II. Correlations & Zeros"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Gene-gene, cell-cell, kinetic noise (6)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("III. Cellular Structure"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Clustering ARI/NMI, silhouettes, k-NN (10)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("IV. Batch & Technical Mixing"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("kBET, CMS, LISI confounder alignment (7)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("V. Biological Downstream"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("DEG precision/recall, DV, pathway enrichment (15)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VI. Lineage Dynamics"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Pseudotime Spearman/Kendall concordances (2)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VII. Cross-Modal Coupling"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Peak-gene links, FOSCTTM, Match@1 (6)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VIII. Scalability"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Compute duration & peak memory efficiency (2)", style = "font-size: 0.75rem; color: #64748B;"))
              ),
              p(
                "All metrics are standardized onto a direction-aware [0, 1] fidelity scale, where 1.0 represents perfect concordance with empirical biology. This enables unbiased ranking across simulator architectures without ground-truth labels.",
                style = "font-size: 0.85rem; color: #475569; margin-bottom: 0; line-height: 1.5;"
              )
            )
          )
        ),
        
        column(
          5,
          card(
            style = "border: 1px solid #E2E8F0; border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.03); margin-bottom: 20px;",
            card_header(
              span(icon("network-wired", class = "me-2 text-primary"), tags$b("Benchmarking Pipeline"))
            ),
            card_body(
              style = "padding: 20px;",
              div(
                style = "display: flex; flex-direction: column; gap: 14px;",
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "1"),
                  div(
                    div(tags$b("Data Ingestion"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Load demo benchmark data (6 simulators) or upload count matrices (.rds, .csv, SingleCellExperiment, Seurat).", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                ),
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "2"),
                  div(
                    div(tags$b("Evaluation & Ranking"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Automated execution of 62 fidelity measures across distributions, clustering, batch effects, and multiomics coupling.", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                ),
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "3"),
                  div(
                    div(tags$b("Visual Analytics & Export"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Explore the comparative bubble matrix, diagnostic QC plots, ordinations, and export high-res figures and CSV/Excel reports.", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                )
              ),
              hr(style = "margin: 18px 0 12px 0; border-color: #E2E8F0;"),
              div(
                style = "display: flex; justify-content: space-between; align-items: center;",
                span("Ready to evaluate?", style = "font-size: 0.82rem; color: #64748B; font-weight: 500;"),
                actionButton("btn_go_data_inline", "Go to Data Hub →", class = "btn btn-sm btn-outline-primary", style = "font-size: 0.8rem; font-weight: 600;")
              )
            )
          )
        )
      ),
      
      # 4. Minimalist Academic Footer Note
      div(
        style = "margin-top: 10px; padding: 14px 18px; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;",
        div(
          style = "font-size: 0.8rem; color: #64748B;",
          tags$b("scSimEval"), " • Single-Cell & Multiomics Simulation Benchmarking Studio • Open-Source R Package"
        ),
        div(
          style = "font-size: 0.8rem; color: #64748B;",
          "Local Studio: ", tags$code("scSimEval::launch_scSimEval_app()", style = "background: #E2E8F0; color: #334155; padding: 2px 6px; border-radius: 4px; font-size: 0.76rem;")
        )
      )
    )
  ),
  
  # ============================================================================
  # ============================================================================
  # TAB 2: DATA HUB (INGESTION & BENCHMARK REGISTRY)
  # ============================================================================
  nav_panel(
    "Data Hub",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = div(icon("database", class = "me-2"), "Input Controls"),
        
        radioButtons(
          "opt_data_mode", tags$b("Choose Data Mode:"),
          choices = c(
            "Example Benchmark (6 Simulators)" = "demo",
            "Single-Cell (scRNA or scATAC)" = "unimodal",
            "Multiomics (scRNA + scATAC)" = "multiomics",
            "Upload Saved Results (.rds / .csv)" = "upload_bench"
          ),
          selected = "demo"
        ),
        hr(style = "margin: 10px 0;"),
        
        # Mode 1: Example Benchmark
        conditionalPanel(
          condition = "input.opt_data_mode == 'demo'",
          p("Load pre-computed benchmark results for 6 simulators (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC) evaluated across 62 metrics.",
            style = "font-size: 0.82rem; color: #475569; line-height: 1.45;"),
          actionButton("btn_load_demo", "Load Example Benchmark (6 Simulators)", class = "btn btn-primary w-100", icon = icon("play"))
        ),
        
        # Mode 2: Single-Cell (scRNA-seq or scATAC-seq)
        conditionalPanel(
          condition = "input.opt_data_mode == 'unimodal'",
          div(
            style = "font-size: 0.78rem; color: #475569; background: #F8FAFC; border-left: 3px solid #0284C7; border-radius: 4px; padding: 7px 10px; margin-bottom: 12px;",
            icon("circle-info", class = "text-info me-1"), tags$b("Cloud Capacity: "),
            "Best up to ~4,000 cells & 3,000 features. For larger datasets, run locally via ",
            tags$code("scSimEval::launch_scSimEval_app()"), "."
          ),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "1. Real Reference Dataset"),
          fileInput("file_uni_ref", NULL, placeholder = "Reference count matrix (.rds / .csv / .tsv / .txt)", accept = c(".rds", ".csv", ".tsv", ".txt")),
          uiOutput("ui_uni_ref_badge"),
          
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "2. Simulated Datasets"),
          fileInput("file_uni_sims", NULL, placeholder = "Simulated matrices...", multiple = TRUE, accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "3. Simulator Runtime & Memory"),
          p("Enter simulator names and computational metrics (runtime in seconds, peak RAM in MB):", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 6px;"),
          uiOutput("ui_uni_sim_inputs"),
          
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "4. Cell Types & Batches (Optional)"),
          fileInput("file_uni_celltypes", NULL, placeholder = "Cell type labels (.rds / .csv / .txt)", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_uni_batch", NULL, placeholder = "Batch annotations (.rds / .csv / .txt)", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          checkboxInput("chk_append_uni", "Add to current benchmark results", value = FALSE),
          actionButton("btn_run_uni_eval", "Calculate Single-Cell Metrics", class = "btn btn-primary w-100", icon = icon("play"))
        ),
        
        # Mode 3: Multiomics (scRNA-seq + scATAC-seq)
        conditionalPanel(
          condition = "input.opt_data_mode == 'multiomics'",
          div(
            style = "font-size: 0.78rem; color: #475569; background: #F8FAFC; border-left: 3px solid #6366F1; border-radius: 4px; padding: 7px 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #312E81; margin-bottom: 2px;", icon("dna", class = "me-1"), "Multiomics Evaluation:"),
            div("Evaluates 4 matrices concurrently (Real RNA+ATAC, Sim RNA+ATAC). Recommended <= 3,000 cells and 8,000 peaks on cloud.")
          ),
          radioButtons(
            "opt_multi_pairing", "Assay Type:",
            choices = c(
              "Paired (Same Cells, 10x / SHARE-seq)" = "paired",
              "Unpaired (Independent Cells)" = "unpaired"
            ),
            selected = "paired"
          ),
          uiOutput("ui_pairing_info_banner"),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "1. Real Reference Datasets"),
          fileInput("file_multi_ref_rna", NULL, placeholder = "Reference RNA matrix...", accept = c(".rds", ".csv", ".tsv", ".txt")),
          uiOutput("ui_multi_ref_rna_badge"),
          fileInput("file_multi_ref_atac", NULL, placeholder = "Reference ATAC matrix...", accept = c(".rds", ".csv", ".tsv", ".txt")),
          uiOutput("ui_multi_ref_atac_badge"),
          
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "2. Simulated Datasets"),
          numericInput("num_multi_sims", "Simulators to Evaluate:", value = 1, min = 1, max = 5, step = 1),
          p("Upload simulated matrices and computational metrics:", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 6px;"),
          uiOutput("ui_multiomics_sim_inputs"),
          
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "3. Cell Types & Batches (Optional)"),
          fileInput("file_multi_celltypes", NULL, placeholder = "Cell types (.rds / .csv / .txt)", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_multi_batch", NULL, placeholder = "Batches (.rds / .csv / .txt)", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          checkboxInput("chk_append_multi", "Add to current benchmark results", value = FALSE),
          actionButton("btn_run_multi_eval", "Calculate Multiomics Metrics", class = "btn btn-primary w-100", icon = icon("play"))
        ),
        
        # Mode 4: Saved Results
        conditionalPanel(
          condition = "input.opt_data_mode == 'upload_bench'",
          p("Load a previously saved benchmark result file (.rds or .csv) from scSimEval.", style = "font-size: 0.82rem; color: #475569;"),
          fileInput("file_bench_upload", "Select Saved File (.rds / .csv):", accept = c(".rds", ".csv")),
          actionButton("btn_load_uploaded_bench", "Load Saved Results", class = "btn btn-secondary w-100", icon = icon("folder-open"))
        )
      ),
      
      uiOutput("ui_datahub_main_panel")
    )
  ),
  # ============================================================================
  # TAB 3: COMPARATIVE BUBBLE MATRIX (FLAGSHIP FIGURE)
  # ============================================================================
  nav_panel(
    "Comparative Bubble Matrix",
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        title = div(icon("sliders", class = "me-2 text-primary"), tags$strong("Display & Filters")),
        
        # 1. Filter Mode
        radioButtons(
          "rad_bubble_filter_mode", tags$b("Filter Mode:"),
          choices = c(
            "All 62 Measures"  = "all",
            "Custom Selection" = "custom"
          ),
          selected = "all"
        ),
        
        # 2. Custom Selection Panel
        conditionalPanel(
          condition = "input.rad_bubble_filter_mode == 'custom'",
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; padding: 12px; margin-top: 4px; margin-bottom: 12px;",
            uiOutput("ui_custom_category_picker"),
            div(
              class = "d-flex justify-content-between align-items-center mt-2 mb-1",
              tags$b("2. Select Metrics:", style = "font-size: 0.88rem; color: #1E293B;"),
              div(
                actionLink("btn_custom_metrics_all", "All", style = "font-size: 0.78rem; margin-right: 6px; font-weight: 600; text-decoration: none;"),
                actionLink("btn_custom_metrics_none", "Clear", style = "font-size: 0.78rem; color: #dc3545; text-decoration: none;")
              )
            ),
            div(
              style = "max-height: 190px; overflow-y: auto; background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 6px 10px;",
              uiOutput("ui_custom_metrics_checklist")
            ),
            div(
              style = "margin-top: 5px; font-size: 0.76rem; color: #64748B; text-align: right;",
              textOutput("txt_custom_metrics_count", inline = TRUE)
            )
          )
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # 3. Simulator Methods Selection
        uiOutput("ui_bubble_method_picker"),
        
        hr(style = "margin: 10px 0;"),
        
        # 4. Collapsible Dimension Sliders
        accordion(
          open = FALSE,
          accordion_panel(
            title = span(icon("up-right-and-down-left-from-center"), " Display & Export Sizing"),
            value = "panel_dimensions",
            sliderInput("sld_bubble_width", "Matrix Width (px):", min = 800, max = 3200, value = 2200, step = 50),
            sliderInput("sld_bubble_height", "Matrix Height (px):", min = 350, max = 1200, value = 680, step = 20),
            p("Adjust plot dimensions for screen display or file download.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
          )
        ),
        
        hr(style = "margin: 12px 0;"),
        
        # 5. Clean download section
        div(
          class = "d-grid gap-2",
          downloadButton("download_bubble_jpeg", "Download JPEG (600 DPI)", class = "btn btn-primary btn-sm"),
          downloadButton("download_bubble_pdf", "Download Vector PDF", class = "btn btn-outline-secondary btn-sm")
        )
      ),
      
      navset_pill(
        id = "bubble_matrix_subtabs",
        selected = "Comparative Bubble Matrix Figure",
        
        # Subtab 1: Figure
        nav_panel(
          "Comparative Bubble Matrix Figure",
          card(
            card_header(
              div(
                class = "d-flex justify-content-between align-items-center",
                span(icon("chart-dot", class = "me-2 text-primary"), tags$b("Simulator Comparison Bubble Matrix")),
                uiOutput("ui_bubble_active_badge")
              )
            ),
            card_body(
              p("Each column represents a simulation method; each row represents an evaluation metric. Larger bubbles show higher scores (closer match to real data). Colors indicate metric categories.", style = "font-size: 0.84rem; color: #64748B; margin-bottom: 12px;"),
              div(
                class = "bubble-scroll-container",
                uiOutput("ui_bubble_plot_render")
              )
            )
          )
        ),
        
        # Subtab 2: Fidelity Leaderboard
        nav_panel(
          "Fidelity Leaderboard",
          card(
            card_header(
              class = "d-flex justify-content-between align-items-center py-2",
              tags$span(icon("ranking-star", class = "me-2 text-primary"), tags$b("Simulator Rankings")),
              downloadButton("download_leaderboard_csv", "Download Leaderboard (CSV)", class = "btn btn-sm btn-outline-primary")
            ),
            card_body(
              p(
                "Overall simulator rankings averaged across all evaluation metrics. Metrics where lower values are better (such as error, runtime, and memory) are inverted, and all scores are scaled between 0.00 and 1.00:",
                style = "font-size: 0.84rem; color: #64748B; margin-bottom: 12px;"
              ),
              DTOutput("table_leaderboard_dt")
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 4: DIAGNOSTIC VISUALIZATIONS (8 SUB-PANELS)
  # ============================================================================
  nav_panel(
    "Visualizations",
    navset_pill(
      id = "viz_subtabs",
      
      # Sub-panel 1: Evaluation Summary
      nav_panel(
        "1. Evaluation Summary",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("chart-pie", class = "me-2 text-primary"), tags$b("Category Evaluation Summary")),
              span(class = "badge bg-light text-secondary border", "Radar & Bar Chart")
            )
          ),
          card_body(
            fluidRow(
              column(4, checkboxInput("chk_sum_labels", "Show Score Labels", value = TRUE)),
              column(4, checkboxInput("chk_sum_norm", "Normalize Scores [0, 1]", value = TRUE)),
              column(4,
                     downloadButton("download_sum_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_sum_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_eval_summary", height = "540px")
          )
        )
      ),
      
      # Sub-panel 2: Distribution QC
      nav_panel(
        "2. Distribution QC",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("chart-line", class = "me-2 text-primary"), tags$b("Real vs Simulated Data Distributions")),
              span(class = "badge bg-light text-secondary border", "Density Curves")
            )
          ),
          card_body(
            fluidRow(
              column(4, selectInput("sel_dist_layout", "QC Layout:", choices = c("Comprehensive" = "comprehensive", "Density Curves Only" = "density"), selected = "comprehensive")),
              column(4, p("Compares gene expression, library sizes, and zeros between real reference data and simulated cells.", style = "font-size: 0.82rem; color: #64748B; margin-bottom: 0;")),
              column(4,
                     downloadButton("download_dist_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_dist_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            uiOutput("ui_dist_qc_plot")
          )
        )
      ),
      
      # Sub-panel 3: Scalability Benchmark
      nav_panel(
        "3. Scalability Benchmark",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("microchip", class = "me-2 text-primary"), tags$b("Runtime and Memory Usage")),
              span(class = "badge bg-light text-secondary border", "Computational Cost")
            )
          ),
          card_body(
            fluidRow(
              column(6,
                     selectInput(
                       "sel_scale_type", "Scalability View:",
                       choices = c(
                         "All 4 Panels Combined" = "composite",
                         "Runtime (Seconds)" = "runtime",
                         "Peak RAM Memory (MB)" = "memory",
                         "Runtime vs Memory Trade-Off" = "tradeoff",
                         "Overall Resource Cost" = "cost"
                       ),
                       selected = "composite"
                     )
              ),
              column(6,
                     downloadButton("download_scale_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_scale_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_scale_bench", height = "620px")
          )
        )
      ),
      
      # Sub-panel 4: Metric Plots
      nav_panel(
        "4. Metric Plots",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("chart-column", class = "me-2 text-primary"), tags$b("Metric Plots by Category")),
              span(class = "badge bg-light text-secondary border", "Boxplots & Bars")
            )
          ),
          card_body(
            fluidRow(
              column(12,
                     radioButtons(
                       "opt_box_view_mode", "View Mode:",
                       choices = c("View Individual Metric" = "individual", "View by Category Group" = "category"),
                       selected = "individual", inline = TRUE
                     )
              )
            ),
            fluidRow(
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'individual'",
                column(3,
                       selectInput(
                         "sel_box_cat_first", "1. Choose Category First:",
                         choices = c(
                           "(I) Distributional Properties",
                           "(II) Correlations & Zero-Inflation",
                           "(III) Cellular Structure & Concordance",
                           "(IV) Batch Effects & Confounder Mixing",
                           "(V) Biological Signal & Downstream Fidelity",
                           "(VI) Trajectory & Lineage Dynamics",
                           "(VII) Cross-Modal Coupling & Modularity",
                           "(VIII) Computational Scalability"
                         ),
                         selected = "(I) Distributional Properties"
                       )
                ),
                column(3, uiOutput("ui_box_metric_picker")),
                column(3,
                       selectInput(
                         "sel_box_indiv_score", "3. Score Type:",
                         choices = c("Normalized Score [0, 1]" = "normalized", "Raw Metric Value" = "raw"),
                         selected = "normalized"
                       )
                ),
                column(3,
                       downloadButton("download_box_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                       downloadButton("download_box_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
                )
              ),
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'category'",
                column(4,
                       selectInput(
                         "sel_box_cat_group", "Choose Category:",
                         choices = c(
                           "All 8 Categories" = "all",
                           "(I) Distributional Properties",
                           "(II) Correlations & Zero-Inflation",
                           "(III) Cellular Structure & Concordance",
                           "(IV) Batch Effects & Confounder Mixing",
                           "(V) Biological Signal & Downstream Fidelity",
                           "(VI) Trajectory & Lineage Dynamics",
                           "(VII) Cross-Modal Coupling & Modularity",
                           "(VIII) Computational Scalability"
                         ),
                         selected = "all"
                       )
                ),
                column(4, selectInput("sel_box_score_type", "Score Type:", choices = c("Normalized Score [0, 1]" = "normalized", "Raw Value" = "raw"), selected = "normalized")),
                column(4,
                       downloadButton("download_box_cat_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                       downloadButton("download_box_cat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
                )
              )
            ),
            hr(),
            uiOutput("ui_plot_metric_boxes")
          )
        )
      ),
      
      # Sub-panel 5: Metric Heatmap
      nav_panel(
        "5. Metric Heatmap",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("table-cells", class = "me-2 text-primary"), tags$b("Performance Heatmap")),
              span(class = "badge bg-light text-secondary border", "Overview Table")
            )
          ),
          card_body(
            fluidRow(
              column(3,
                     selectInput(
                       "sel_heat_cat", "Category Filter:",
                       choices = c(
                         "All 8 Categories Combined" = "all",
                         "(I) Distributional Properties",
                         "(II) Correlations & Zero-Inflation",
                         "(III) Cellular Structure & Concordance",
                         "(IV) Batch Effects & Confounder Mixing",
                         "(V) Biological Signal & Downstream Fidelity",
                         "(VI) Trajectory & Lineage Dynamics",
                         "(VII) Cross-Modal Coupling & Modularity",
                         "(VIII) Computational Scalability"
                       ),
                       selected = "all"
                     )
              ),
              column(3,
                     sliderInput("sld_heat_height", "Heatmap Height (px):", min = 400, max = 1500, value = 1100, step = 50)
              ),
              column(3,
                     sliderInput("sld_heat_width", "Heatmap Width (px):", min = 600, max = 1300, value = 880, step = 20)
              ),
              column(3,
                     downloadButton("download_heat_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_heat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_heat")
            )
          )
        )
      ),
      
      # Sub-panel 6: PCA Ordination
      nav_panel(
        "6. PCA Ordination",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("compass", class = "me-2 text-primary"), tags$b("PCA Plot of Simulators")),
              span(class = "badge bg-light text-secondary border", "Metric Space")
            )
          ),
          card_body(
            fluidRow(
              column(3,
                     selectInput(
                       "sel_pca_cat", "Category Filter:",
                       choices = c(
                         "All Categories Combined" = "all",
                         "(I) Distributional Properties" = "(I) Distributional Properties",
                         "(II) Correlations & Zero-Inflation" = "(II) Correlations & Zero-Inflation",
                         "(III) Cellular Structure & Concordance" = "(III) Cellular Structure & Concordance",
                         "(IV) Batch Effects & Confounder Mixing" = "(IV) Batch Effects & Confounder Mixing",
                         "(V) Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
                         "(VII) Cross-Modal Coupling & Modularity" = "(VII) Cross-Modal Coupling & Modularity"
                       ),
                       selected = "all"
                     )
              ),
              column(3, selectInput("sel_pca_panel", "PCA Panel:", choices = c("All Panels (Biplot + Loadings + Scores + Variance)" = "all", "Biplot Only" = "biplot", "Metric Loadings Only" = "loadings", "Simulator Scores Only" = "scores", "Variance Explained Only" = "scree"), selected = "all")),
              column(3, numericInput("num_pca_top_metrics", "Top Metrics to Label:", value = 8, min = 3, max = 25, step = 1)),
              column(3,
                     downloadButton("download_pca_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_pca_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            fluidRow(
              column(6,
                     sliderInput("sld_pca_height", "PCA Height (px):", min = 500, max = 1600, value = 1100, step = 50)
              ),
              column(6,
                     sliderInput("sld_pca_width", "PCA Width (px):", min = 650, max = 1400, value = 950, step = 25)
              )
            ),
            div(
              style = "margin-bottom: 8px; font-size: 0.80rem; color: #64748B; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px 12px;",
              icon("circle-info", class = "text-info me-1"),
              tags$b("Tip: "),
              "Arrows show the metrics that most strongly separate high-performing simulators from lower-performing ones."
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_pca")
            )
          )
        )
      ),
      
      # Sub-panel 7: MDS Metric Space
      nav_panel(
        "7. MDS Metric Space",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("circle-nodes", class = "me-2 text-primary"), tags$b("MDS Distance Plot")),
              span(class = "badge bg-light text-secondary border", "Simulator Similarity")
            )
          ),
          card_body(
            fluidRow(
              column(3,
                     selectInput(
                       "sel_mds_cat", "Category Filter:",
                       choices = c(
                         "All Categories Combined" = "all",
                         "(I) Distributional Properties" = "(I) Distributional Properties",
                         "(II) Correlations & Zero-Inflation" = "(II) Correlations & Zero-Inflation",
                         "(III) Cellular Structure & Concordance" = "(III) Cellular Structure & Concordance",
                         "(IV) Batch Effects & Confounder Mixing" = "(IV) Batch Effects & Confounder Mixing",
                         "(V) Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
                         "(VII) Cross-Modal Coupling & Modularity" = "(VII) Cross-Modal Coupling & Modularity"
                       ),
                       selected = "all"
                     )
              ),
              column(3, selectInput("sel_mds_by", "Compare By:", choices = c("By Simulators" = "simulators", "By Metric Summaries" = "summaries"), selected = "simulators")),
              column(3,
                     sliderInput("sld_mds_height", "MDS Height (px):", min = 400, max = 1300, value = 620, step = 20)
              ),
              column(3,
                     sliderInput("sld_mds_width", "MDS Width (px):", min = 600, max = 1400, value = 900, step = 20)
              )
            ),
            fluidRow(
              column(12,
                     downloadButton("download_mds_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_mds_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_mds")
            )
          )
        )
      ),
      
      # Sub-panel 8: Cell Embeddings (t-SNE & UMAP)
      nav_panel(
        "8. Cell Embeddings (t-SNE & UMAP)",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("project-diagram", class = "me-2 text-primary"), tags$b("Cell Projections (UMAP, t-SNE, PCA)")),
              span(class = "badge bg-light text-secondary border", "Cell Clusters")
            )
          ),
          card_body(
            fluidRow(
              column(3,
                     radioButtons(
                       "sel_emb_reduction", "1. Method:",
                       choices = c("UMAP" = "umap", "t-SNE" = "tsne", "PCA" = "pca"),
                       selected = "umap", inline = TRUE
                     ),
                     selectInput(
                       "sel_emb_layout", "2. Layout:",
                       choices = c(
                         "All Simulators Grid" = "facet",
                         "Side-by-Side (Reference vs Single Simulator)" = "side_by_side"
                       ),
                       selected = "facet"
                     ),
                     conditionalPanel(
                       condition = "input.sel_emb_layout == 'side_by_side'",
                       selectInput("sel_emb_single_sim", "Select Simulator to Compare:", choices = NULL)
                     )
              ),
              column(3,
                     uiOutput("ui_emb_sim_picker"),
                     selectInput(
                       "sel_emb_color", "Color Cells By:",
                       choices = c(
                         "Cell Type / Group" = "cell_type",
                         "Cluster" = "cluster",
                         "Library Size (Depth)" = "library_size",
                         "Batch" = "batch",
                         "Dataset Source" = "dataset"
                       ),
                       selected = "cell_type"
                     )
              ),
              column(3,
                     sliderInput("sld_emb_pcs", "Number of PCs:", min = 5, max = 50, value = 20, step = 5),
                     conditionalPanel(
                       condition = "input.sel_emb_reduction == 'tsne'",
                       sliderInput("sld_emb_perp", "t-SNE Perplexity:", min = 5, max = 50, value = 15, step = 5)
                     ),
                     conditionalPanel(
                       condition = "input.sel_emb_reduction == 'umap'",
                       sliderInput("sld_emb_neighbors", "UMAP Neighbors:", min = 5, max = 50, value = 15, step = 5)
                     )
              ),
              column(3,
                     fluidRow(
                       column(6, sliderInput("sld_emb_pt_size", "Point Size:", min = 0.2, max = 3.0, value = 1.0, step = 0.1)),
                       column(6, sliderInput("sld_emb_alpha", "Alpha:", min = 0.2, max = 1.0, value = 0.8, step = 0.05))
                     ),
                     div(class = "mt-2",
                         downloadButton("download_emb_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                         downloadButton("download_emb_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary mb-1")
                     )
              )
            ),
            hr(),
            div(
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_cell_embeddings")
            ),
            hr(style = "margin: 20px 0; border-color: #CBD5E1;"),
            
            # Quantitative Metrics Box
            div(
              class = "p-3", style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px;",
              div(
                class = "d-flex justify-content-between align-items-center mb-2",
                h6(icon("chart-line", class = "text-primary me-2"), tags$b("Cluster Separation and Quality Metrics"), style = "color: #0F172A; margin: 0; font-size: 0.92rem;"),
                span(class = "badge bg-light text-secondary border", "Silhouette • ARI • Library Depth")
              ),
              p("Shows how well cell clusters and library sizes are preserved between real and simulated data:", style = "font-size: 0.82rem; color: #64748B; margin-bottom: 10px;"),
              DTOutput("table_emb_quality_metrics")
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 5: DOWNLOAD RESULTS
  # ============================================================================
  nav_panel(
    "Download Results",
    fluidRow(
      column(
        6,
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("file-zipper", class = "me-2 text-success"), tags$b("Complete Results Archive (.zip)")),
              span(class = "badge bg-success-subtle text-success border border-success-subtle", ".zip")
            )
          ),
          card_body(
            p("Download all benchmark tables, figures, and reports in a single zip file:", style = "font-size: 0.86rem; color: #475569; margin-bottom: 14px;"),
            div(
              class = "row g-2 mb-3",
              div(class = "col-6",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 12px; height: 100%;",
                  div(icon("file-excel", class = "text-success me-1"), tags$b("Excel Workbook", style = "font-size: 0.82rem; color: #1E293B;")),
                  p("Full results with rankings, properties & category breakdown.", style = "font-size: 0.76rem; color: #64748B; margin: 2px 0 0 0;")
                )
              ),
              div(class = "col-6",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 12px; height: 100%;",
                  div(icon("file-pdf", class = "text-danger me-1"), tags$b("Multi-Page PDF", style = "font-size: 0.82rem; color: #1E293B;")),
                  p("Compiled report containing all 11 evaluation figures.", style = "font-size: 0.76rem; color: #64748B; margin: 2px 0 0 0;")
                )
              ),
              div(class = "col-6",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 12px; height: 100%;",
                  div(icon("table", class = "text-primary me-1"), tags$b("Master Tables", style = "font-size: 0.82rem; color: #1E293B;")),
                  p("All 62 evaluation metrics in .csv and .tsv formats.", style = "font-size: 0.76rem; color: #64748B; margin: 2px 0 0 0;")
                )
              ),
              div(class = "col-6",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 12px; height: 100%;",
                  div(icon("code", class = "text-secondary me-1"), tags$b("R Object (.rds)", style = "font-size: 0.82rem; color: #1E293B;")),
                  p("R data file for custom plotting and downstream analysis.", style = "font-size: 0.76rem; color: #64748B; margin: 2px 0 0 0;")
                )
              )
            ),
            downloadButton("download_complete_zip", "Download Complete Results Archive (.zip)", class = "btn btn-success w-100 py-2", icon = icon("file-zipper"))
          )
        )
      ),
      column(
        6,
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("table-cells", class = "me-2 text-primary"), tags$b("Data Tables & Files")),
              span(class = "badge bg-light text-secondary border", "Tables & RDS")
            )
          ),
          card_body(
            p("Export benchmark score tables for Excel, Google Sheets, or R:", style = "font-size: 0.86rem; color: #475569; margin-bottom: 14px;"),
            downloadButton("download_excel", "Download Excel Workbook (.xlsx)", class = "btn btn-primary w-100 mb-2", icon = icon("file-excel")),
            downloadButton("download_csv", "Download CSV Table (.csv)", class = "btn btn-outline-primary w-100 mb-2", icon = icon("file-csv")),
            downloadButton("download_txt", "Download Text Table (.txt)", class = "btn btn-outline-secondary w-100 mb-2", icon = icon("file-lines")),
            downloadButton("download_rds", "Download R Data File (.rds)", class = "btn btn-outline-secondary w-100", icon = icon("code"))
          )
        )
      )
    ),
    fluidRow(
      style = "margin-top: 18px;",
      column(
        6,
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("file-pdf", class = "me-2 text-danger"), tags$b("Full PDF Report")),
              span(class = "badge bg-danger-subtle text-danger border border-danger-subtle", "11 Figures")
            )
          ),
          card_body(
            p("A multi-page PDF document with all 11 benchmark figures:", style = "font-size: 0.86rem; color: #475569; margin-bottom: 12px;"),
            div(
              style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 14px; margin-bottom: 14px; font-size: 0.78rem; color: #475569;",
              div(class = "d-flex flex-wrap gap-2",
                span(class = "badge bg-light text-secondary border", "Bubble Matrix"),
                span(class = "badge bg-light text-secondary border", "Evaluation Summary"),
                span(class = "badge bg-light text-secondary border", "Scalability Benchmark"),
                span(class = "badge bg-light text-secondary border", "Metric Distributions"),
                span(class = "badge bg-light text-secondary border", "Performance Heatmap"),
                span(class = "badge bg-light text-secondary border", "PCA Plot"),
                span(class = "badge bg-light text-secondary border", "MDS Plot"),
                span(class = "badge bg-light text-secondary border", "Distribution QC"),
                span(class = "badge bg-light text-secondary border", "UMAP / t-SNE / PCA Embeddings")
              )
            ),
            downloadButton("download_all_plots_pdf", "Download Full PDF Report (.pdf)", class = "btn btn-danger w-100 py-2", icon = icon("file-pdf"))
          )
        )
      ),
      column(
        6,
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("image", class = "me-2 text-info"), tags$b("Export Single Figure")),
              span(class = "badge bg-light text-secondary border", "PDF & 600 DPI JPEG")
            )
          ),
          card_body(
            p("Choose any figure to download at publication quality:", style = "font-size: 0.86rem; color: #475569; margin-bottom: 10px;"),
            selectInput(
              "sel_export_figure_type",
              "Select Figure:",
              choices = list(
                "Benchmark Overviews" = c(
                  "Comparative Bubble Matrix" = "bubble",
                  "Overall Evaluation Summary" = "summary",
                  "Scalability Benchmark" = "scalability",
                  "Metric Plots by Category (Boxplots / Barplots)" = "boxplots",
                  "Performance Heatmap" = "heatmap",
                  "Simulator PCA Ordination (Metrics)" = "pca_metric",
                  "Simulator MDS Ordination" = "mds_metric",
                  "Distribution QC Curves" = "dist_qc"
                ),
                "Cell Embeddings (2D Projections)" = c(
                  "UMAP Embeddings: All Simulators Grid" = "emb_umap",
                  "t-SNE Embeddings: All Simulators Grid" = "emb_tsne",
                  "PCA Embeddings: All Simulators Grid" = "emb_pca",
                  "UMAP Embeddings: Pairwise 1-to-1 Comparison" = "emb_compare_umap",
                  "t-SNE Embeddings: Pairwise 1-to-1 Comparison" = "emb_compare_tsne",
                  "PCA Embeddings: Pairwise 1-to-1 Comparison" = "emb_compare_pca"
                )
              ),
              selected = "emb_umap"
            ),
            conditionalPanel(
              condition = "input.sel_export_figure_type.indexOf('compare') !== -1",
              selectInput("sel_export_compare_sim", "Select Simulator for 1-to-1 Comparison:", choices = NULL)
            ),
            div(
              class = "d-flex gap-2 mt-3",
              downloadButton("download_selected_plot_pdf", "Vector PDF", class = "btn btn-outline-primary flex-fill", icon = icon("file-pdf")),
              downloadButton("download_selected_plot_jpeg", "JPEG (600 DPI)", class = "btn btn-primary flex-fill", icon = icon("file-image"))
            )
          )
        )
      )
    ),
    fluidRow(
      style = "margin-top: 18px;",
      column(
        12,
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              span(icon("table-list", class = "me-2 text-primary"), tags$b("Benchmark Results Table")),
              span(class = "badge bg-light text-secondary border", "Search & Filter")
            )
          ),
          card_body(
            p("Search, sort, and filter benchmark scores across simulators and metric categories:", style = "font-size: 0.84rem; color: #64748B; margin-bottom: 12px;"),
            DTOutput("table_master_export")
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 6: HELP & GETTING STARTED
  # ============================================================================
  nav_panel(
    "Help & Getting Started",
    fluidRow(
      column(
        12,
        card(
          card_body(
            # Introduction Hero Callout
            div(
              style = "background: #F8FAFC; border-left: 4px solid #1B4F72; border-radius: 6px; padding: 18px 22px; margin-bottom: 24px;",
              h4("Unified Evaluation & Benchmarking for Single-Cell Simulations", style = "font-weight: 800; color: #1B4F72; margin-top: 0; font-size: 1.18rem;"),
              p(
                strong("scSimEval"), " is a scientific framework to test how well synthetic single-cell datasets reproduce real biological experiments. ",
                "It provides a ground-truth-free evaluation pipeline with ", strong("62 evaluation metrics organized into 8 biological and computational categories"), 
                " across single-cell RNA-seq, ATAC-seq, and paired multiomics."
              ),
              p(
                "This web application lets you load example benchmarks or upload your own datasets, calculate scores, compare simulators, and download publication-ready figures.",
                style = "margin-bottom: 0; color: #475569; font-size: 0.88rem;"
              ),
              div(
                style = "margin-top: 14px; padding-top: 12px; border-top: 1px solid #CBD5E1; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px;",
                span(
                  tags$i(class = "fa fa-book", style = "margin-right: 6px; color: #1B4F72;"),
                  strong("Online Documentation & Vignettes: "),
                  tags$a(
                    href = "https://kabilanbio.github.io/scSimEval",
                    target = "_blank",
                    rel = "noopener noreferrer",
                    "https://kabilanbio.github.io/scSimEval",
                    style = "color: #0284C7; text-decoration: underline; font-weight: 700; font-size: 0.90rem; margin-left: 4px;"
                  )
                ),
                tags$a(
                  href = "https://kabilanbio.github.io/scSimEval",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  class = "btn btn-sm btn-primary",
                  style = "font-weight: 600; border-radius: 6px;",
                  tags$i(class = "fa fa-external-link-alt", style = "margin-right: 5px;"),
                  "Open Online Documentation"
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 1: Workflow Architecture
            # ------------------------------------------------------------------
            h4("1. Workflow Architecture & Pipeline", style = "font-weight: 700; color: #1B4F72; margin-top: 10px; font-size: 1.08rem;"),
            p("The benchmarking process follows a standardized four-step pipeline using real reference and simulated count matrices:", style = "color: #475569; font-size: 0.88rem;"),
            
            div(
              style = "text-align: center; margin: 20px auto 16px auto; max-width: 980px;",
              tags$img(
                src = "scfigures/workflow_diagram.png",
                alt = "scSimEval Benchmarking Workflow Architecture",
                style = "width: 100%; max-width: 960px; height: auto; display: block; margin: 0 auto; border-radius: 8px; border: 1px solid #E2E8F0;"
              ),
              p(
                tags$b("Figure 1 | The scSimEval Benchmarking Workflow. "),
                "Step 1: Input real reference and simulated count matrices with metadata. Step 2: Calculate 62 evaluation metrics across 8 categories. Step 3: Standardize and consolidate all scores into a benchmark table. Step 4: Explore simulator rankings, bubble matrices, and high-resolution figures.",
                style = "font-size: 0.84rem; color: #64748B; margin-top: 10px; max-width: 960px; margin-left: auto; margin-right: auto;"
              )
            ),
            
            tags$div(
              class = "row g-3 my-2",
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 1: Input Real & Simulated Data"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("Provide real experimental reference matrices and simulated count matrices (genes x cells for scRNA-seq, peaks x cells for scATAC-seq), along with cell type labels, batch labels, and computational resource records (runtime and peak RAM).", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 2: Calculate Evaluation Metrics"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("scSimEval calculates 62 quantitative measures across 8 biological and computational categories directly comparing simulated data to real experimental data.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 3: Run Consolidated Benchmark Pipeline"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("The benchmark engine standardizes and consolidates evaluation scores across all evaluated simulators into a tidy summary table, mapping each score to its category.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 4: Review Scores, Rankings & Figures"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("Inspect the bubble matrix, simulator rankings, distribution curves, PCA biplots, and MDS spaces. Export publication figures at 600 DPI alongside Excel/CSV data archives.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 2: The Eight Evaluation Categories
            # ------------------------------------------------------------------
            h4("2. The Eight Evaluation Categories", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("The 62 evaluation criteria in scSimEval are structured across eight foundational categories covering statistical, cellular, molecular, and computational dimensions:", style = "color: #475569; font-size: 0.88rem;"),
            
            tags$div(
              class = "table-responsive",
              tags$table(
                class = "table table-hover table-bordered table-sm",
                style = "font-size: 0.84rem;",
                tags$thead(
                  class = "table-light",
                  tags$tr(
                    tags$th(style = "width: 5%; text-align: center;", "#"),
                    tags$th(style = "width: 25%;", "Category"),
                    tags$th(style = "width: 15%; text-align: center;", "Focus"),
                    tags$th(style = "width: 55%;", "Representative Metrics & Scope")
                  )
                ),
                tags$tbody(
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "I"),
                    tags$td(tags$b("Distributional Properties")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-primary border", "Statistical")),
                    tags$td("Mean expression, variance, dispersion (CV, Fano factor), library size, detected feature rate, and distribution distance (Wasserstein, KS statistic).")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "II"),
                    tags$td(tags$b("Correlations & Zero-Inflation")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-info border", "Co-expression")),
                    tags$td("Gene-gene Spearman/Pearson correlation networks, dropout rate curves, zero proportion concordance, and mean-variance trend fidelity.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "III"),
                    tags$td(tags$b("Cellular Structure & Concordance")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-success border", "Clustering")),
                    tags$td("Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), Silhouette width, cluster purity, and cell-cell distance matrix concordance.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "IV"),
                    tags$td(tags$b("Batch Effects & Confounder Mixing")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-warning border", "Technical")),
                    tags$td("kBET (batch test), LISI (Local Inverse Simpson Index), Average Silhouette Width for Batch (ASW-batch), and principal component batch correlation.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "V"),
                    tags$td(tags$b("Biological Signal & Downstream Fidelity")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-danger border", "Phenotype")),
                    tags$td("Differentially Expressed Gene (DEG) concordance (Jaccard index, sensitivity, specificity, logFC correlation), marker gene ranking, and pathway enrichment overlap.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VI"),
                    tags$td(tags$b("Trajectory & Lineage Dynamics")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-secondary border", "Pseudotime")),
                    tags$td("Pseudotime correlation with reference, trajectory topology preservation (graph alignment), and gene expression patterns along lineages.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VII"),
                    tags$td(tags$b("Cross-Modal Coupling & Modularity")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-indigo border", "Multiomics")),
                    tags$td("Peak-to-gene linkage concordance, cross-modal cell pairing accuracy (FOSCTTM, Match@1), TF regulatory consistency, and multi-layer cluster concordance.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VIII"),
                    tags$td(tags$b("Computational Scalability")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-dark border", "Efficiency")),
                    tags$td("Wall-clock runtime (seconds), peak memory usage (MB RAM), cell-scaling rate, and feature-scaling efficiency.")
                  )
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 3: How Scores are Calculated and Normalized (RESTORED FIGURE 2)
            # ------------------------------------------------------------------
            h4("3. How Scores are Calculated and Normalized", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("To fairly compare different metrics with different units and scales, scSimEval standardizes all scores using a two-step normalization pipeline:", style = "color: #475569; font-size: 0.88rem;"),
            
            div(
              style = "text-align: center; margin: 20px auto 16px auto; max-width: 950px;",
              tags$img(
                src = "scfigures/score_normalization_workflow.png",
                alt = "Two-Step Score Normalization & Visual Mapping Pipeline",
                style = "width: 100%; max-width: 920px; height: auto; display: block; margin: 0 auto; border-radius: 8px; border: 1px solid #E2E8F0;"
              ),
              p(
                tags$b("Figure 2 | Score Normalization Pipeline. "),
                "Metrics where lower values indicate better results (such as error, runtime, and memory) are inverted so that higher is always better. All values are then scaled between 0.00 and 1.00 so simulators can be ranked and compared fairly in the bubble matrix and summary charts.",
                style = "font-size: 0.84rem; color: #64748B; margin-top: 10px; max-width: 920px; margin-left: auto; margin-right: auto;"
              )
            ),
            
            tags$div(
              class = "row g-3 my-2",
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("1. Direction Inversion"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("Metrics where lower values are better (such as error, distance, runtime, or memory) are inverted:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Inverted = Max - Value"), style = "text-align: center; font-weight: bold; font-size: 0.82rem;"),
                  p("This ensures higher numbers always indicate better simulation quality.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("2. Scale Between 0 and 1"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("All metrics are rescaled across simulation methods onto a 0.0 to 1.0 scale:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Score = (Value - Min) / (Max - Min)"), style = "text-align: center; font-weight: bold; font-size: 0.82rem;"),
                  p("0.0 represents the lowest performing simulator, and 1.0 represents the best agreement with real data.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("3. Average Overall Score"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("Overall simulator fidelity is calculated by taking the average of all standardized scores:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Fidelity = mean(Scores)"), style = "text-align: center; font-weight: bold; font-size: 0.82rem;"),
                  p("Shown as a score from 0 to 1 and as a percentage (100% means top performance across all metrics).", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("4. How Plots Display Scores"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  tags$ul(
                    style = "font-size: 0.78rem; color: #475569; padding-left: 14px; margin-bottom: 0;",
                    tags$li(tags$b("Bubble Size: "), "Larger bubbles indicate higher scores."),
                    tags$li(tags$b("Top Scores: "), "Scores >= 0.96 are highlighted with square markers."),
                    tags$li(tags$b("Colors: "), "Each metric category has its own distinct color."),
                    tags$li(tags$b("Ranking: "), "Simulators are ordered by their average score.")
                  )
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 4: Guide to the App Tabs
            # ------------------------------------------------------------------
            h4("4. Guide to the App Tabs", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("Navigate through the application tabs to perform simulation evaluations:", style = "color: #475569; font-size: 0.88rem;"),
            
            tags$div(
              class = "accordion my-3", id = "accordionGuide",
              
              # Module A: Home
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingOne",
                  tags$button(
                    class = "accordion-button", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseOne", "aria-expanded" = "true", "aria-controls" = "collapseOne",
                    tags$b("Tab 1: Home — Framework Overview & Quick Launch")
                  )
                ),
                tags$div(
                  id = "collapseOne", class = "accordion-collapse collapse show", "aria-labelledby" = "headingOne",
                  tags$div(
                    class = "accordion-body",
                    p("Presents the executive summary of scSimEval, metric counts across all 8 categories, and direct buttons to jump to each module.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              ),
              
              # Module B: Data Hub
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingTwo",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseTwo", "aria-expanded" = "false", "aria-controls" = "collapseTwo",
                    tags$b("Tab 2: Data Hub — Load Examples or Upload Data")
                  )
                ),
                tags$div(
                  id = "collapseTwo", class = "accordion-collapse collapse", "aria-labelledby" = "headingTwo",
                  tags$div(
                    class = "accordion-body",
                    tags$ul(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 16px; margin: 0;",
                      tags$li(tags$b("Example Benchmark (6 Simulators): "), "Click to instantly load precomputed evaluations across Splatter, scDesign3, SCRIP, SymSim, dyngen, and simATAC."),
                      tags$li(tags$b("Single-Cell (scRNA / scATAC): "), "Upload real reference matrices and simulated datasets (.rds, .csv, .tsv, .txt) to automatically calculate evaluation metrics."),
                      tags$li(tags$b("Multiomics (scRNA + scATAC): "), "Evaluate paired or unpaired datasets assessing cross-modal coupling alongside single-cell properties."),
                      tags$li(tags$b("Upload Saved Results (.rds / .csv): "), "Restore saved evaluation files from previous runs.")
                    )
                  )
                )
              ),
              
              # Module C: Bubble Matrix
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingThree",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseThree", "aria-expanded" = "false", "aria-controls" = "collapseThree",
                    tags$b("Tab 3: Comparative Bubble Matrix & Leaderboard")
                  )
                ),
                tags$div(
                  id = "collapseThree", class = "accordion-collapse collapse", "aria-labelledby" = "headingThree",
                  tags$div(
                    class = "accordion-body",
                    p("Interactive visualization where bubble sizes show standardized scores and colors show categories. Switch to the Fidelity Leaderboard subtab to view composite rankings.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              ),
              
              # Module D: Visualizations
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingFour",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseFour", "aria-expanded" = "false", "aria-controls" = "collapseFour",
                    tags$b("Tab 4: Visualizations — 8 Diagnostic Panels")
                  )
                ),
                tags$div(
                  id = "collapseFour", class = "accordion-collapse collapse", "aria-labelledby" = "headingFour",
                  tags$div(
                    class = "accordion-body",
                    p("Explore eight dedicated diagnostic subtabs: Category Evaluation Summary, Distribution Quality QC curves, Computational Scalability, Metric Boxplots & Barplots, Performance Heatmap, PCA plots, MDS plots, and Cell Projections (UMAP, t-SNE, PCA).", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              ),
              
              # Module E: Download Results
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingFive",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseFive", "aria-expanded" = "false", "aria-controls" = "collapseFive",
                    tags$b("Tab 5: Download Results — Reports & Data Files")
                  )
                ),
                tags$div(
                  id = "collapseFive", class = "accordion-collapse collapse", "aria-labelledby" = "headingFive",
                  tags$div(
                    class = "accordion-body",
                    p("Export all benchmarking results as a complete .zip archive, Excel workbook (.xlsx), CSV/TXT tables, serialized RDS object, full multi-page PDF report, or individual publication figures.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 5: Online Documentation Box
            # ------------------------------------------------------------------
            div(
              style = "margin-top: 24px; padding: 18px 20px; background: #F0F9FF; border: 1px solid #BAE6FD; border-radius: 8px;",
              div(
                class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
                div(
                  h6(tags$b("Online Documentation, Tutorials & Guides"), style = "color: #0369A1; margin: 0; font-size: 0.95rem;"),
                  p("Access comprehensive Vignettes, tutorials, and function reference at GitHub Pages:", style = "font-size: 0.82rem; color: #0284C7; margin: 2px 0 0 0;")
                ),
                tags$a(
                  href = "https://kabilanbio.github.io/scSimEval",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  class = "btn btn-sm btn-outline-primary",
                  style = "font-weight: 600;",
                  tags$i(class = "fa fa-external-link-alt me-1"), "https://kabilanbio.github.io/scSimEval"
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 6: Computational Sizing & Guidelines
            # ------------------------------------------------------------------
            div(
              style = "margin-top: 22px; padding: 20px 22px; background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
              div(
                class = "d-flex justify-content-between align-items-center mb-2",
                h5(tags$b("Recommended Dataset Sizes"), style = "color: #1E3A8A; margin: 0; font-size: 1.00rem;"),
                span(class = "badge bg-light text-primary border", "Cloud vs Local Workstation")
              ),
              p(
                "Choose the right environment depending on your dataset size (number of cells and genes):",
                style = "color: #475569; font-size: 0.85rem; margin-bottom: 12px;"
              ),
              tags$div(
                class = "table-responsive",
                tags$table(
                  class = "table table-bordered table-sm",
                  style = "font-size: 0.84rem; vertical-align: middle;",
                  tags$thead(
                    style = "background-color: #F8FAFC;",
                    tags$tr(
                      tags$th(style = "width: 22%;", "Environment"),
                      tags$th(style = "width: 18%;", "Hardware"),
                      tags$th(style = "width: 32%;", "Recommended Size"),
                      tags$th(style = "width: 28%;", "Best For")
                    )
                  ),
                  tags$tbody(
                    tags$tr(
                      tags$td(tags$b("Hosted Web App (Cloud)")),
                      tags$td("4 GB RAM, 1 vCPU"),
                      tags$td(HTML("scRNA: &le; 4,000 cells &times; 3,000 genes<br>Multiomics: &le; 3,000 cells &times; 8,000 peaks")),
                      tags$td("Quick exploratory analysis & example benchmark without installing R.")
                    ),
                    tags$tr(
                      tags$td(tags$b("Local R Package")),
                      tags$td("Your Workstation"),
                      tags$td("Scalable to 50,000+ cells and full peaksets"),
                      tags$td("Large multiomics datasets and high-throughput benchmarking.")
                    )
                  )
                )
              ),
              tags$div(
                style = "margin-top: 12px; background: #F8FAFC; border-radius: 6px; padding: 10px 14px; border: 1px solid #E2E8F0;",
                span(tags$b("Launch locally in R: "), style = "font-size: 0.82rem; color: #334155;"),
                tags$code("remotes::install_github('kabilanbio/scSimEval'); scSimEval::launch_scSimEval_app()", style = "font-size: 0.82rem; background: #FFFFFF; padding: 2px 6px; border: 1px solid #CBD5E1; border-radius: 4px; margin-left: 6px;")
              ),
              div(
                style = "margin-top: 12px; padding: 10px 14px; background: #F8FAFC; border-left: 3px solid #6366F1; border-radius: 4px; font-size: 0.82rem; color: #475569;",
                tags$b("Multiomics Note: "),
                "Multiomics evaluates 4 count matrices concurrently (Real RNA+ATAC, Sim RNA+ATAC). On the cloud, filtering peaksets to the top 5,000 - 8,000 peaks is recommended for best speed."
              )
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 7: CONTACT (SINGLE ICAR LOGO, PRINCIPAL SCIENTIST, KABILAN SAKTHIVEL)
  # ============================================================================
  nav_panel(
    "Contact",
    fluidRow(
      column(
        12,
        card(
          card_body(
            # Institute Header Section (Single ICAR Logo)
            div(
              style = "margin-bottom: 22px; padding: 18px 24px; background: #FFFFFF; border-radius: 10px; border: 1px solid #E2E8F0; box-shadow: 0 1px 4px rgba(0,0,0,0.03);",
              div(
                class = "d-flex align-items-center justify-content-center gap-4 flex-wrap",
                tags$img(src = "scfigures/icar_logo.jpg", height = 90, width = 90, style = "border-radius: 10px; object-fit: contain;"),
                div(
                  class = "text-center text-md-start",
                  h4("Division of Agricultural Bioinformatics", style = "font-size: 1.35rem; color: #1E3A8A; margin-bottom: 3px; font-weight: 700;"),
                  h5("ICAR - Indian Agricultural Statistics Research Institute (IASRI)", style = "font-size: 1.10rem; color: #334155; margin-bottom: 3px; font-weight: 600;"),
                  p("New Delhi, India", style = "font-size: 0.95rem; color: #64748B; margin-bottom: 0; font-weight: 500;")
                )
              )
            ),
            
            # Author Team Grid: Row 1 (3 Authors)
            div(
              class = "row justify-content-center g-3 my-2",
              
              # 1. Kabilan Sakthivel
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/kabilan.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Kabilan Sakthivel", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Ph.D. Bioinformatics", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:kabilan151414@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "kabilan151414@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 2. Dr Dwijesh Chandra Mishra
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/mishra_sir.jpg", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Dwijesh Chandra Mishra", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:dwij.mishra@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "dwij.mishra@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 3. Dr Shashi Bhushan Lal
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sb_lal_sir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Shashi Bhushan Lal", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:sblall16@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "sblall16@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              )
            ),
            
            # Author Team Grid: Row 2 (3 Authors)
            div(
              class = "row justify-content-center g-3 my-2",
              
              # 4. Dr Sudhir Srivastava
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sudhir_sir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Sudhir Srivastava", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Senior Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:sudhir0401bm@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "sudhir0401bm@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 5. Dr Krishna Kumar Chaturvedi
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/kkcsir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Krishna Kumar Chaturvedi", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:kkcchaturvedi@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "kkcchaturvedi@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 6. Dr Sharanbasappa
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sharan_photo.jpg", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Sharanbasappa", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:smadival509@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "smadival509@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              )
            ),
            
            # Feedback & GitHub Box
            div(
              style = "text-align: center; font-size: 0.88rem; color: #334155; margin: 26px auto 16px auto; max-width: 850px; padding: 14px 18px; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px;",
              icon("github", class = "me-2 text-dark"),
              strong("For feedback, questions, or bug reports: "),
              "visit our GitHub repository at ",
              tags$a(href = "https://github.com/kabilanbio/scSimEval", target = "_blank", style = "color: #0284C7; text-decoration: none; font-weight: 600;", "github.com/kabilanbio/scSimEval")
            ),
            
            # Footer Note
            div(
              style = "font-size: 0.80rem; text-align: center; color: #64748B; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px; margin-top: 10px;",
              span(icon("microscope", class = "me-1 text-primary"), tags$strong("Division of Agricultural Bioinformatics, ICAR-IASRI, New Delhi, India."))
            )
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER DEFINITION
# ==============================================================================
server <- function(input, output, session) {
  
  # Reactive state container
  # Reactive state container (Data loaded on-demand via options on sidebar)
  rv <- reactiveValues(
    benchmark_df = NULL,
    methods = NULL,
    toy_ref = NULL,
    toy_sim = NULL,
    sim_matrices = list(),
    cell_types = NULL,
    batch = NULL,
    source_name = "No Data Loaded",
    dataset_summary_df = NULL
  )
  
  # Navigation triggers
  observeEvent(input$btn_go_data, { nav_select("nav_active", "Data Hub") })
  observeEvent(input$btn_go_bubble, { nav_select("nav_active", "Comparative Bubble Matrix") })
  observeEvent(input$btn_go_viz, { nav_select("nav_active", "Visualizations") })
  observeEvent(input$btn_go_download, { nav_select("nav_active", "Download Results") })
  observeEvent(input$btn_go_help, { nav_select("nav_active", "Help & Getting Started") })
  observeEvent(input$btn_go_contact, { nav_select("nav_active", "Contact") })
  
  # ----------------------------------------------------------------------------
  # Google Analytics: Live / Cached Page Views Counter
  # ----------------------------------------------------------------------------
  output$total_pageviews <- renderText({
    # 1. Baseline cumulative all-time views (from launch through today)
    baseline_views <- 18
    
    # Global process environment to track cumulative sessions across container lifetime
    if (!exists(".scSimEval_views_env", envir = .GlobalEnv)) {
      assign(".scSimEval_views_env", new.env(parent = emptyenv()), envir = .GlobalEnv)
      .scSimEval_views_env$session_count <- 0
      .scSimEval_views_env$ga_views <- 0
      .scSimEval_views_env$last_ga_check <- 0
    }
    
    # Increment session count within this running application instance
    .scSimEval_views_env$session_count <- .scSimEval_views_env$session_count + 1
    
    # Read persisted cache file (if present)
    cached_file_views <- baseline_views
    cache_path <- "pageviews_cache.rds"
    if (file.exists(cache_path)) {
      try({
        cached <- readRDS(cache_path)
        if (is.numeric(cached) && cached >= baseline_views) {
          cached_file_views <- cached
        }
      }, silent = TRUE)
    }
    
    # 2. Query Google Analytics 4 (throttled to once every 5 minutes to avoid UI lag)
    now_ts <- as.numeric(Sys.time())
    if ((now_ts - .scSimEval_views_env$last_ga_check) > 300) {
      .scSimEval_views_env$last_ga_check <- now_ts
      
      ga_key_content <- Sys.getenv("GA_KEY_JSON", "")
      ga_key <- ""
      if (nzchar(ga_key_content)) {
        tmp_key <- tempfile(fileext = ".json")
        try(writeLines(ga_key_content, tmp_key), silent = TRUE)
        ga_key <- tmp_key
      } else {
        ga_key <- Sys.getenv("GA_AUTH_FILE", "")
        if (!nzchar(ga_key)) {
          candidates <- c(
            "google_key.json",
            file.path("www", "google_key.json"),
            file.path("inst", "shiny", "scSimEvalApp", "google_key.json"),
            file.path("inst", "shiny", "scSimEvalApp", "www", "google_key.json"),
            file.path("..", "google_key.json"),
            file.path("..", "..", "google_key.json")
          )
          for (cand in candidates) {
            if (file.exists(cand)) {
              ga_key <- cand
              break
            }
          }
        }
      }
      
      prop_id <- Sys.getenv("GA_PROPERTY_ID", "557610038")
      
      if (nzchar(ga_key) && file.exists(ga_key) && requireNamespace("googleAnalyticsR", quietly = TRUE)) {
        tryCatch({
          googleAnalyticsR::ga_auth(json_file = ga_key)
          if (!nzchar(prop_id)) {
            accs <- tryCatch(googleAnalyticsR::ga_account_list("ga4"), error = function(e) NULL)
            if (!is.null(accs) && nrow(accs) > 0 && "propertyId" %in% colnames(accs)) {
              prop_id <- as.character(accs$propertyId[1])
            }
          }
          if (nzchar(prop_id)) {
            df <- googleAnalyticsR::ga_data(
              propertyId = prop_id,
              date_range = c("2024-01-01", "today"),
              metrics = "screenPageViews"
            )
            if (!is.null(df) && nrow(df) > 0 && "screenPageViews" %in% colnames(df)) {
              val <- sum(as.numeric(df$screenPageViews), na.rm = TRUE)
              if (val > 0) {
                .scSimEval_views_env$ga_views <- val
              }
            }
          }
        }, error = function(e) NULL)
      }
    }
    
    # 3. Monotonically increasing cumulative total (never resets on container recycling)
    total_views <- max(
      baseline_views + .scSimEval_views_env$session_count - 1,
      cached_file_views + .scSimEval_views_env$session_count - 1,
      .scSimEval_views_env$ga_views + .scSimEval_views_env$session_count - 1
    )
    
    # Persist updated count to disk cache
    try(saveRDS(total_views, cache_path), silent = TRUE)
    
    formatC(total_views, format = "d", big.mark = ",")
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 1 - Load Demo Benchmark
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_load_demo, {
    if (!is.null(initial_demo)) {
      rv$benchmark_df <- initial_demo$benchmark_summary_table
      rv$methods <- initial_demo$methods
      rv$toy_ref <- initial_demo$toy_data$ref
      rv$toy_sim <- initial_demo$toy_data$sim
      rv$sim_matrices <- init_demo_sim_matrices(initial_demo)
      rv$cell_types <- initial_demo$toy_data$cell_types
      rv$batch <- initial_demo$toy_data$batch
      rv$source_name <- "Demo Benchmark (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC)"
      rv$dataset_summary_df <- init_demo_summary(initial_demo)
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      showNotification("Demo benchmark loaded successfully!", type = "message")
    } else {
      showNotification("Demo benchmark file not found on disk.", type = "warning")
    }
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 2 Dynamic Inputs (Single-Cell: scRNA-seq / scATAC-seq)
  # ----------------------------------------------------------------------------
  output$ui_uni_sim_inputs <- renderUI({
    n_sims <- if (!is.null(input$num_uni_sims)) as.integer(input$num_uni_sims) else 1
    if (is.na(n_sims) || n_sims < 1) n_sims <- 1
    if (n_sims > 6) n_sims <- 6
    
    inputs_list <- lapply(seq_len(n_sims), function(i) {
      default_name <- paste("Simulator", i)
      div(
        class = "sim-input-card",
        tags$b(paste0("Simulator ", i, ": "), style = "font-size: 0.9rem; color: #1B4F72;"),
        textInput(paste0("uni_sim_name_", i), "Simulator Name:", value = default_name),
        fileInput(paste0("file_uni_sim_", i), "Simulated Count Matrix (.rds / .csv / .tsv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
        fluidRow(
          column(6, numericInput(paste0("uni_sim_time_", i), "Elapsed Time (s):", value = round(25 + i * 10, 1), min = 0.1, step = 0.5)),
          column(6, numericInput(paste0("uni_sim_mem_", i), "Peak RAM (MB):", value = round(650 + i * 150, 0), min = 1, step = 10))
        )
      )
    })
    tagList(inputs_list)
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 2 Evaluation Trigger (Single-Cell: 1 or Multiple Simulators)
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_run_uni_eval, {
    req(input$file_uni_ref)
    n_sims <- if (!is.null(input$num_uni_sims)) as.integer(input$num_uni_sims) else 1
    if (is.na(n_sims) || n_sims < 1) n_sims <- 1
    
    withProgress(message = "Single-Cell Evaluation", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Loading biological reference matrix...")
        ref_mat <- read_uploaded_matrix(input$file_uni_ref$datapath, input$file_uni_ref$name)
        
        cell_types_vec <- read_uploaded_labels(input$file_uni_celltypes$datapath, input$file_uni_celltypes$name)
        batch_vec <- read_uploaded_labels(input$file_uni_batch$datapath, input$file_uni_batch$name)
        
        results_list <- list()
        summary_rows <- list()
        all_sim_mats <- list()
        first_sim_mat <- NULL
        
        summary_rows[[1]] <- extract_dataset_summary(
          ref_mat, role = "Biological Reference", method_name = "Empirical Reference",
          modality = "Single-Cell (Counts)", cell_types = cell_types_vec, batch_info = batch_vec
        )
        
        for (i in seq_len(n_sims)) {
          sim_file <- input[[paste0("file_uni_sim_", i)]]
          if (is.null(sim_file)) {
            stop(sprintf("Please upload the simulated count matrix for Simulator %d.", i))
          }
          
          sim_name <- input[[paste0("uni_sim_name_", i)]]
          if (is.null(sim_name) || trimws(sim_name) == "") {
            sim_name <- tools::file_path_sans_ext(sim_file$name)
            if (is.null(sim_name) || trimws(sim_name) == "") sim_name <- paste("Simulator", i)
          }
          
          sim_time <- as.numeric(input[[paste0("uni_sim_time_", i)]])
          if (is.null(sim_time) || is.na(sim_time)) sim_time <- 30.0
          
          sim_mem <- as.numeric(input[[paste0("uni_sim_mem_", i)]])
          if (is.null(sim_mem) || is.na(sim_mem)) sim_mem <- 800.0
          
          incProgress(0.7 / n_sims, detail = sprintf("Evaluating [%d/%d]: %s", i, n_sims, sim_name))
          
          sim_mat <- read_uploaded_matrix(sim_file$datapath, sim_file$name)
          if (i == 1) first_sim_mat <- sim_mat
          all_sim_mats[[sim_name]] <- sim_mat
          
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            sim_mat, role = "Simulated", method_name = sim_name,
            modality = "Single-Cell (Counts)", cell_types = cell_types_vec, batch_info = batch_vec
          )
          
          res_i <- evaluate_simulation_accuracy(
            ref_data = ref_mat,
            sim_data = sim_mat,
            memory_mb = sim_mem,
            elapsed_time = sim_time,
            compute_bivariate = FALSE,
            verbose = FALSE
          )
          
          tbl_i <- res_i$metrics_summary_table
          tbl_i$Method <- sim_name
          if (!"Score" %in% colnames(tbl_i) && "Value" %in% colnames(tbl_i)) {
            tbl_i$Score <- tbl_i$Value
          }
          tbl_i <- standardize_benchmark_categories(tbl_i)
          results_list[[i]] <- tbl_i
        }
        
        incProgress(0.1, detail = "Consolidating evaluated datasets...")
        combined_df <- do.call(rbind, results_list)
        
        if (isTRUE(input$chk_append_uni) && !is.null(rv$benchmark_df)) {
          existing_clean <- rv$benchmark_df[!rv$benchmark_df$Method %in% unique(combined_df$Method), , drop = FALSE]
          rv$benchmark_df <- rbind(existing_clean, combined_df[, intersect(colnames(existing_clean), colnames(combined_df))])
        } else {
          rv$benchmark_df <- combined_df
        }
        
        new_summary_df <- do.call(rbind, summary_rows)
        if (isTRUE(input$chk_append_uni) && !is.null(rv$dataset_summary_df)) {
          existing_sum_clean <- rv$dataset_summary_df[!rv$dataset_summary_df[["Dataset / Simulator"]] %in% new_summary_df[["Dataset / Simulator"]], , drop = FALSE]
          rv$dataset_summary_df <- rbind(existing_sum_clean, new_summary_df)
        } else {
          rv$dataset_summary_df <- new_summary_df
        }
        
        rv$methods <- unique(rv$benchmark_df$Method)
        rv$toy_ref <- ref_mat
        rv$toy_sim <- first_sim_mat
        rv$sim_matrices <- all_sim_mats
        rv$cell_types <- cell_types_vec
        rv$batch <- batch_vec
        rv$source_name <- sprintf("Single-Cell Benchmark (%d Simulators)", length(rv$methods))
        
        updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
        showNotification(sprintf("Single-cell evaluation complete! Evaluated %d simulator(s).", n_sims), type = "message")
      }, error = function(e) {
        showNotification(paste("Evaluation error:", e$message), type = "error")
      })
    })
  })
  
  output$ui_pairing_info_banner <- renderUI({
    if (identical(input$opt_multi_pairing, "unpaired")) {
      div(
        class = "alert alert-warning py-2 px-3 mb-3",
        style = "font-size: 0.82rem; border-left: 4px solid #D97706;",
        tags$b(icon("info-circle"), " Unpaired Multiomics Mode Active:"),
        p("Evaluates unimodal fidelity per layer plus population-level cross-modal metrics (label transfer, peak co-accessibility, and gene co-expression modularity). Cell-level pairing metrics (FOSCTTM, Match@1, Cross-Modal Generation, Peak-to-Gene Coupling) are excluded as cell barcodes are unlinked.", style = "margin-bottom: 0;")
      )
    } else {
      div(
        class = "alert alert-info py-2 px-3 mb-3",
        style = "font-size: 0.82rem; border-left: 4px solid #2563EB;",
        tags$b(icon("check-circle"), " Paired Multiomics Mode Active:"),
        p("Simultaneous co-assay from identical cells. Evaluates all 62 measures across all 8 canonical categories, including FOSCTTM, Match@1, Cross-Modal Generation, and Direct Regulatory Coupling.", style = "margin-bottom: 0;")
      )
    }
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 3 Dynamic Inputs (Multiomics: scRNA-seq + scATAC-seq)
  # ----------------------------------------------------------------------------
  output$ui_multiomics_sim_inputs <- renderUI({
    n_sims <- if (!is.null(input$num_multi_sims)) as.integer(input$num_multi_sims) else 1
    if (is.na(n_sims) || n_sims < 1) n_sims <- 1
    if (n_sims > 5) n_sims <- 5
    
    inputs_list <- lapply(seq_len(n_sims), function(i) {
      default_name <- paste("MultiSimulator", i)
      div(
        class = "sim-input-card",
        tags$b(paste0("Simulator ", i, ": "), style = "font-size: 0.9rem; color: #1B4F72;"),
        textInput(paste0("multi_sim_name_", i), "Simulator Name:", value = default_name),
        fileInput(paste0("file_multi_sim_rna_", i), "Simulated RNA Matrix (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
        fileInput(paste0("file_multi_sim_atac_", i), "Simulated ATAC Matrix (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
        fluidRow(
          column(6, numericInput(paste0("multi_sim_time_", i), "Elapsed Time (s):", value = round(45 + i * 15, 1), min = 0.1, step = 0.5)),
          column(6, numericInput(paste0("multi_sim_mem_", i), "Peak RAM (MB):", value = round(1100 + i * 250, 0), min = 1, step = 10))
        )
      )
    })
    tagList(inputs_list)
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 3 Evaluation Trigger (Multiomics Simulators)
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_run_multi_eval, {
    req(input$file_multi_ref_rna, input$file_multi_ref_atac)
    n_sims <- if (!is.null(input$num_multi_sims)) as.integer(input$num_multi_sims) else 1
    if (is.na(n_sims) || n_sims < 1) n_sims <- 1
    
    withProgress(message = "Multiomics Evaluation", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Loading biological reference RNA and ATAC...")
        ref_rna <- read_uploaded_matrix(input$file_multi_ref_rna$datapath, input$file_multi_ref_rna$name)
        ref_atac <- read_uploaded_matrix(input$file_multi_ref_atac$datapath, input$file_multi_ref_atac$name)
        
        cell_types_vec <- read_uploaded_labels(input$file_multi_celltypes$datapath, input$file_multi_celltypes$name)
        batch_vec <- read_uploaded_labels(input$file_multi_batch$datapath, input$file_multi_batch$name)
        
        results_list <- list()
        summary_rows <- list()
        all_sim_rna_mats <- list()
        first_sim_rna <- NULL
        
        summary_rows[[1]] <- extract_dataset_summary(
          ref_rna, role = "Biological Reference", method_name = "Empirical Reference (RNA)",
          modality = "scRNA-seq", cell_types = cell_types_vec, batch_info = batch_vec
        )
        summary_rows[[2]] <- extract_dataset_summary(
          ref_atac, role = "Biological Reference", method_name = "Empirical Reference (ATAC)",
          modality = "scATAC-seq", cell_types = cell_types_vec, batch_info = batch_vec
        )
        
        for (i in seq_len(n_sims)) {
          rna_file <- input[[paste0("file_multi_sim_rna_", i)]]
          atac_file <- input[[paste0("file_multi_sim_atac_", i)]]
          
          if (is.null(rna_file) || is.null(atac_file)) {
            stop(sprintf("Please upload both simulated RNA and ATAC files for Simulator %d.", i))
          }
          
          sim_name <- input[[paste0("multi_sim_name_", i)]]
          if (is.null(sim_name) || trimws(sim_name) == "") {
            sim_name <- paste("MultiSimulator", i)
          }
          
          sim_time <- as.numeric(input[[paste0("multi_sim_time_", i)]])
          if (is.null(sim_time) || is.na(sim_time)) sim_time <- 45.0
          
          sim_mem <- as.numeric(input[[paste0("multi_sim_mem_", i)]])
          if (is.null(sim_mem) || is.na(sim_mem)) sim_mem <- 1200.0
          
          incProgress(0.7 / n_sims, detail = sprintf("Evaluating multiomics [%d/%d]: %s", i, n_sims, sim_name))
          
          sim_rna <- read_uploaded_matrix(rna_file$datapath, rna_file$name)
          sim_atac <- read_uploaded_matrix(atac_file$datapath, atac_file$name)
          if (i == 1) first_sim_rna <- sim_rna
          all_sim_rna_mats[[sim_name]] <- sim_rna
          
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            sim_rna, role = "Simulated", method_name = paste0(sim_name, " (RNA)"),
            modality = "scRNA-seq", cell_types = cell_types_vec, batch_info = batch_vec
          )
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            sim_atac, role = "Simulated", method_name = paste0(sim_name, " (ATAC)"),
            modality = "scATAC-seq", cell_types = cell_types_vec, batch_info = batch_vec
          )
          
          pairing_mode <- if (!is.null(input$opt_multi_pairing)) input$opt_multi_pairing else "paired"
          
          res_i <- evaluate_multiomics_accuracy(
            ref_multi = list(rna = ref_rna, atac = ref_atac),
            sim_multi = list(rna = sim_rna, atac = sim_atac),
            cell_types = cell_types_vec,
            batch_info = batch_vec,
            pairing = pairing_mode,
            memory_mb = sim_mem,
            elapsed_time = sim_time,
            compute_bivariate = FALSE,
            verbose = FALSE
          )
          
          tbl_i <- res_i$benchmark_summary_table
          tbl_i$Method <- sim_name
          if (!"Score" %in% colnames(tbl_i) && "Value" %in% colnames(tbl_i)) {
            tbl_i$Score <- tbl_i$Value
          }
          tbl_i <- standardize_benchmark_categories(tbl_i)
          results_list[[i]] <- tbl_i
        }
        
        incProgress(0.1, detail = "Consolidating multiomics benchmark...")
        combined_df <- do.call(rbind, results_list)
        
        if (isTRUE(input$chk_append_multi) && !is.null(rv$benchmark_df)) {
          existing_clean <- rv$benchmark_df[!rv$benchmark_df$Method %in% unique(combined_df$Method), , drop = FALSE]
          rv$benchmark_df <- rbind(existing_clean, combined_df[, intersect(colnames(existing_clean), colnames(combined_df))])
        } else {
          rv$benchmark_df <- combined_df
        }
        
        new_summary_df <- do.call(rbind, summary_rows)
        if (isTRUE(input$chk_append_multi) && !is.null(rv$dataset_summary_df)) {
          existing_sum_clean <- rv$dataset_summary_df[!rv$dataset_summary_df[["Dataset / Simulator"]] %in% new_summary_df[["Dataset / Simulator"]], , drop = FALSE]
          rv$dataset_summary_df <- rbind(existing_sum_clean, new_summary_df)
        } else {
          rv$dataset_summary_df <- new_summary_df
        }
        
        rv$methods <- unique(rv$benchmark_df$Method)
        rv$toy_ref <- ref_rna
        rv$toy_sim <- first_sim_rna
        rv$sim_matrices <- all_sim_rna_mats
        rv$cell_types <- cell_types_vec
        rv$batch <- batch_vec
        rv$source_name <- sprintf("Multiomics Benchmark (%d Simulators, %s)", length(rv$methods), ifelse(pairing_mode == "paired", "Paired", "Unpaired"))
        
        updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
        showNotification(sprintf("Multiomics evaluation complete! Evaluated %d simulator(s).", n_sims), type = "message")
      }, error = function(e) {
        showNotification(paste("Multiomics evaluation error:", e$message), type = "error")
      })
    })
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 4 - Load Saved Benchmark File
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_load_uploaded_bench, {
    req(input$file_bench_upload)
    tryCatch({
      ext <- tolower(tools::file_ext(input$file_bench_upload$name))
      if (ext == "rds") {
        obj <- readRDS(input$file_bench_upload$datapath)
        if (is.data.frame(obj)) {
          df <- obj
        } else if (is.list(obj) && !is.null(obj$benchmark_summary_table)) {
          df <- obj$benchmark_summary_table
        } else if (is.list(obj) && !is.null(obj$metrics_summary_table)) {
          df <- obj$metrics_summary_table
        } else if (is.list(obj) && !is.null(obj$summary_table)) {
          df <- obj$summary_table
        } else {
          stop("Unrecognized RDS format. Must contain benchmark summary table.")
        }
        
        if (is.list(obj) && !is.null(obj$dataset_summary)) {
          rv$dataset_summary_df <- obj$dataset_summary
        } else if (is.list(obj) && !is.null(obj$dataset_properties)) {
          rv$dataset_summary_df <- obj$dataset_properties
        } else {
          rv$dataset_summary_df <- NULL
        }
      } else if (ext == "csv") {
        df <- utils::read.csv(input$file_bench_upload$datapath, check.names = FALSE)
        rv$dataset_summary_df <- NULL
      } else {
        stop("Unsupported file type. Please upload .rds or .csv.")
      }
      
      if (!"Method" %in% colnames(df) && "Simulator" %in% colnames(df)) df$Method <- df$Simulator
      if (!"Score" %in% colnames(df) && "Value" %in% colnames(df)) df$Score <- df$Value
      df <- standardize_benchmark_categories(df)
      
      rv$benchmark_df <- df
      rv$methods <- unique(df$Method)
      rv$source_name <- paste0("Uploaded File: ", input$file_bench_upload$name)
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      showNotification("Benchmark results loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: File Upload Badges & Dataset Properties Summary
  # ----------------------------------------------------------------------------
  output$ui_uni_ref_badge <- renderUI({
    req(input$file_uni_ref)
    tryCatch({
      mat <- read_uploaded_matrix(input$file_uni_ref$datapath, input$file_uni_ref$name)
      if (!is.null(mat)) {
        nc <- ncol(mat); nr <- nrow(mat)
        sp_pct <- if (inherits(mat, "sparseMatrix")) (1 - (Matrix::nnzero(mat) / (as.numeric(nc) * as.numeric(nr)))) * 100 else (sum(mat == 0) / (as.numeric(nc) * as.numeric(nr))) * 100
        div(class = "alert alert-light py-2 px-3 mb-2 border", style = "font-size: 0.82rem; background: #F8FAFC;",
            tags$b(icon("check-circle", class = "text-success"), " Uploaded Reference: "),
            sprintf("%s cells × %s features (Sparsity: %.1f%%)",
                    formatC(nc, format = "d", big.mark = ","),
                    formatC(nr, format = "d", big.mark = ","),
                    sp_pct))
      }
    }, error = function(e) NULL)
  })
  
  output$ui_multi_ref_rna_badge <- renderUI({
    req(input$file_multi_ref_rna)
    tryCatch({
      mat <- read_uploaded_matrix(input$file_multi_ref_rna$datapath, input$file_multi_ref_rna$name)
      if (!is.null(mat)) {
        nc <- ncol(mat); nr <- nrow(mat)
        sp_pct <- if (inherits(mat, "sparseMatrix")) (1 - (Matrix::nnzero(mat) / (as.numeric(nc) * as.numeric(nr)))) * 100 else (sum(mat == 0) / (as.numeric(nc) * as.numeric(nr))) * 100
        div(class = "alert alert-light py-2 px-3 mb-2 border", style = "font-size: 0.82rem; background: #F8FAFC;",
            tags$b(icon("check-circle", class = "text-success"), " Uploaded RNA: "),
            sprintf("%s cells × %s genes (Sparsity: %.1f%%)",
                    formatC(nc, format = "d", big.mark = ","),
                    formatC(nr, format = "d", big.mark = ","),
                    sp_pct))
      }
    }, error = function(e) NULL)
  })
  
  output$ui_multi_ref_atac_badge <- renderUI({
    req(input$file_multi_ref_atac)
    tryCatch({
      mat <- read_uploaded_matrix(input$file_multi_ref_atac$datapath, input$file_multi_ref_atac$name)
      if (!is.null(mat)) {
        nc <- ncol(mat); nr <- nrow(mat)
        sp_pct <- if (inherits(mat, "sparseMatrix")) (1 - (Matrix::nnzero(mat) / (as.numeric(nc) * as.numeric(nr)))) * 100 else (sum(mat == 0) / (as.numeric(nc) * as.numeric(nr))) * 100
        div(class = "alert alert-light py-2 px-3 mb-2 border", style = "font-size: 0.82rem; background: #F8FAFC;",
            tags$b(icon("check-circle", class = "text-success"), " Uploaded ATAC: "),
            sprintf("%s cells × %s peaks (Sparsity: %.1f%%)",
                    formatC(nc, format = "d", big.mark = ","),
                    formatC(nr, format = "d", big.mark = ","),
                    sp_pct))
      }
    }, error = function(e) NULL)
  })
  output$ui_dataset_summary_download_btn <- renderUI({
    if (!is.null(rv$dataset_summary_df) && nrow(rv$dataset_summary_df) > 0) {
      downloadButton("download_dataset_summary_csv", "Export Properties (CSV)", class = "btn btn-sm btn-outline-secondary", icon = icon("file-csv"))
    }
  })
  
  output$ui_dataset_summary_kpis <- renderUI({
    if (is.null(rv$dataset_summary_df) || nrow(rv$dataset_summary_df) == 0) {
      return(p("Dataset properties will appear here once matrices are ingested.", style = "font-size: 0.85rem; color: #94A3B8; font-style: italic;"))
    }
    df <- rv$dataset_summary_df
    ref_idx <- which(df$Role == "Biological Reference" | df$Role == "Reference")
    row_pick <- if (length(ref_idx) > 0) df[ref_idx[1], ] else df[1, ]
    
    feats_str <- if (length(ref_idx) >= 2) {
      paste0(df[ref_idx[1], "Features (P)"], " (RNA) + ", df[ref_idx[2], "Features (P)"], " (ATAC)")
    } else {
      row_pick[["Features (P)"]]
    }
    
    fluidRow(
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #1E3A8A; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Cells (N)"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Cells (N)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", feats_str),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Features (P)"))),
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Sparsity"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Sparsity (S%)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #6366F1; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Cell Types (Groups)"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Biological Groups (K)"))),
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0284C7; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Batches"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Batches (B)")))
    )
  })
  
  output$table_dataset_summary <- renderDT({
    req(rv$dataset_summary_df)
    df <- rv$dataset_summary_df
    bold_cols <- intersect(c("Dataset / Simulator", "Dataset", "Role", "Modality"), colnames(df))
    mono_cols <- intersect(c("Cells (N)", "Features (P)", "Sparsity", "Cell Types (Groups)", "Batches", "Median Lib Size", "Median Detected Features", "Mean Expression"), colnames(df))
    
    dt <- datatable(
      df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "t",
        autoWidth = TRUE
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    )
    if (length(bold_cols) > 0) {
      dt <- dt %>% formatStyle(columns = bold_cols, fontWeight = "bold")
    }
    if (length(mono_cols) > 0) {
      dt <- dt %>% formatStyle(columns = mono_cols, fontFamily = "monospace")
    }
    dt
  })
  
  output$download_dataset_summary_csv <- downloadHandler(
    filename = function() { paste0("dataset_properties_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv") },
    content = function(file) {
      req(rv$dataset_summary_df)
      utils::write.csv(rv$dataset_summary_df, file, row.names = FALSE)
    }
  )
  
  # ----------------------------------------------------------------------------
  # Data Hub: Main Panel UI & Dynamic Explorer
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_jump_to_bubble, {
    nav_select("nav_active", "Comparative Bubble Matrix")
  })

  output$ui_datahub_main_panel <- renderUI({
    if (is.null(rv$benchmark_df) && is.null(rv$toy_ref)) {
      # Clean Minimalist Empty State
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center",
            span(icon("database", class = "me-2 text-primary"), tags$strong("Data Hub")),
            span(class = "badge bg-light text-secondary border", "Choose an option on the left")
          )
        ),
        card_body(
          style = "padding: 38px 24px; text-align: center;",
          div(
            style = "max-width: 620px; margin: 0 auto;",
            div(
              icon("folder-open", class = "text-primary mb-3", style = "font-size: 2.8rem; opacity: 0.85;")
            ),
            h4("Welcome to Data Hub", style = "font-weight: 700; color: #1E293B; margin-bottom: 8px; font-size: 1.15rem;"),
            p(
              "Load our example benchmark to explore all features with 6 pre-evaluated simulators, or upload your own real and simulated count matrices using the controls on the left.",
              style = "font-size: 0.90rem; color: #64748B; line-height: 1.6; margin-bottom: 22px;"
            ),
            div(
              class = "d-flex justify-content-center gap-3 flex-wrap",
              actionButton("btn_empty_load_demo", "Load Example Benchmark (6 Simulators)", class = "btn btn-primary px-3 py-2", icon = icon("play")),
              actionButton("btn_empty_help", "Open User Guide", class = "btn btn-outline-secondary px-3 py-2", icon = icon("book-open"))
            )
          )
        )
      )
    } else {
      # Active Ingestion: Scientific Tabset Layout
      navset_card_tab(
        title = div(
          class = "d-flex align-items-center justify-content-between w-100",
          span(icon("database", class = "me-2 text-primary"), tags$strong("Loaded Datasets")),
          span(class = "badge bg-success-subtle text-success border border-success-subtle", style = "font-size: 0.78rem; font-weight: 600;", "Data Loaded")
        ),
        selected = "Dataset Properties",
        
        # Subtab 1: Dataset Properties
        nav_panel(
          title = span(icon("table-cells", class = "me-1"), "Dataset Properties"),
          div(
            style = "padding: 8px 4px;",
            uiOutput("ui_status_banner"),
            div(
              class = "dataset-summary-box mb-3",
              div(
                class = "d-flex justify-content-between align-items-center mt-2 mb-2",
                div(
                  h6(tags$b("Dataset Properties Summary"), style = "color: #0F172A; margin: 0; font-size: 0.95rem;"),
                  p("Number of cells, features, sparsity, and sequencing depth across real and simulated datasets:",
                    style = "font-size: 0.82rem; color: #64748B; margin: 0;")
                ),
                uiOutput("ui_dataset_summary_download_btn")
              ),
              uiOutput("ui_dataset_summary_kpis"),
              div(style = "margin-top: 14px; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px; background: #FFFFFF;",
                  DTOutput("table_dataset_summary"))
            ),
            uiOutput("ui_datahub_annotations_breakdown")
          )
        ),
        
        # Subtab 2: View Count Matrices
        nav_panel(
          title = span(icon("magnifying-glass", class = "me-1"), "View Count Matrices"),
          div(
            style = "padding: 8px 4px;",
            div(
              style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 14px; margin-bottom: 14px;",
              fluidRow(
                column(4, uiOutput("ui_inspect_dataset_selector")),
                column(
                  5,
                  radioButtons(
                    "opt_inspect_view_mode", "View:",
                    choices = c(
                      "Raw Counts" = "matrix_head",
                      "Gene Statistics" = "feature_stats",
                      "Cell Statistics" = "cell_stats"
                    ),
                    selected = "matrix_head",
                    inline = TRUE
                  )
                ),
                column(
                  3,
                  div(
                    style = "margin-top: 24px; text-align: right;",
                    downloadButton("download_inspect_data_csv", "Export Table (CSV)", class = "btn btn-sm btn-outline-secondary", icon = icon("file-csv"))
                  )
                )
              )
            ),
            uiOutput("ui_inspect_matrix_kpis"),
            div(
              style = "margin-top: 14px; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px; background: #FFFFFF;",
              DTOutput("table_inspect_matrix_view")
            )
          )
        ),
        
        # Subtab 3: Benchmark Scores Preview
        nav_panel(
          title = span(icon("chart-simple", class = "me-1"), "Benchmark Scores Preview"),
          div(
            style = "padding: 8px 4px;",
            div(
              class = "d-flex justify-content-between align-items-center mb-3",
              div(
                h6(tags$b("Benchmark Scores Summary"), style = "color: #0F172A; margin-bottom: 2px; font-size: 0.95rem;"),
                p("Fidelity scores between 0 and 1 for each evaluated simulation method.",
                  style = "font-size: 0.82rem; color: #64748B; margin: 0;")
              ),
              actionButton("btn_jump_to_bubble", "View Bubble Matrix (Tab 3) →", class = "btn btn-outline-primary btn-sm", icon = icon("arrow-right"))
            ),
            div(
              style = "border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px; background: #FFFFFF;",
              DTOutput("table_active_data_preview")
            )
          )
        )
      )
    }
  })

  output$ui_status_banner <- renderUI({
    if (is.null(rv$benchmark_df) && is.null(rv$toy_ref)) {
      return(NULL)
    }
    n_methods <- if (!is.null(rv$benchmark_df)) length(unique(rv$benchmark_df$Method)) else length(rv$methods)
    n_metrics <- if (!is.null(rv$benchmark_df)) length(unique(rv$benchmark_df$Metric)) else 0
    n_records <- if (!is.null(rv$benchmark_df)) nrow(rv$benchmark_df) else 0
    methods_str <- if (!is.null(rv$methods)) paste(rv$methods, collapse = ", ") else "None"
    
    div(
      style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-left: 4px solid #1E3A8A; border-radius: 6px; padding: 10px 14px; margin-bottom: 14px;",
      div(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        div(
          tags$span(style = "font-weight: 700; color: #0F172A; font-size: 0.92rem;", icon("circle-check", class = "text-success me-1"), "Active Benchmark: ", rv$source_name),
          tags$span(style = "font-size: 0.82rem; color: #64748B; margin-left: 12px;",
                    sprintf("Simulators (%d): %s", n_methods, methods_str))
        ),
        div(
          span(class = "badge bg-white text-secondary border", style = "font-size: 0.78rem;", sprintf("%s Records", formatC(n_records, format = "d", big.mark = ","))),
          span(class = "badge bg-white text-primary border ms-1", style = "font-size: 0.78rem;", sprintf("%d Evaluated Metrics", n_metrics))
        )
      )
    )
  })

  output$ui_datahub_annotations_breakdown <- renderUI({
    has_ct <- !is.null(rv$cell_types) && length(rv$cell_types) > 0
    has_bt <- !is.null(rv$batch) && length(rv$batch) > 0
    if (!has_ct && !has_bt) return(NULL)
    
    fluidRow(
      style = "margin-top: 14px;",
      if (has_ct) {
        ct_tbl <- as.data.frame(table(Cell_Type = rv$cell_types))
        ct_tbl$Percentage <- sprintf("%.1f%%", ct_tbl$Freq / sum(ct_tbl$Freq) * 100)
        colnames(ct_tbl) <- c("Cell Type / Cluster", "Cell Count", "Percentage")
        column(
          if (has_bt) 6 else 12,
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 12px;",
            h6(tags$b(icon("layer-group", class = "me-1"), " Biological Group Distribution"), style = "color: #1E3A8A; margin-bottom: 8px; font-size: 0.88rem;"),
            renderTable(ct_tbl, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")
          )
        )
      },
      if (has_bt) {
        bt_tbl <- as.data.frame(table(Batch = rv$batch))
        bt_tbl$Percentage <- sprintf("%.1f%%", bt_tbl$Freq / sum(bt_tbl$Freq) * 100)
        colnames(bt_tbl) <- c("Technical Batch", "Cell Count", "Percentage")
        column(
          if (has_ct) 6 else 12,
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 12px;",
            h6(tags$b(icon("boxes-stacked", class = "me-1"), " Technical Batch Distribution"), style = "color: #0D9488; margin-bottom: 8px; font-size: 0.88rem;"),
            renderTable(bt_tbl, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")
          )
        )
      }
    )
  })

  output$ui_inspect_dataset_selector <- renderUI({
    choices <- list()
    if (!is.null(rv$toy_ref)) {
      choices[["Empirical Reference (Real Cells)"]] <- "ref"
    }
    if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      for (sname in names(rv$sim_matrices)) {
        choices[[paste0("Simulated: ", sname)]] <- paste0("sim_", sname)
      }
    } else if (!is.null(rv$toy_sim)) {
      choices[["Simulated Matrix"]] <- "toy_sim"
    }
    if (length(choices) == 0) choices[["No Matrix Available"]] <- "none"
    selectInput("sel_inspect_matrix_choice", "Select Dataset / Matrix:", choices = choices, selected = choices[[1]])
  })

  active_inspect_data <- reactive({
    choice <- input$sel_inspect_matrix_choice
    mode <- input$opt_inspect_view_mode
    
    mat <- NULL
    if (identical(choice, "ref") && !is.null(rv$toy_ref)) {
      mat <- rv$toy_ref
    } else if (identical(choice, "toy_sim") && !is.null(rv$toy_sim)) {
      mat <- rv$toy_sim
    } else if (!is.null(choice) && startsWith(choice, "sim_") && !is.null(rv$sim_matrices)) {
      sname <- sub("^sim_", "", choice)
      mat <- rv$sim_matrices[[sname]]
    } else if (!is.null(rv$toy_ref)) {
      mat <- rv$toy_ref
    } else if (!is.null(rv$toy_sim)) {
      mat <- rv$toy_sim
    }
    if (is.null(mat)) return(NULL)
    
    if (identical(mode, "matrix_head")) {
      nr <- min(50, nrow(mat))
      nc <- min(25, ncol(mat))
      sub_mat <- as.matrix(mat[1:nr, 1:nc, drop = FALSE])
      df <- as.data.frame(sub_mat)
      feat_names <- rownames(sub_mat)
      if (is.null(feat_names)) feat_names <- paste0("Feature_", seq_len(nr))
      df <- cbind(Feature_ID = feat_names, df)
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "matrix_head")
    } else if (identical(mode, "feature_stats")) {
      mean_expr <- Matrix::rowMeans(mat)
      det_rate <- Matrix::rowMeans(mat > 0)
      if (inherits(mat, "dgCMatrix")) {
        sq_means <- Matrix::rowMeans(mat^2)
        vars <- (sq_means - mean_expr^2) * (ncol(mat) / max(1, ncol(mat) - 1))
        vars[vars < 0] <- 0
      } else {
        vars <- apply(mat, 1, stats::var)
      }
      feat_names <- rownames(mat)
      if (is.null(feat_names)) feat_names <- paste0("Feature_", seq_len(nrow(mat)))
      df <- data.frame(
        Feature = feat_names,
        Mean_Expression = round(as.numeric(mean_expr), 4),
        Variance = round(as.numeric(vars), 4),
        Detection_Rate_Pct = round(as.numeric(det_rate) * 100, 2),
        Dropout_Rate_Pct = round((1 - as.numeric(det_rate)) * 100, 2),
        stringsAsFactors = FALSE
      )
      df <- df[order(-df$Mean_Expression), ]
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "feature_stats")
    } else {
      lib_sizes <- Matrix::colSums(mat)
      det_genes <- Matrix::colSums(mat > 0)
      cell_names <- colnames(mat)
      if (is.null(cell_names)) cell_names <- paste0("Cell_", seq_len(ncol(mat)))
      df <- data.frame(
        Cell_ID = cell_names,
        Library_Size = as.numeric(lib_sizes),
        Detected_Features = as.integer(det_genes),
        Sparsity_Pct = round((1 - as.numeric(det_genes) / max(1, nrow(mat))) * 100, 2),
        stringsAsFactors = FALSE
      )
      if (!is.null(rv$cell_types) && length(rv$cell_types) == ncol(mat)) {
        df$Cell_Type <- rv$cell_types
      }
      if (!is.null(rv$batch) && length(rv$batch) == ncol(mat)) {
        df$Batch <- rv$batch
      }
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "cell_stats")
    }
  })

  output$ui_inspect_matrix_kpis <- renderUI({
    res <- active_inspect_data()
    if (is.null(res) || is.null(res$mat)) return(NULL)
    mat <- res$mat
    n_cells <- ncol(mat)
    n_feats <- nrow(mat)
    sparsity_val <- round(mean(mat == 0, na.rm = TRUE) * 100, 1)
    med_depth <- round(stats::median(Matrix::colSums(mat)))
    
    fluidRow(
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #1E3A8A; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(n_cells, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Inspected Cells (N)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(n_feats, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Features (P)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", sprintf("%.1f%%", sparsity_val)),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Matrix Sparsity (S%)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #6366F1; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(med_depth, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Median Library Depth (L~)")))
    )
  })

  output$table_inspect_matrix_view <- renderDT({
    res <- active_inspect_data()
    req(res, res$df)
    dt <- datatable(
      res$df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "ftip"
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    )
    num_cols <- which(sapply(res$df, is.numeric))
    if (length(num_cols) > 0) {
      dt <- dt %>% formatStyle(num_cols, fontFamily = "monospace", textAlign = "right")
    }
    dt
  })

  output$download_inspect_data_csv <- downloadHandler(
    filename = function() {
      mode <- input$opt_inspect_view_mode %||% "view"
      paste0("scSimEval_inspect_", mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      res <- active_inspect_data()
      req(res, res$df)
      utils::write.csv(res$df, file, row.names = FALSE)
    }
  )

  output$table_active_data_preview <- renderDT({
    req(rv$benchmark_df)
    datatable(
      head(rv$benchmark_df, 100),
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
      rownames = FALSE,
      class = "compact stripe hover border"
    ) %>% formatRound(columns = which(sapply(head(rv$benchmark_df, 100), is.numeric)), digits = 4)
  })

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # Tab 3: Comparative Bubble Matrix & Leaderboard (Minimalist & User-Friendly)
  # ----------------------------------------------------------------------------
  output$ui_bubble_active_badge <- renderUI({
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      return(span(class = "badge bg-secondary", "No Metrics Selected"))
    }
    n_m <- length(unique(df$Metric))
    n_s <- length(unique(df$Method))
    mode <- input$rad_bubble_filter_mode %||% "all"
    lbl <- if (identical(mode, "custom")) "Custom Selection" else "All 62 Measures"
    span(class = "badge bg-primary", style = "font-size: 0.8rem; font-weight: 500;",
         sprintf("%s: %d Metrics across %d Simulators", lbl, n_m, n_s))
  })

  output$ui_bubble_method_picker <- renderUI({
    req(rv$methods)
    div(
      div(
        class = "d-flex justify-content-between align-items-center mb-1",
        tags$b("Simulators to Compare:"),
        div(
          actionLink("btn_bubble_methods_all", "All", style = "font-size: 0.78rem; margin-right: 6px; text-decoration: none; font-weight: 600;"),
          actionLink("btn_bubble_methods_none", "Clear", style = "font-size: 0.78rem; text-decoration: none; color: #dc3545;")
        )
      ),
      checkboxGroupInput(
        "sel_bubble_methods", NULL,
        choices = rv$methods,
        selected = rv$methods
      )
    )
  })

  observeEvent(input$btn_bubble_methods_all, {
    req(rv$methods)
    updateCheckboxGroupInput(session, "sel_bubble_methods", selected = rv$methods)
  })

  observeEvent(input$btn_bubble_methods_none, {
    updateCheckboxGroupInput(session, "sel_bubble_methods", selected = character(0))
  })

  output$ui_custom_category_picker <- renderUI({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cats <- if ("Category" %in% colnames(df)) {
      avail <- unique(df$Category)
      avail <- avail[!is.na(avail) & trimws(avail) != ""]
      canon_order <- c(
        "(I) Distributional Properties",
        "(II) Correlations & Zero-Inflation",
        "(III) Cellular Structure & Concordance",
        "(IV) Batch Effects & Confounder Mixing",
        "(V) Biological Signal & Downstream Fidelity",
        "(VI) Trajectory & Lineage Dynamics",
        "(VII) Cross-Modal Coupling & Modularity",
        "(VIII) Computational Scalability"
      )
      intersect(canon_order, avail)
    } else {
      character(0)
    }
    choices <- c("All Categories (All Measures)" = "all", cats)
    selectInput("sel_custom_cat", "1. Filter by Category:", choices = choices, selected = "all")
  })

  output$ui_custom_metrics_checklist <- renderUI({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) {
      df
    } else {
      df[df$Category == cat_chosen, , drop = FALSE]
    }
    
    avail_metrics <- unique(as.character(sub_df$Metric))
    avail_metrics <- avail_metrics[!is.na(avail_metrics) & trimws(avail_metrics) != ""]
    
    checkboxGroupInput(
      "chk_custom_metrics", NULL,
      choices = avail_metrics,
      selected = avail_metrics
    )
  })

  output$txt_custom_metrics_count <- renderText({
    n_sel <- length(input$chk_custom_metrics)
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) df else df[df$Category == cat_chosen, , drop = FALSE]
    n_cat <- length(unique(sub_df$Metric))
    sprintf("%d / %d metrics selected", n_sel, n_cat)
  })

  observeEvent(input$btn_custom_metrics_all, {
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) df else df[df$Category == cat_chosen, , drop = FALSE]
    avail_metrics <- unique(as.character(sub_df$Metric))
    updateCheckboxGroupInput(session, "chk_custom_metrics", selected = avail_metrics)
  })

  observeEvent(input$btn_custom_metrics_none, {
    updateCheckboxGroupInput(session, "chk_custom_metrics", selected = character(0))
  })

  filtered_bubble_data <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    # 1. Simulator Methods Filter
    if (!is.null(input$sel_bubble_methods)) {
      if (length(input$sel_bubble_methods) == 0) return(NULL)
      df <- df[df$Method %in% input$sel_bubble_methods, , drop = FALSE]
    }
    
    # 2. Mode Filter: 'all' vs 'custom'
    mode <- input$rad_bubble_filter_mode %||% "all"
    
    if (identical(mode, "custom")) {
      cat_chosen <- input$sel_custom_cat %||% "all"
      if (!identical(cat_chosen, "all") && "Category" %in% colnames(df)) {
        df <- df[df$Category == cat_chosen, , drop = FALSE]
      }
      
      if (!is.null(input$chk_custom_metrics)) {
        if (length(input$chk_custom_metrics) == 0) return(NULL)
        if ("Metric" %in% colnames(df)) {
          df <- df[df$Metric %in% input$chk_custom_metrics, , drop = FALSE]
        }
      }
    }
    
    if (nrow(df) == 0) return(NULL)
    df
  })

  bubble_plot_reactive <- reactive({
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      p_empty <- ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 1, y = 1,
          label = "No metrics or simulators selected.\nPlease select at least one simulator or check metrics in Custom Selection.",
          size = 5.2, color = "#64748B", fontface = "italic"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#F8FAFC", color = "#E2E8F0"))
      return(p_empty)
    }
    
    df <- standardize_benchmark_categories(df)
    
    p <- plot_benchmark_bubble_matrix(
      data              = df,
      base_size         = 11,
      compact_strips    = TRUE,
      show_missing_dots = FALSE,
      normalize_scores  = TRUE
    )
    
    n_metrics <- length(unique(df$Metric))
    n_cats <- length(unique(df$Category))
    n_sims <- length(unique(df$Method))
    mode <- input$rad_bubble_filter_mode %||% "all"
    mode_desc <- if (identical(mode, "custom")) "Custom Selection Benchmark" else "Full Benchmark (All 62 Canonical Measures)"
    
    p + ggplot2::labs(
      subtitle = sprintf("%s: %d Metric(s) across %d Category(ies) for %d Simulator(s)", mode_desc, n_metrics, n_cats, n_sims),
      caption = paste0(
        "Circle: standard performance (< 0.96)  |  Square: top performer (>= 0.96).\n",
        "All ", n_metrics, " metrics direction-normalized: for error/distance metrics, scores are inverted as 1 - norm(x) so 1.0 indicates closest agreement to empirical reference."
      )
    )
  })

  output$ui_bubble_plot_render <- renderUI({
    df <- filtered_bubble_data()
    n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
    default_h <- max(380, min(1200, 160 + n_m * 22))
    default_w <- if (n_m <= 10) 1400 else 2200
    h_val <- if (!is.null(input$sld_bubble_height)) input$sld_bubble_height else default_h
    w_val <- if (!is.null(input$sld_bubble_width)) input$sld_bubble_width else default_w
    plotOutput("plot_bubble_matrix", width = paste0(w_val, "px"), height = paste0(h_val, "px"))
  })

  output$plot_bubble_matrix <- renderPlot({
    bubble_plot_reactive()
  })

  leaderboard_reactive <- reactive({
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(
        Overall_Rank = integer(0),
        Method = character(0),
        Average_Fidelity = character(0),
        Fidelity_Score = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    compute_method_leaderboard(df)
  })

  output$table_leaderboard_dt <- renderDT({
    lb <- leaderboard_reactive()
    if (nrow(lb) == 0) {
      return(datatable(data.frame(Notice = "No data selected for ranking. Please select simulators or check metrics from the sidebar.")))
    }
    display_df <- lb
    colnames(display_df) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
    
    datatable(
      display_df,
      options = list(
        dom = "t",
        pageLength = 25,
        ordering = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    ) %>%
      formatStyle(
        "Overall Rank",
        fontWeight = "bold",
        backgroundColor = styleEqual(
          c(1, 2, 3),
          c("#FFF9DB", "#F1F3F5", "#FFF4E6")
        )
      ) %>%
      formatStyle(
        "Average Fidelity",
        fontWeight = "bold",
        color = "#1D72B8"
      ) %>%
      formatStyle(
        "Fidelity Score",
        fontFamily = "monospace",
        fontWeight = "bold"
      )
  })

  output$download_leaderboard_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_leaderboard_", Sys.Date(), ".csv") },
    content = function(file) {
      lb <- leaderboard_reactive()
      colnames(lb) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
      utils::write.csv(lb, file, row.names = FALSE)
    }
  )

  output$download_leaderboard_csv_side <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_leaderboard_", Sys.Date(), ".csv") },
    content = function(file) {
      lb <- leaderboard_reactive()
      colnames(lb) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
      utils::write.csv(lb, file, row.names = FALSE)
    }
  )

  output$download_bubble_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", Sys.Date(), ".jpeg") },
    content = function(file) {
      df <- filtered_bubble_data()
      n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
      h_in <- max(4.5, min(14, 2.5 + n_m * 0.16))
      w_in <- if (!is.null(input$sld_bubble_width)) max(10, input$sld_bubble_width / 95) else 23
      export_single_jpeg(file, bubble_plot_reactive(), width = w_in, height = h_in, dpi = 600)
    }
  )

  output$download_bubble_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", Sys.Date(), ".pdf") },
    content = function(file) {
      df <- filtered_bubble_data()
      n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
      h_in <- max(4.5, min(14, 2.5 + n_m * 0.16))
      w_in <- if (!is.null(input$sld_bubble_width)) max(10, input$sld_bubble_width / 95) else 23
      grDevices::pdf(file, width = w_in, height = h_in)
      print(bubble_plot_reactive())
      grDevices::dev.off()
    }
  )
  # Tab 4: Diagnostic Visualizations (7 Panels in 1 Row)
  # ----------------------------------------------------------------------------
  
  # 1. Evaluation Summary
  eval_summary_reactive <- reactive({
    req(rv$benchmark_df)
    plot_evaluation_summary(
      data = rv$benchmark_df,
      show_labels = input$chk_sum_labels,
      normalize_scores = input$chk_sum_norm,
      base_size = 14
    )
  })
  output$plot_eval_summary <- renderPlot({ eval_summary_reactive() })
  output$download_sum_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, eval_summary_reactive(), width = 13, height = 7.5, dpi = 600) }
  )
  output$download_sum_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7.5); print(eval_summary_reactive()); grDevices::dev.off() }
  )
  
  # 2. Distribution QC
  output$ui_dist_qc_plot <- renderUI({
    plot_h <- if (identical(input$sel_dist_layout, "comprehensive")) "850px" else "550px"
    plotOutput("plot_dist_qc", height = plot_h)
  })
  dist_qc_reactive <- reactive({
    req(rv$toy_ref, rv$toy_sim)
    plot_distribution_qc(
      ref_data = rv$toy_ref,
      sim_data = rv$toy_sim,
      layout   = input$sel_dist_layout,
      base_size = 13
    )
  })
  output$plot_dist_qc <- renderPlot({ dist_qc_reactive() })
  output$download_dist_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, dist_qc_reactive(), width = 14, height = 9, dpi = 600) }
  )
  output$download_dist_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 14, height = 9); print(dist_qc_reactive()); grDevices::dev.off() }
  )
  
  # 3. Scalability Benchmark
  scale_bench_reactive <- reactive({
    req(rv$benchmark_df)
    plot_scalability_benchmark(
      benchmark_data = rv$benchmark_df,
      type = input$sel_scale_type,
      base_size = 13
    )
  })
  output$plot_scale_bench <- renderPlot({ scale_bench_reactive() })
  output$download_scale_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, scale_bench_reactive(), width = 13, height = 8, dpi = 600) }
  )
  output$download_scale_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(scale_bench_reactive()); grDevices::dev.off() }
  )
  
  # 4. Metric Plots
  output$ui_box_metric_picker <- renderUI({
    req(rv$benchmark_df, input$sel_box_cat_first)
    sub_df <- rv$benchmark_df[rv$benchmark_df$Category == input$sel_box_cat_first, , drop = FALSE]
    avail_metrics <- sort(unique(sub_df$Metric))
    selectInput("sel_box_metric_single", "2. Choose Metric Name:", choices = avail_metrics, selected = avail_metrics[1])
  })
  
  metric_box_reactive <- reactive({
    req(rv$benchmark_df)
    if (identical(input$opt_box_view_mode, "individual")) {
      req(input$sel_box_metric_single)
      score_type_val <- if (!is.null(input$sel_box_indiv_score)) input$sel_box_indiv_score else "normalized"
      plot_individual_metric_bar(
        benchmark_data = rv$benchmark_df,
        metric         = input$sel_box_metric_single,
        score_type     = score_type_val,
        base_size      = 12
      )
    } else {
      if (identical(input$sel_box_cat_group, "all")) {
        plot_metric_boxplots(
          benchmark_data = rv$benchmark_df,
          categories     = NULL,
          score_type     = input$sel_box_score_type,
          facet_by       = "category",
          base_size      = 11
        )
      } else {
        plot_category_metric_bars(
          benchmark_data = rv$benchmark_df,
          category       = input$sel_box_cat_group,
          score_type     = input$sel_box_score_type,
          base_size      = 11
        )
      }
    }
  })

  output$ui_plot_metric_boxes <- renderUI({
    h <- 560
    if (identical(input$opt_box_view_mode, "individual")) {
      h <- 520
    } else {
      if (!identical(input$sel_box_cat_group, "all") && !is.null(rv$benchmark_df)) {
        sub_df <- rv$benchmark_df[rv$benchmark_df$Category == input$sel_box_cat_group, , drop = FALSE]
        n_m <- length(unique(sub_df$Metric))
        if (n_m <= 4) {
          h <- 440
        } else if (n_m <= 8) {
          h <- 640
        } else if (n_m <= 12) {
          h <- 820
        } else {
          h <- 960
        }
      } else {
        h <- 580
      }
    }
    plotOutput("plot_metric_boxes", height = paste0(h, "px"))
  })

  output$plot_metric_boxes <- renderPlot({ metric_box_reactive() })

  output$download_box_jpeg <- downloadHandler(
    filename = function() {
      if (identical(input$opt_box_view_mode, "individual")) {
        paste0("scSimEval_metric_bar_", gsub("[^A-Za-z0-9_-]", "_", input$sel_box_metric_single), "_", Sys.Date(), ".jpeg")
      } else {
        paste0("scSimEval_metric_plot_", Sys.Date(), ".jpeg")
      }
    },
    content = function(file) {
      w <- if (identical(input$opt_box_view_mode, "individual")) 10 else 13
      h <- if (identical(input$opt_box_view_mode, "individual")) 6.5 else 7.5
      export_single_jpeg(file, metric_box_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_box_pdf <- downloadHandler(
    filename = function() {
      if (identical(input$opt_box_view_mode, "individual")) {
        paste0("scSimEval_metric_bar_", gsub("[^A-Za-z0-9_-]", "_", input$sel_box_metric_single), "_", Sys.Date(), ".pdf")
      } else {
        paste0("scSimEval_metric_plot_", Sys.Date(), ".pdf")
      }
    },
    content = function(file) {
      w <- if (identical(input$opt_box_view_mode, "individual")) 10 else 13
      h <- if (identical(input$opt_box_view_mode, "individual")) 6.5 else 7.5
      grDevices::pdf(file, width = w, height = h)
      print(metric_box_reactive())
      grDevices::dev.off()
    }
  )

  get_box_cat_dims <- function() {
    if (identical(input$sel_box_cat_group, "all") || is.null(rv$benchmark_df)) {
      return(list(w = 13, h = 7.5))
    }
    sub_df <- rv$benchmark_df[rv$benchmark_df$Category == input$sel_box_cat_group, , drop = FALSE]
    n_m <- length(unique(sub_df$Metric))
    if (n_m <= 4) {
      list(w = 11, h = 5.5)
    } else if (n_m <= 8) {
      list(w = 13, h = 7.5)
    } else {
      list(w = 14, h = 10.5)
    }
  }

  output$download_box_cat_jpeg <- downloadHandler(
    filename = function() {
      if (identical(input$sel_box_cat_group, "all")) {
        paste0("scSimEval_category_boxplots_", Sys.Date(), ".jpeg")
      } else {
        cat_slug <- gsub("[^A-Za-z0-9_-]", "_", input$sel_box_cat_group)
        paste0("scSimEval_category_barplots_", cat_slug, "_", Sys.Date(), ".jpeg")
      }
    },
    content = function(file) {
      dims <- get_box_cat_dims()
      export_single_jpeg(file, metric_box_reactive(), width = dims$w, height = dims$h, dpi = 600)
    }
  )
  output$download_box_cat_pdf <- downloadHandler(
    filename = function() {
      if (identical(input$sel_box_cat_group, "all")) {
        paste0("scSimEval_category_boxplots_", Sys.Date(), ".pdf")
      } else {
        cat_slug <- gsub("[^A-Za-z0-9_-]", "_", input$sel_box_cat_group)
        paste0("scSimEval_category_barplots_", cat_slug, "_", Sys.Date(), ".pdf")
      }
    },
    content = function(file) {
      dims <- get_box_cat_dims()
      grDevices::pdf(file, width = dims$w, height = dims$h)
      print(metric_box_reactive())
      grDevices::dev.off()
    }
  )
  
  # 5. Metric Performance Heatmap (Interactive dimensions & Category filtering)
  observeEvent(input$sel_heat_cat, {
    if (identical(input$sel_heat_cat, "all")) {
      updateSliderInput(session, "sld_heat_height", value = 1100)
      updateSliderInput(session, "sld_heat_width", value = 880)
    } else {
      updateSliderInput(session, "sld_heat_height", value = 520)
      updateSliderInput(session, "sld_heat_width", value = 850)
    }
  }, ignoreInit = TRUE)

  output$ui_plot_metric_heat <- renderUI({
    w <- if (!is.null(input$sld_heat_width)) input$sld_heat_width else 880
    h <- if (!is.null(input$sld_heat_height)) input$sld_heat_height else 1100
    plotOutput("plot_metric_heat", width = paste0(w, "px"), height = paste0(h, "px"))
  })

  metric_heat_reactive <- reactive({
    req(rv$benchmark_df)
    cat_sel <- if (identical(input$sel_heat_cat, "all")) NULL else input$sel_heat_cat
    b_size <- if (is.null(cat_sel)) 9.5 else 11
    plot_metric_heatmap(
      benchmark_data = rv$benchmark_df,
      category = cat_sel,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      base_size = b_size
    )
  })
  output$plot_metric_heat <- renderPlot({ metric_heat_reactive() })
  output$download_heat_jpeg <- downloadHandler(
    filename = function() {
      cat_suffix <- if (identical(input$sel_heat_cat, "all")) "all_categories" else gsub("[^A-Za-z0-9_-]", "_", input$sel_heat_cat)
      paste0("scSimEval_metric_heatmap_", cat_suffix, "_", Sys.Date(), ".jpeg")
    },
    content = function(file) {
      is_all <- is.null(input$sel_heat_cat) || identical(input$sel_heat_cat, "all")
      w <- if (is_all) 9.5 else 9.0
      h <- if (is_all) 15.0 else 7.5
      export_single_jpeg(file, metric_heat_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_heat_pdf <- downloadHandler(
    filename = function() {
      cat_suffix <- if (identical(input$sel_heat_cat, "all")) "all_categories" else gsub("[^A-Za-z0-9_-]", "_", input$sel_heat_cat)
      paste0("scSimEval_metric_heatmap_", cat_suffix, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      is_all <- is.null(input$sel_heat_cat) || identical(input$sel_heat_cat, "all")
      w <- if (is_all) 9.5 else 9.0
      h <- if (is_all) 15.0 else 7.5
      grDevices::pdf(file, width = w, height = h)
      print(metric_heat_reactive())
      grDevices::dev.off()
    }
  )
  
  # 6. PCA Ordination (6 Category-wise options & All Categories Combined)
  observeEvent(input$sel_pca_panel, {
    if (identical(input$sel_pca_panel, "both")) {
      updateSliderInput(session, "sld_pca_height", value = 1100)
    } else {
      updateSliderInput(session, "sld_pca_height", value = 600)
    }
  }, ignoreInit = TRUE)

  output$ui_plot_metric_pca <- renderUI({
    w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width else 950
    h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height else 1100
    plotOutput("plot_metric_pca_out", width = paste0(w, "px"), height = paste0(h, "px"))
  })

  metric_pca_reactive <- reactive({
    req(rv$benchmark_df)
    cat_sel <- if (identical(input$sel_pca_cat, "all")) NULL else input$sel_pca_cat
    n_top <- if (!is.null(input$num_pca_top_metrics) && is.finite(input$num_pca_top_metrics)) input$num_pca_top_metrics else 14
    plot_metric_pca(
      benchmark_data = rv$benchmark_df,
      category = cat_sel,
      panel = input$sel_pca_panel,
      top_n_loadings = n_top,
      base_size = 12
    )
  })
  output$plot_metric_pca_out <- renderPlot({ metric_pca_reactive() })
  output$download_pca_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".jpeg") },
    content = function(file) {
      w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width / 90 else 10.5
      h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height / 90 else (if (identical(input$sel_pca_panel, "both")) 12.5 else 7.5)
      export_single_jpeg(file, metric_pca_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_pca_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".pdf") },
    content = function(file) {
      w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width / 90 else 10.5
      h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height / 90 else (if (identical(input$sel_pca_panel, "both")) 12.5 else 7.5)
      grDevices::pdf(file, width = w, height = h)
      print(metric_pca_reactive())
      grDevices::dev.off()
    }
  )
  
  # 7. MDS Metric Space (6 Category-wise options)
  output$ui_plot_metric_mds <- renderUI({
    w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width else 900
    h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height else 620
    plotOutput("plot_metric_mds_out", width = paste0(w, "px"), height = paste0(h, "px"))
  })

  metric_mds_reactive <- reactive({
    req(rv$benchmark_df)
    cat_sel <- if (identical(input$sel_mds_cat, "all")) NULL else input$sel_mds_cat
    plot_metric_mds(
      benchmark_data = rv$benchmark_df,
      category = cat_sel,
      ordination_by = input$sel_mds_by,
      base_size = 13
    )
  })
  output$plot_metric_mds_out <- renderPlot({ metric_mds_reactive() })
  output$download_mds_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".jpeg") },
    content = function(file) {
      w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width / 80 else 12
      h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height / 80 else 7.5
      export_single_jpeg(file, metric_mds_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_mds_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".pdf") },
    content = function(file) {
      w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width / 80 else 12
      h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height / 80 else 7.5
      grDevices::pdf(file, width = w, height = h)
      print(metric_mds_reactive())
      grDevices::dev.off()
    }
  )
  
  # ----------------------------------------------------------------------------
  # 8. Cell Embeddings (UMAP, t-SNE, PCA) & Quantitative Quality Metrics
  # ----------------------------------------------------------------------------
  output$ui_emb_sim_picker <- renderUI({
    sim_names <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      names(rv$sim_matrices)
    } else if (!is.null(rv$methods)) {
      rv$methods
    } else {
      "Simulated"
    }
    
    tagList(
      div(
        class = "d-flex justify-content-between align-items-center mb-1",
        tags$b("Choose Simulators to Include:", style = "font-size: 0.85rem; color: #1B4F72;"),
        div(
          actionLink("link_emb_select_all", "All", style = "font-size: 0.75rem; margin-right: 6px;"),
          actionLink("link_emb_clear_all", "Clear", style = "font-size: 0.75rem;")
        )
      ),
      checkboxGroupInput(
        "chk_emb_methods", label = NULL,
        choices = sim_names,
        selected = sim_names
      )
    )
  })
  
  observeEvent(input$link_emb_select_all, {
    sim_names <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) names(rv$sim_matrices) else rv$methods
    updateCheckboxGroupInput(session, "chk_emb_methods", selected = sim_names)
  })
  observeEvent(input$link_emb_clear_all, {
    updateCheckboxGroupInput(session, "chk_emb_methods", selected = character(0))
  })
  
  observe({
    sim_names <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) names(rv$sim_matrices) else rv$methods
    if (!is.null(sim_names) && length(sim_names) > 0) {
      updateSelectInput(session, "sel_emb_single_sim", choices = sim_names, selected = sim_names[1])
    }
  })
  
  # Reactive cell embeddings calculation
  reactive_cell_embeddings <- reactive({
    req(rv$toy_ref)
    ref_mat <- rv$toy_ref
    sim_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      rv$sim_matrices
    } else if (!is.null(rv$toy_sim)) {
      list("Simulated" = rv$toy_sim)
    } else {
      NULL
    }
    req(sim_list)
    
    red <- if (!is.null(input$sel_emb_reduction)) input$sel_emb_reduction else "umap"
    n_pcs <- if (!is.null(input$sld_emb_pcs)) input$sld_emb_pcs else 20
    perp <- if (!is.null(input$sld_emb_perp)) input$sld_emb_perp else 15
    n_neigh <- if (!is.null(input$sld_emb_neighbors)) input$sld_emb_neighbors else 15
    
    compute_dataset_embeddings(
      reference = ref_mat,
      simulated = sim_list,
      reduction = red,
      n_pcs = n_pcs,
      perplexity = perp,
      n_neighbors = n_neigh,
      cell_types = rv$cell_types,
      batch = rv$batch
    )
  })
  
  # Reactive embedding plot
  reactive_emb_plot <- reactive({
    emb_data <- reactive_cell_embeddings()
    req(emb_data, nrow(emb_data) > 0)
    
    red <- if (!is.null(input$sel_emb_reduction)) input$sel_emb_reduction else "umap"
    lay <- if (!is.null(input$sel_emb_layout)) input$sel_emb_layout else "facet"
    col_by <- if (!is.null(input$sel_emb_color)) input$sel_emb_color else "cell_type"
    pt_size <- if (!is.null(input$sld_emb_pt_size)) input$sld_emb_pt_size else 1.0
    alpha <- if (!is.null(input$sld_emb_alpha)) input$sld_emb_alpha else 0.8
    
    chosen_sims <- if (lay == "side_by_side") {
      input$sel_emb_single_sim
    } else {
      input$chk_emb_methods
    }
    
    plot_dataset_embeddings(
      embedding_data = emb_data,
      reduction = red,
      layout = lay,
      color_by = col_by,
      selected_methods = chosen_sims,
      pt_size = pt_size,
      alpha = alpha,
      base_size = 12
    )
  })
  
  output$ui_plot_cell_embeddings <- renderUI({
    lay <- if (!is.null(input$sel_emb_layout)) input$sel_emb_layout else "facet"
    n_methods <- if (lay == "side_by_side") 2 else max(2, length(input$chk_emb_methods) + 1)
    
    h_px <- if (lay == "side_by_side" || n_methods <= 3) 520 else if (n_methods <= 6) 720 else 960
    w_px <- if (lay == "side_by_side") 960 else 1150
    
    plotOutput("plot_cell_embeddings_out", width = paste0(w_px, "px"), height = paste0(h_px, "px"))
  })
  
  output$plot_cell_embeddings_out <- renderPlot({
    reactive_emb_plot()
  })
  
  output$table_emb_quality_metrics <- renderDT({
    emb_data <- reactive_cell_embeddings()
    req(emb_data, nrow(emb_data) > 0)
    
    q_df <- compute_embedding_quality_metrics(emb_data)
    
    datatable(
      q_df,
      rownames = FALSE,
      options = list(
        dom = "t",
        pageLength = 20,
        scrollX = TRUE
      ),
      class = "compact stripe hover"
    ) %>%
      formatStyle(
        "Role",
        backgroundColor = styleEqual(c("Reference", "Simulated"), c("#EBF5FB", "#FEF9E7")),
        fontWeight = "bold"
      )
  })
  
  output$download_emb_jpeg <- downloadHandler(
    filename = function() {
      paste0("scSimEval_cell_embeddings_", if (!is.null(input$sel_emb_reduction)) input$sel_emb_reduction else "umap", "_", Sys.Date(), ".jpeg")
    },
    content = function(file) {
      export_single_jpeg(file, reactive_emb_plot(), width = 14, height = 9, dpi = 600)
    }
  )
  
  output$download_emb_pdf <- downloadHandler(
    filename = function() {
      paste0("scSimEval_cell_embeddings_", if (!is.null(input$sel_emb_reduction)) input$sel_emb_reduction else "umap", "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      grDevices::pdf(file, width = 14, height = 9)
      print(reactive_emb_plot())
      grDevices::dev.off()
    }
  )
  
  # ----------------------------------------------------------------------------
  # Tab 5: Download Results (Excel, CSV, RDS, PDF, and Complete ZIP)
  # ----------------------------------------------------------------------------
  
  output$download_excel <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".xlsx") },
    content = function(file) {
      req(rv$benchmark_df)
      export_excel_workbook(file, rv$benchmark_df, leaderboard_reactive(), rv$dataset_summary_df)
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.csv(rv$benchmark_df, file, row.names = FALSE)
    }
  )
  
  output$download_txt <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".txt") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.table(rv$benchmark_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  output$download_rds <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".rds") },
    content = function(file) {
      req(rv$benchmark_df)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leaderboard_reactive(),
        dataset_summary = rv$dataset_summary_df
      ), file)
    }
  )
  
  output$download_all_plots_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_all_plots_report_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(rv$benchmark_df)
      generate_all_plots_pdf(file, rv$benchmark_df, rv$toy_ref, rv$toy_sim, rv$sim_matrices, rv$cell_types, rv$batch)
    }
  )
  
  output$download_complete_zip <- downloadHandler(
    filename = function() { paste0("scSimEval_complete_benchmark_results_", Sys.Date(), ".zip") },
    content = function(file) {
      req(rv$benchmark_df)
      
      tmp_dir <- file.path(tempdir(), paste0("scSimEval_bundle_", as.integer(Sys.time())))
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Structured subdirectories
      metrics_dir <- file.path(tmp_dir, "metrics")
      fig_dir <- file.path(tmp_dir, "figures")
      fig_grouped_dir <- file.path(fig_dir, "grouped")
      fig_indiv_dir <- file.path(fig_dir, "individual")
      dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(fig_grouped_dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(fig_indiv_dir, showWarnings = FALSE, recursive = TRUE)
      
      leader_df <- leaderboard_reactive()
      
      # 1. Excel Workbook (.xlsx)
      export_excel_workbook(file.path(tmp_dir, "scSimEval_benchmark_results.xlsx"), rv$benchmark_df, leader_df, rv$dataset_summary_df)
      try(file.copy(file.path(tmp_dir, "scSimEval_benchmark_results.xlsx"), file.path(metrics_dir, "all_benchmark_metrics.xlsx")), silent = TRUE)
      
      # 2. Master CSV Table (.csv)
      utils::write.csv(rv$benchmark_df, file.path(tmp_dir, "scSimEval_benchmark_results.csv"), row.names = FALSE)
      utils::write.csv(rv$benchmark_df, file.path(metrics_dir, "all_benchmark_metrics.csv"), row.names = FALSE)
      
      # 2b. Master TXT Table (.txt - tab delimited)
      utils::write.table(rv$benchmark_df, file.path(tmp_dir, "scSimEval_benchmark_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
      utils::write.table(rv$benchmark_df, file.path(metrics_dir, "all_benchmark_metrics.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      # 2c. Method Rankings Leaderboard (.csv & .txt)
      if (!is.null(leader_df)) {
        utils::write.csv(leader_df, file.path(tmp_dir, "scSimEval_method_rankings.csv"), row.names = FALSE)
        utils::write.csv(leader_df, file.path(metrics_dir, "method_rankings_leaderboard.csv"), row.names = FALSE)
        utils::write.table(leader_df, file.path(tmp_dir, "scSimEval_method_rankings.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        utils::write.table(leader_df, file.path(metrics_dir, "method_rankings_leaderboard.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
      }
      
      # 2d. Dataset Properties Summary (.csv & .txt)
      if (!is.null(rv$dataset_summary_df)) {
        utils::write.csv(rv$dataset_summary_df, file.path(tmp_dir, "scSimEval_dataset_properties_summary.csv"), row.names = FALSE)
        utils::write.csv(rv$dataset_summary_df, file.path(metrics_dir, "dataset_properties_summary.csv"), row.names = FALSE)
        utils::write.table(rv$dataset_summary_df, file.path(tmp_dir, "scSimEval_dataset_properties_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        utils::write.table(rv$dataset_summary_df, file.path(metrics_dir, "dataset_properties_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
      }
      
      # 3. RDS Object (.rds)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leader_df,
        dataset_summary = rv$dataset_summary_df
      ), file.path(tmp_dir, "scSimEval_benchmark_results.rds"))
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leader_df,
        dataset_summary = rv$dataset_summary_df
      ), file.path(metrics_dir, "all_benchmark_metrics.rds"))
      
      # 4. Multi-Page PDF Report (Compiled)
      generate_all_plots_pdf(file.path(tmp_dir, "scSimEval_all_plots_report.pdf"), rv$benchmark_df, rv$toy_ref, rv$toy_sim, rv$sim_matrices, rv$cell_types, rv$batch)
      try(file.copy(file.path(tmp_dir, "scSimEval_all_plots_report.pdf"), file.path(fig_grouped_dir, "all_figures_compiled_report.pdf")), silent = TRUE)
      
      # 5. Grouped Figures (Both 600 DPI Publication JPEG & Vectorized PDF)
      p_bubble <- bubble_plot_reactive()
      p_sum <- eval_summary_reactive()
      p_scale <- scale_bench_reactive()
      p_box <- metric_box_reactive()
      p_heat <- metric_heat_reactive()
      p_pca <- metric_pca_reactive()
      p_mds <- metric_mds_reactive()
      
      try({
        export_single_jpeg(file.path(fig_dir, "01_bubble_matrix.jpeg"), p_bubble, width = 23, height = 7, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "01_bubble_matrix_comparative.jpeg"), p_bubble, width = 23, height = 7, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "01_bubble_matrix_comparative.pdf"), p_bubble, width = 23, height = 7)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "02_evaluation_summary.jpeg"), p_sum, width = 13, height = 7.5, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "02_evaluation_summary_grouped.jpeg"), p_sum, width = 13, height = 7.5, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "02_evaluation_summary_grouped.pdf"), p_sum, width = 13, height = 7.5)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "03_scalability_benchmark.jpeg"), p_scale, width = 13, height = 8, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "03_scalability_benchmark_grouped.jpeg"), p_scale, width = 13, height = 8, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "03_scalability_benchmark_grouped.pdf"), p_scale, width = 13, height = 8)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "04_metric_boxplots.jpeg"), p_box, width = 13, height = 7.5, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "04_metric_boxplots_by_category.jpeg"), p_box, width = 13, height = 7.5, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "04_metric_boxplots_by_category.pdf"), p_box, width = 13, height = 7.5)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "05_metric_heatmap.jpeg"), p_heat, width = 14, height = 12, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "05_metric_performance_heatmap.jpeg"), p_heat, width = 14, height = 12, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "05_metric_performance_heatmap.pdf"), p_heat, width = 14, height = 12)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "06_metric_pca.jpeg"), p_pca, width = 13, height = 7.5, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "06_simulator_pca_ordination.jpeg"), p_pca, width = 13, height = 7.5, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "06_simulator_pca_ordination.pdf"), p_pca, width = 13, height = 7.5)
      }, silent = TRUE)
      
      try({
        export_single_jpeg(file.path(fig_dir, "07_metric_mds.jpeg"), p_mds, width = 13, height = 7.5, dpi = 600)
        export_single_jpeg(file.path(fig_grouped_dir, "07_simulator_mds_ordination.jpeg"), p_mds, width = 13, height = 7.5, dpi = 600)
        export_single_pdf(file.path(fig_grouped_dir, "07_simulator_mds_ordination.pdf"), p_mds, width = 13, height = 7.5)
      }, silent = TRUE)
      
      if (!is.null(rv$toy_ref) && !is.null(rv$toy_sim)) {
        try({
          p_dist <- dist_qc_reactive()
          export_single_jpeg(file.path(fig_dir, "08_distribution_qc.jpeg"), p_dist, width = 14, height = 9, dpi = 600)
          export_single_jpeg(file.path(fig_grouped_dir, "08_comparative_distribution_qc.jpeg"), p_dist, width = 14, height = 9, dpi = 600)
          export_single_pdf(file.path(fig_grouped_dir, "08_comparative_distribution_qc.pdf"), p_dist, width = 14, height = 9)
        }, silent = TRUE)
      }
      
      # Cell Embeddings Comparison Grids (UMAP, t-SNE, PCA)
      if (!is.null(rv$toy_ref) && (!is.null(rv$sim_matrices) || !is.null(rv$toy_sim))) {
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        for (red in c("umap", "tsne", "pca")) {
          try({
            emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = red, cell_types = rv$cell_types, batch = rv$batch)
            p_red <- plot_dataset_embeddings(emb_obj, reduction = red, layout = "grid", base_size = 11)
            prefix <- if (red == "umap") "09_cell_embeddings_umap_grid" else if (red == "tsne") "10_cell_embeddings_tsne_grid" else "11_cell_embeddings_pca_grid"
            export_single_jpeg(file.path(fig_dir, paste0(prefix, ".jpeg")), p_red, width = 14, height = 9, dpi = 600)
            export_single_jpeg(file.path(fig_grouped_dir, paste0(prefix, ".jpeg")), p_red, width = 14, height = 9, dpi = 600)
            export_single_pdf(file.path(fig_grouped_dir, paste0(prefix, ".pdf")), p_red, width = 14, height = 9)
            
            # Individual simulator 1-to-1 comparison embeddings
            for (sim_nm in names(s_list)) {
              p_sim <- plot_dataset_embeddings(emb_obj, reduction = red, layout = "compare", compare_sim = sim_nm, base_size = 11)
              safe_sim <- gsub("[^A-Za-z0-9_-]", "_", sim_nm)
              export_single_jpeg(file.path(fig_indiv_dir, paste0("cell_embeddings_compare_", red, "_", safe_sim, ".jpeg")), p_sim, width = 11, height = 5.5, dpi = 300)
              export_single_pdf(file.path(fig_indiv_dir, paste0("cell_embeddings_compare_", red, "_", safe_sim, ".pdf")), p_sim, width = 11, height = 5.5)
            }
          }, silent = TRUE)
        }
      }
      
      # 6. Individual Figures (Individual Metric Barplots)
      if ("Metric" %in% colnames(rv$benchmark_df)) {
        u_metrics <- unique(as.character(rv$benchmark_df$Metric))
        for (m in u_metrics) {
          try({
            p_m <- plot_individual_metric_bar(rv$benchmark_df, metric = m, score_type = "normalized", base_size = 11)
            safe_m <- gsub("[^A-Za-z0-9_-]", "_", m)
            export_single_jpeg(file.path(fig_indiv_dir, paste0("metric_bar_", safe_m, ".jpeg")), p_m, width = 10, height = 6.5, dpi = 300)
            export_single_pdf(file.path(fig_indiv_dir, paste0("metric_bar_", safe_m, ".pdf")), p_m, width = 10, height = 6.5)
          }, silent = TRUE)
        }
      }
      
      zip_files <- list.files(tmp_dir, full.names = FALSE, recursive = TRUE)
      zip::zip(file, files = zip_files, root = tmp_dir)
      unlink(tmp_dir, recursive = TRUE)
    }
  )

  # Dynamic update for simulator choice in individual figure export
  observe({
    req(rv$benchmark_df)
    sims <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      names(rv$sim_matrices)
    } else if ("Method" %in% colnames(rv$benchmark_df)) {
      unique(as.character(rv$benchmark_df$Method))
    } else {
      character(0)
    }
    updateSelectInput(session, "sel_export_compare_sim", choices = sims, selected = sims[1])
  })

  # Reactive plot generator for Tab 5 Individual Figure Export
  selected_export_plot_reactive <- reactive({
    req(rv$benchmark_df, input$sel_export_figure_type)
    type <- input$sel_export_figure_type
    
    switch(
      type,
      "bubble" = bubble_plot_reactive(),
      "summary" = eval_summary_reactive(),
      "scalability" = scale_bench_reactive(),
      "boxplots" = metric_box_reactive(),
      "heatmap" = metric_heat_reactive(),
      "pca_metric" = metric_pca_reactive(),
      "mds_metric" = metric_mds_reactive(),
      "dist_qc" = {
        req(rv$toy_ref, rv$toy_sim)
        dist_qc_reactive()
      },
      "emb_umap" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "umap", cell_types = rv$cell_types, batch = rv$batch)
        plot_dataset_embeddings(emb_obj, reduction = "umap", layout = "grid", base_size = 11)
      },
      "emb_tsne" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "tsne", cell_types = rv$cell_types, batch = rv$batch)
        plot_dataset_embeddings(emb_obj, reduction = "tsne", layout = "grid", base_size = 11)
      },
      "emb_pca" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "pca", cell_types = rv$cell_types, batch = rv$batch)
        plot_dataset_embeddings(emb_obj, reduction = "pca", layout = "grid", base_size = 11)
      },
      "emb_compare_umap" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "umap", cell_types = rv$cell_types, batch = rv$batch)
        sim_choice <- input$sel_export_compare_sim
        if (is.null(sim_choice) || !sim_choice %in% names(s_list)) sim_choice <- names(s_list)[1]
        plot_dataset_embeddings(emb_obj, reduction = "umap", layout = "compare", compare_sim = sim_choice, base_size = 11)
      },
      "emb_compare_tsne" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "tsne", cell_types = rv$cell_types, batch = rv$batch)
        sim_choice <- input$sel_export_compare_sim
        if (is.null(sim_choice) || !sim_choice %in% names(s_list)) sim_choice <- names(s_list)[1]
        plot_dataset_embeddings(emb_obj, reduction = "tsne", layout = "compare", compare_sim = sim_choice, base_size = 11)
      },
      "emb_compare_pca" = {
        req(rv$toy_ref)
        s_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) rv$sim_matrices else list("Simulated" = rv$toy_sim)
        emb_obj <- compute_dataset_embeddings(rv$toy_ref, s_list, reduction = "pca", cell_types = rv$cell_types, batch = rv$batch)
        sim_choice <- input$sel_export_compare_sim
        if (is.null(sim_choice) || !sim_choice %in% names(s_list)) sim_choice <- names(s_list)[1]
        plot_dataset_embeddings(emb_obj, reduction = "pca", layout = "compare", compare_sim = sim_choice, base_size = 11)
      },
      NULL
    )
  })

  # Individual figure PDF download
  output$download_selected_plot_pdf <- downloadHandler(
    filename = function() {
      paste0("scSimEval_", input$sel_export_figure_type, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- selected_export_plot_reactive()
      req(p)
      type <- input$sel_export_figure_type
      w <- if (type == "bubble") 23 else if (type == "heatmap") 14 else if (grepl("compare", type)) 11 else if (grepl("emb", type)) 14 else 13
      h <- if (type == "bubble") 7 else if (type == "heatmap") 12 else if (grepl("compare", type)) 5.5 else if (grepl("emb", type)) 9 else 7.5
      export_single_pdf(file, p, width = w, height = h)
    }
  )

  # Individual figure JPEG (600 DPI) download
  output$download_selected_plot_jpeg <- downloadHandler(
    filename = function() {
      paste0("scSimEval_", input$sel_export_figure_type, "_", Sys.Date(), ".jpeg")
    },
    content = function(file) {
      p <- selected_export_plot_reactive()
      req(p)
      type <- input$sel_export_figure_type
      w <- if (type == "bubble") 23 else if (type == "heatmap") 14 else if (grepl("compare", type)) 11 else if (grepl("emb", type)) 14 else 13
      h <- if (type == "bubble") 7 else if (type == "heatmap") 12 else if (grepl("compare", type)) 5.5 else if (grepl("emb", type)) 9 else 7.5
      export_single_jpeg(file, p, width = w, height = h, dpi = 600)
    }
  )
  
  # Searchable Master Table with Factor Filter Dropdowns
  output$table_master_export <- renderDT({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    # Convert text columns to factors for automatic dropdown filtering
    if ("Method" %in% colnames(df)) df$Method <- as.factor(df$Method)
    if ("Category" %in% colnames(df)) df$Category <- as.factor(df$Category)
    if ("Metric" %in% colnames(df)) df$Metric <- as.factor(df$Metric)
    
    datatable(
      df,
      filter = list(position = "top", clear = FALSE),
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE,
        searchHighlight = TRUE
      ),
      rownames = FALSE,
      class = "compact stripe hover"
    ) %>% formatRound(columns = which(sapply(df, is.numeric)), digits = 4)
  })
}

# ==============================================================================
# Run Shiny App
# ==============================================================================
shinyApp(ui = ui, server = server)
