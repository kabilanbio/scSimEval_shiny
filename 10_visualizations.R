# ============================================================================
# GLOBAL VARIABLE DECLARATIONS (for R CMD check NOTE suppression)
# ============================================================================
utils::globalVariables(c(
  ".data", "Dataset", "Data_Name", "Modality", "Category", "Property",
  "Metric", "Value", "Property_Clean", "Method", "Score", "Score_Norm",
  "Method_Class", "Is_Top", "Display_Metric", "Display_Category", "Score_Raw",
  "Feature_Pair", "Correlation", "Simulator", "Ref_Quantile", "Sim_Quantile",
  "ZeroFrac", "Mean", "label_txt", "x_pos", "y_pos",
  "MDS1", "MDS2", "PC1", "PC2", "PC1_scaled", "PC2_scaled", "mag",
  "Dim1", "Dim2", "Facet_Label",
  "Time_Sec", "Time_Type", "resource_cost", "throughput_cells_sec",
  "cpu_efficiency", "peak_memory_mb", "elapsed_time_seconds", "cpu_time_seconds",
  "Score_Label", "Panel", "Label", "ymin", "ymax", "is_even", "Is_Ref"
))

# ============================================================================
# CANONICAL METRIC -> CATEGORY MAPPING TABLES
# These tables form the scientific backbone of the visualization suite.
# They map internal metric/category strings emitted by the evaluation engine
# to the canonical eight-group evaluation framework.
# ============================================================================

#' @keywords internal
#' @keywords internal
.METRIC_DISPLAY_LABELS <- c(
  # (I) Distributional Properties (14)
  "KS"                               = "KS Distance",
  "Wasserstein"                      = "Wasserstein Dist.",
  "MAD"                              = "MAD",
  "MAE"                              = "MAE",
  "RMSE"                             = "RMSE",
  "OV"                               = "Overlap Coeff.",
  "Bhattacharyya"                    = "Bhattacharyya Dist.",
  "ECDF_DiffArea"                    = "ECDF Diff. Area",
  "Runs_Statistic"                   = "Runs Statistic",
  "Fasano_Franceschini_2D_KS"        = "Fasano-Franceschini 2D",
  "Peacock_2D_KS"                    = "Peacock 2D",
  "KDE_Bivariate_zstat"              = "2D KDE z-stat",
  "MMD"                              = "MMD",
  "Frechet_SingleCell_Distance"      = "Frechet SC Dist.",

  # (II) Correlations & Zero-Inflation (6)
  "bcv_discrepancy"                  = "Biological CV (BCV)",
  "dropout_curve_midpoint_x0"        = "Dropout Midpoint (x0)",
  "dropout_curve_slope_k"            = "Dropout Slope (k)",
  "mean_variance_r2"                 = "Mean-Variance R2",
  "cell_pearson_cor"                 = "Cell Pearson Corr.",
  "feature_pearson_cor"              = "Feature Pearson Corr.",

  # (III) Cellular Structure & Concordance (10)
  "silhouette_sim"                   = "Silhouette ASW",
  "dunn_sim"                         = "Dunn Index",
  "davies_bouldin_sim"               = "Davies-Bouldin",
  "calinski_harabasz_sim"            = "Calinski-Harabasz",
  "ari"                              = "ARI",
  "nmi"                              = "NMI",
  "ami"                              = "AMI",
  "v_measure"                        = "V-Measure",
  "neighborhood_purity"              = "Neighborhood Purity",
  "Generative_Precision"             = "Generative Precision",

  # (IV) Batch Effects & Confounder Mixing (7)
  "batch_silhouette"                 = "Batch Silhouette",
  "shannon_entropy"                  = "Shannon Entropy",
  "pcr_r2"                           = "PC Regression R2",
  "cms"                              = "CMS Score",
  "isi"                              = "LISI / ISI",
  "seurat_mixing_metric"             = "Seurat Mixing",
  "ldf_diff"                         = "ldfDiff",

  # (V) Biological Signal & Downstream Fidelity (15)
  "SimBench_SMAPE"                   = "SimBench SMAPE",
  "SimBench_DE_Fidelity_Score"       = "DE Fidelity (1-SMAPE)",
  "Log2FC_Pearson_Corr"              = "Log2FC Pearson Corr",
  "Log2FC_Spearman_Corr"             = "Log2FC Spearman Corr",
  "Top_DEG_Jaccard_Overlap"          = "Top DEG Jaccard",
  "DEG_Ratio"                        = "DEG Ratio (Sim/Real)",
  "PValue_Uniformity_Chisq"          = "P-Val Unif. (Chisq)",
  "Distribution_Score"               = "Distribution Score",
  "Classifier_Accuracy"              = "Classifier Accuracy",
  "Classifier_Macro_F1"              = "Classifier Macro-F1",
  "Classifier_Macro_Recall"          = "Classifier Macro-Recall",
  "Silhouette_Sim"                   = "Group Silhouette",
  "Silhouette_Discrepancy"           = "Silhouette Discrepancy",
  "PVE_Group_Sim"                    = "Group PVE",
  "PVE_Discrepancy"                  = "Group PVE Discrepancy",

  # (VI) Trajectory & Lineage Dynamics (2)
  "pseudotime_correlation"           = "Pseudotime Corr. (rho)",
  "tree_height_rmse"                 = "Lineage Tree RMSE",

  # (VII) Cross-Modal Coupling & Modularity (6)
  "cross_modal_accuracy"             = "Cross-Modal Accuracy",
  "cross_modal_F1"                   = "Cross-Modal F1",
  "FOSCTTM"                          = "FOSCTTM",
  "match_at_1"                       = "Match@1",
  "rv_coefficient"                   = "RV Coefficient",
  "module_correlation_r"             = "Module Corr. (r)",

  # (VIII) Computational Scalability (2)
  "elapsed_time_seconds"             = "Elapsed Time (s)",
  "peak_memory_mb"                   = "Peak RAM (MiB)"
)

#' @keywords internal
.METRIC_CATEGORY_MAP <- c(
  # (I) Distributional Properties (14)
  "KS"                               = "(I) Distributional Properties",
  "Wasserstein"                      = "(I) Distributional Properties",
  "MAD"                              = "(I) Distributional Properties",
  "MAE"                              = "(I) Distributional Properties",
  "RMSE"                             = "(I) Distributional Properties",
  "OV"                               = "(I) Distributional Properties",
  "Bhattacharyya"                    = "(I) Distributional Properties",
  "ECDF_DiffArea"                    = "(I) Distributional Properties",
  "Runs_Statistic"                   = "(I) Distributional Properties",
  "Fasano_Franceschini_2D_KS"        = "(I) Distributional Properties",
  "Peacock_2D_KS"                    = "(I) Distributional Properties",
  "KDE_Bivariate_zstat"              = "(I) Distributional Properties",
  "MMD"                              = "(I) Distributional Properties",
  "Frechet_SingleCell_Distance"      = "(I) Distributional Properties",

  # (II) Correlations & Zero-Inflation (6)
  "bcv_discrepancy"                  = "(II) Correlations & Zero-Inflation",
  "dropout_curve_midpoint_x0"        = "(II) Correlations & Zero-Inflation",
  "dropout_curve_slope_k"            = "(II) Correlations & Zero-Inflation",
  "mean_variance_r2"                 = "(II) Correlations & Zero-Inflation",
  "cell_pearson_cor"                 = "(II) Correlations & Zero-Inflation",
  "feature_pearson_cor"              = "(II) Correlations & Zero-Inflation",

  # (III) Cellular Structure & Concordance (10)
  "silhouette_sim"                   = "(III) Cellular Structure & Concordance",
  "dunn_sim"                         = "(III) Cellular Structure & Concordance",
  "davies_bouldin_sim"               = "(III) Cellular Structure & Concordance",
  "calinski_harabasz_sim"            = "(III) Cellular Structure & Concordance",
  "ari"                              = "(III) Cellular Structure & Concordance",
  "nmi"                              = "(III) Cellular Structure & Concordance",
  "ami"                              = "(III) Cellular Structure & Concordance",
  "v_measure"                        = "(III) Cellular Structure & Concordance",
  "neighborhood_purity"              = "(III) Cellular Structure & Concordance",
  "Generative_Precision"             = "(III) Cellular Structure & Concordance",

  # (IV) Batch Effects & Confounder Mixing (7)
  "batch_silhouette"                 = "(IV) Batch Effects & Confounder Mixing",
  "shannon_entropy"                  = "(IV) Batch Effects & Confounder Mixing",
  "pcr_r2"                           = "(IV) Batch Effects & Confounder Mixing",
  "cms"                              = "(IV) Batch Effects & Confounder Mixing",
  "isi"                              = "(IV) Batch Effects & Confounder Mixing",
  "seurat_mixing_metric"             = "(IV) Batch Effects & Confounder Mixing",
  "ldf_diff"                         = "(IV) Batch Effects & Confounder Mixing",

  # (V) Biological Signal & Downstream Fidelity (15)
  "SimBench_SMAPE"                   = "(V) Biological Signal & Downstream Fidelity",
  "SimBench_DE_Fidelity_Score"       = "(V) Biological Signal & Downstream Fidelity",
  "Log2FC_Pearson_Corr"              = "(V) Biological Signal & Downstream Fidelity",
  "Log2FC_Spearman_Corr"             = "(V) Biological Signal & Downstream Fidelity",
  "Top_DEG_Jaccard_Overlap"          = "(V) Biological Signal & Downstream Fidelity",
  "DEG_Ratio"                        = "(V) Biological Signal & Downstream Fidelity",
  "PValue_Uniformity_Chisq"          = "(V) Biological Signal & Downstream Fidelity",
  "Distribution_Score"               = "(V) Biological Signal & Downstream Fidelity",
  "Classifier_Accuracy"              = "(V) Biological Signal & Downstream Fidelity",
  "Classifier_Macro_F1"              = "(V) Biological Signal & Downstream Fidelity",
  "Classifier_Macro_Recall"          = "(V) Biological Signal & Downstream Fidelity",
  "Silhouette_Sim"                   = "(V) Biological Signal & Downstream Fidelity",
  "Silhouette_Discrepancy"           = "(V) Biological Signal & Downstream Fidelity",
  "PVE_Group_Sim"                    = "(V) Biological Signal & Downstream Fidelity",
  "PVE_Discrepancy"                  = "(V) Biological Signal & Downstream Fidelity",

  # (VI) Trajectory & Lineage Dynamics (2)
  "pseudotime_correlation"           = "(VI) Trajectory & Lineage Dynamics",
  "tree_height_rmse"                 = "(VI) Trajectory & Lineage Dynamics",

  # (VII) Cross-Modal Coupling & Modularity (6)
  "cross_modal_accuracy"             = "(VII) Cross-Modal Coupling & Modularity",
  "cross_modal_F1"                   = "(VII) Cross-Modal Coupling & Modularity",
  "FOSCTTM"                          = "(VII) Cross-Modal Coupling & Modularity",
  "match_at_1"                       = "(VII) Cross-Modal Coupling & Modularity",
  "rv_coefficient"                   = "(VII) Cross-Modal Coupling & Modularity",
  "module_correlation_r"             = "(VII) Cross-Modal Coupling & Modularity",

  # (VIII) Computational Scalability (2)
  "elapsed_time_seconds"             = "(VIII) Computational Scalability",
  "peak_memory_mb"                   = "(VIII) Computational Scalability"
)

.CANONICAL_CATEGORY_ORDER <- c(
  "(I) Distributional Properties",
  "(II) Correlations & Zero-Inflation",
  "(III) Cellular Structure & Concordance",
  "(IV) Batch Effects & Confounder Mixing",
  "(V) Biological Signal & Downstream Fidelity",
  "(VI) Trajectory & Lineage Dynamics",
  "(VII) Cross-Modal Coupling & Modularity",
  "(VIII) Computational Scalability"
)

#' @keywords internal
.CATEGORY_TWO_LINE_MAP <- c(
  "(I) Distributional Properties"              = "(I) Distributional\nProperties",
  "(II) Correlations & Zero-Inflation"         = "(II) Correlations &\nZero-Inflation",
  "(III) Cellular Structure & Concordance"     = "(III) Cellular Structure\n& Concordance",
  "(IV) Batch Effects & Confounder Mixing"     = "(IV) Batch Effects &\nConfounder Mixing",
  "(V) Biological Signal & Downstream Fidelity"= "(V) Biological Signal &\nDownstream Fidelity",
  "(VI) Trajectory & Lineage Dynamics"         = "(VI) Trajectory &\nLineage Dynamics",
  "(VII) Cross-Modal Coupling & Modularity"    = "(VII) Cross-Modal\nCoupling & Modularity",
  "(VIII) Computational Scalability"           = "(VIII) Computational\nScalability"
)

#' @keywords internal
.COMPACT_CATEGORY_STRIPS <- c(
  "(I) Distributional Properties"              = "(I)\nDistribution",
  "(II) Correlations & Zero-Inflation"         = "(II)\nCorrelation",
  "(III) Cellular Structure & Concordance"     = "(III)\nCell Structure",
  "(IV) Batch Effects & Confounder Mixing"     = "(IV)\nBatch Mixing",
  "(V) Biological Signal & Downstream Fidelity"= "(V)\nBio-Signal & DE",
  "(VI) Trajectory & Lineage Dynamics"         = "(VI)\nTraj.",
  "(VII) Cross-Modal Coupling & Modularity"    = "(VII)\nCross-Modal",
  "(VIII) Computational Scalability"           = "(VIII)\nScalability"
)

#' @keywords internal
.CANONICAL_CATEGORY_COLORS <- c(
  "(I) Distributional Properties"              = "#2E86AB",
  "(II) Correlations & Zero-Inflation"         = "#17A589",
  "(III) Cellular Structure & Concordance"     = "#C0392B",
  "(IV) Batch Effects & Confounder Mixing"     = "#D4AC0D",
  "(V) Biological Signal & Downstream Fidelity"= "#CA6F1E",
  "(VI) Trajectory & Lineage Dynamics"         = "#7D3C98",
  "(VII) Cross-Modal Coupling & Modularity"    = "#2E4057",
  "(VIII) Computational Scalability"           = "#1E8449"
)

#' @keywords internal
.SUMMARY_TYPE_MAP <- c(
  # Gene / Feature level (red)
  "feature_pearson_cor"         = "gene",
  "gene_cor"                    = "gene",
  "gene_detection_frequency"    = "gene",
  "zero_fraction_feature"       = "gene",
  "mean_expression"             = "gene",
  "average_logcpm"              = "gene",
  "variance"                    = "gene",
  "variance_logcpm"             = "gene",
  "cv"                          = "gene",
  "coefficient_of_variation"    = "gene",
  "dispersion"                  = "gene",
  "bcv_discrepancy"             = "gene",
  "mean_variance_r2"            = "gene",
  "dropout_curve_midpoint_x0"   = "gene",
  "dropout_curve_slope_k"       = "gene",
  "Log2FC_Pearson_Corr"         = "gene",
  "Log2FC_Spearman_Corr"        = "gene",
  "Top_DEG_Jaccard_Overlap"     = "gene",
  "DEG_Ratio"                   = "gene",
  "SimBench_DE_Fidelity_Score"  = "gene",
  "SimBench_SMAPE"              = "gene",
  "Distribution_Score"          = "gene",

  # Cell level (blue)
  "cell_pearson_cor"            = "cell",
  "cell_cor"                    = "cell",
  "cell_detection_frequency"    = "cell",
  "zero_fraction_cell"          = "cell",
  "library_size"                = "cell",
  "log_library_size"            = "cell",
  "pca_dist"                    = "cell",
  "cell_to_cell_distance"       = "cell",
  "local_density_factor"        = "cell",
  "ldf"                         = "cell",
  "ldf_diff"                    = "cell",
  "knn_hubness"                 = "cell",
  "neighborhood_purity"         = "cell",
  "KNN_occurences"              = "cell",
  "tmm_factor"                  = "cell",
  "effective_library_size"      = "cell",

  # Global / Structure / Batch / Trajectory / Multiomics (green)
  "silhouette_sim"              = "global",
  "Silhouette_Sim"              = "global",
  "silhouette_width"            = "global",
  "cms"                         = "global",
  "cell_specific_mixing_score"  = "global",
  "pcr_r2"                      = "global",
  "percent_variance_explained"  = "global",
  "batch_silhouette"            = "global",
  "shannon_entropy"             = "global",
  "isi"                         = "global",
  "seurat_mixing_metric"        = "global",
  "ari"                         = "global",
  "nmi"                         = "global",
  "ami"                         = "global",
  "v_measure"                   = "global",
  "dunn_sim"                    = "global",
  "davies_bouldin_sim"          = "global",
  "calinski_harabasz_sim"       = "global",
  "pseudotime_correlation"      = "global",
  "tree_height_rmse"            = "global",
  "FOSCTTM"                     = "global",
  "cross_modal_accuracy"        = "global",
  "cross_modal_F1"              = "global",
  "rv_coefficient"              = "global",
  "module_correlation_r"        = "global",
  "PVE_Group_Sim"               = "global",
  "PVE_Discrepancy"             = "global",
  "Classifier_Accuracy"         = "global",
  "elapsed_time_seconds"        = "global",
  "peak_memory_mb"              = "global"
)

#' @keywords internal
.SUMMARY_DISPLAY_NAME_MAP <- c(
  "feature_pearson_cor"         = "gene-to-gene correlation",
  "gene_cor"                    = "gene-to-gene correlation",
  "gene_detection_frequency"    = "gene detection frequency",
  "zero_fraction_feature"       = "gene detection frequency",
  "mean_expression"             = "average of logCPM",
  "average_logcpm"              = "average of logCPM",
  "variance"                    = "variance of logCPM",
  "variance_logcpm"             = "variance of logCPM",
  "cv"                          = "coefficient of variation",
  "coefficient_of_variation"    = "coefficient of variation",
  "dispersion"                  = "dispersion",
  "bcv_discrepancy"             = "BCV discrepancy",
  "mean_variance_r2"            = "mean-variance fit",
  "dropout_curve_midpoint_x0"   = "dropout midpoint",
  "dropout_curve_slope_k"       = "dropout rate",
  "Log2FC_Pearson_Corr"         = "log2FC Pearson cor",
  "Top_DEG_Jaccard_Overlap"     = "DEG Jaccard overlap",
  "Distribution_Score"          = "distribution score",

  "cell_pearson_cor"            = "cell-to-cell correlation",
  "cell_cor"                    = "cell-to-cell correlation",
  "cell_detection_frequency"    = "cell detection frequency",
  "zero_fraction_cell"          = "cell detection frequency",
  "library_size"                = "library size",
  "log_library_size"            = "log-library size",
  "pca_dist"                    = "cell-to-cell distance",
  "cell_to_cell_distance"       = "cell-to-cell distance",
  "local_density_factor"        = "local density factor",
  "ldf"                         = "local density factor",
  "ldf_diff"                    = "local density factor",
  "knn_hubness"                 = "KNN occurences",
  "neighborhood_purity"         = "KNN occurences",
  "KNN_occurences"              = "KNN occurences",

  "silhouette_sim"              = "silhouette width",
  "Silhouette_Sim"              = "silhouette width",
  "silhouette_width"            = "silhouette width",
  "cms"                         = "cell-specific mixing score",
  "cell_specific_mixing_score"  = "cell-specific mixing score",
  "pcr_r2"                      = "percent variance explained",
  "percent_variance_explained"  = "percent variance explained",
  "batch_silhouette"            = "batch silhouette",
  "shannon_entropy"             = "batch entropy",
  "seurat_mixing_metric"        = "mixing metric",
  "pseudotime_correlation"      = "pseudotime correlation",
  "tree_height_rmse"            = "tree height RMSE",
  "FOSCTTM"                     = "FOSCTTM cross-modality",
  "rv_coefficient"              = "RV coefficient",
  "elapsed_time_seconds"        = "elapsed time (s)",
  "peak_memory_mb"              = "peak RAM (MiB)"
)


#' @keywords internal
.LEGACY_CATEGORY_MAP <- c(
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
  "Cross-Modal Relationships"               = "(VII) Cross-Modal Coupling & Modularity",
  "Cross-Modal Coupling & Modularity"       = "(VII) Cross-Modal Coupling & Modularity",
  "Cross-Modal"                             = "(VII) Cross-Modal Coupling & Modularity",
  "Computational Scalability"               = "(VIII) Computational Scalability",
  "Scalability"                             = "(VIII) Computational Scalability"
)

#' @keywords internal
#' @keywords internal
.HIGHER_IS_BETTER_METRICS <- c(
  "OV",
  "mean_variance_r2", "cell_pearson_cor", "feature_pearson_cor",
  "silhouette_sim", "dunn_sim", "calinski_harabasz_sim", "ari", "nmi", "ami", "v_measure", "neighborhood_purity", "Generative_Precision",
  "cms", "isi", "seurat_mixing_metric", "shannon_entropy",
  "SimBench_DE_Fidelity_Score", "Log2FC_Pearson_Corr", "Log2FC_Spearman_Corr", "Top_DEG_Jaccard_Overlap",
  "Distribution_Score", "Classifier_Accuracy", "Classifier_Macro_F1", "Classifier_Macro_Recall",
  "Silhouette_Sim", "PVE_Group_Sim",
  "pseudotime_correlation",
  "cross_modal_accuracy", "cross_modal_F1", "match_at_1", "rv_coefficient", "module_correlation_r"
)

# Shared publication ggplot2 theme
#' @keywords internal
.pub_theme <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major   = ggplot2::element_line(color = "#EAECF0", linewidth = 0.35),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(color = "#BDC3C7", fill = NA, linewidth = 0.6),
      axis.title         = ggplot2::element_text(face = "bold", size = base_size * 0.95, color = "#2C3E50"),
      axis.text          = ggplot2::element_text(size = base_size * 0.88, color = "#2C3E50"),
      axis.ticks         = ggplot2::element_line(color = "#BDC3C7", linewidth = 0.4),
      strip.background   = ggplot2::element_rect(fill = "#F4F6F7", color = "#BDC3C7", linewidth = 0.5),
      strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.9, color = "#1A252F"),
      legend.background  = ggplot2::element_blank(),
      legend.key         = ggplot2::element_blank(),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.88),
      legend.text        = ggplot2::element_text(size = base_size * 0.82),
      legend.position    = "bottom",
      legend.box         = "vertical",
      plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.18,
                                                  color = "#1A252F", hjust = 0.5, margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size * 0.92, color = "#566573",
                                                  hjust = 0.5, margin = ggplot2::margin(b = 8)),
      plot.caption       = ggplot2::element_text(size = base_size * 0.75, color = "#85929E",
                                                  hjust = 0, margin = ggplot2::margin(t = 6)),
      plot.margin        = ggplot2::margin(t = 10, r = 12, b = 8, l = 10)
    )
}


# ============================================================================
# 1. plot_distribution_qc()
# ============================================================================

#' Plot Single-Cell Summary Distribution Quality (Comparative Overlays)
#'
#' Compares empirical reference and simulated datasets across canonical single-cell
#' summary properties: library size, TMM, effective library size, mean expression,
#' feature variance, scaled variance, fraction zero per cell, fraction zero per gene,
#' cell-cell correlation, gene-gene correlation, and bivariate relationships (mean vs. variance,
#' mean vs. zero fraction, library size vs. zero fraction, and biological signals).
#' Designed after landmark benchmarking studies (countsimQC / Duo et al.) with a clean,
#' human-crafted scientific layout featuring solid visible contour lines, semi-transparent
#' color-filled density curves, and simple understated titles.
#'
#' \if{html}{\figure{comparative_distribution_qc.png}{options: width="800" alt="Comparative Single-Cell Distribution QC"}}
#'
#' @param ref_data Numeric count matrix for empirical reference (features x cells).
#' @param sim_data Numeric count matrix for simulated data (features x cells), or a named list
#'   of simulated count matrices representing multiple simulators (e.g. \code{list("scDesign3" = m1, "Splatter" = m2)}).
#' @param layout Character. Visualization layout: \code{"comprehensive"} (default 14-panel benchmarking grid
#'   matching canonical single-cell papers) or \code{"density"} (faceted density curves only).
#' @param properties Optional character vector of specific properties to include if custom selection is desired.
#' @param title Character. Main plot title. Default \code{"Data properties"}. Set to \code{NULL} for no title.
#' @param ref_name Character. Label for the empirical reference in the legend. Default \code{"Original"}.
#' @param sim_name Character. Label for the simulation if a single matrix is provided. Default \code{"Simulation"}.
#' @param palette Optional named character vector of colors. Default uses warm red for reference and classic steel blue for simulation.
#' @param cell_types Optional factor or character vector of cell identities for computing biological signal proportions.
#' @param base_size Numeric. Base font size. Default \code{9}.
#'
#' @return A \code{ggplot} or \code{patchwork} composite object containing the comparative panels.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
#' @examples
#' data(example_scrna)
#' p <- plot_distribution_qc(example_scrna$ref, example_scrna$sim)
#' if (requireNamespace("ggplot2", quietly = TRUE)) print(p)
plot_distribution_qc <- function(
  ref_data,
  sim_data,
  layout        = c("comprehensive", "density"),
  properties    = NULL,
  title         = "Data properties",
  ref_name      = "Original",
  sim_name      = "Simulation",
  palette       = NULL,
  cell_types    = NULL,
  base_size     = 9
) {
  layout <- match.arg(layout)
  ref_mat <- as.matrix(ref_data)

  if (!is.list(sim_data) || is.data.frame(sim_data)) {
    sims <- list()
    sims[[sim_name]] <- as.matrix(sim_data)
  } else {
    sims <- lapply(sim_data, as.matrix)
  }

  sim_names <- names(sims)
  all_datasets <- c(sim_names, ref_name) # Simulators first, then Original

  if (is.null(palette)) {
    curated_sim_colors <- c("#2b5c8f", "#2e7d32", "#d95f02", "#7b1fa2", "#00838f", "#455a64")
    pal <- c()
    for (i in seq_along(sim_names)) {
      pal[sim_names[i]] <- curated_sim_colors[((i - 1) %% length(curated_sim_colors)) + 1]
    }
    pal[ref_name] <- "#b23939" # Warm brick red for Original reference
  } else {
    pal <- palette
  }

  theme_human <- function(bs = base_size) {
    ggplot2::theme_classic(base_size = bs) +
      ggplot2::theme(
        panel.border       = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line          = ggplot2::element_blank(),
        panel.grid.major   = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.ticks         = ggplot2::element_line(color = "black", linewidth = 0.4),
        axis.text          = ggplot2::element_text(color = "black", size = bs * 0.84),
        axis.title         = ggplot2::element_text(color = "black", size = bs * 0.88),
        plot.title         = ggplot2::element_text(color = "black", size = bs * 0.95, face = "plain", hjust = 0.5, margin = ggplot2::margin(b = 4)),
        legend.position    = "none",
        plot.margin        = ggplot2::margin(t = 4, r = 6, b = 4, l = 6)
      )
  }

  # 1. Library size
  lib_r <- colSums(ref_mat)
  df_lib_list <- list(data.frame(val = lib_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_lib_list[[sn]] <- data.frame(val = colSums(sims[[sn]]), Dataset = sn)
  }
  df_lib <- do.call(rbind, df_lib_list)
  df_lib$Dataset <- factor(df_lib$Dataset, levels = all_datasets)
  p_lib <- ggplot2::ggplot(df_lib, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Library size", x = "library size", y = "density") +
    theme_human()

  # 2. TMM Normalization Factor
  if (requireNamespace("edgeR", quietly = TRUE)) {
    tmm_r <- edgeR::normLibSizes(ref_mat)
    tmm_sims <- lapply(sims, edgeR::normLibSizes)
  } else {
    tmm_r <- colSums(ref_mat) / stats::median(colSums(ref_mat))
    tmm_sims <- lapply(sims, function(m) colSums(m) / stats::median(colSums(m)))
  }
  df_tmm_list <- list(data.frame(val = tmm_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_tmm_list[[sn]] <- data.frame(val = tmm_sims[[sn]], Dataset = sn)
  }
  df_tmm <- do.call(rbind, df_tmm_list)
  df_tmm$Dataset <- factor(df_tmm$Dataset, levels = all_datasets)
  p_tmm <- ggplot2::ggplot(df_tmm, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "TMM", x = "TMM", y = "density") +
    theme_human()

  # 3. Effective library size
  df_eff_list <- list(data.frame(val = lib_r * tmm_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_eff_list[[sn]] <- data.frame(val = colSums(sims[[sn]]) * tmm_sims[[sn]], Dataset = sn)
  }
  df_eff <- do.call(rbind, df_eff_list)
  df_eff$Dataset <- factor(df_eff$Dataset, levels = all_datasets)
  p_eff <- ggplot2::ggplot(df_eff, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Effective library size", x = "effective library size", y = "density") +
    theme_human()

  # 4. Mean expression
  m_r <- rowMeans(ref_mat)
  df_mean_list <- list(data.frame(val = m_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_mean_list[[sn]] <- data.frame(val = rowMeans(sims[[sn]]), Dataset = sn)
  }
  df_mean <- do.call(rbind, df_mean_list)
  df_mean$Dataset <- factor(df_mean$Dataset, levels = all_datasets)
  p_mean <- ggplot2::ggplot(df_mean, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Mean expression", x = "mean expression", y = "density") +
    theme_human()

  # 5. Variance
  v_r <- apply(ref_mat, 1, stats::var)
  df_var_list <- list(data.frame(val = v_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_var_list[[sn]] <- data.frame(val = apply(sims[[sn]], 1, stats::var), Dataset = sn)
  }
  df_var <- do.call(rbind, df_var_list)
  df_var$Dataset <- factor(df_var$Dataset, levels = all_datasets)
  p_var <- ggplot2::ggplot(df_var, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Variance", x = "variance", y = "density") +
    theme_human()

  # 6. Scaled variance (dispersion)
  df_svar_list <- list(data.frame(val = v_r / (m_r + 1e-4), Dataset = ref_name))
  for (sn in sim_names) {
    m_s <- rowMeans(sims[[sn]])
    v_s <- apply(sims[[sn]], 1, stats::var)
    df_svar_list[[sn]] <- data.frame(val = v_s / (m_s + 1e-4), Dataset = sn)
  }
  df_svar <- do.call(rbind, df_svar_list)
  df_svar$Dataset <- factor(df_svar$Dataset, levels = all_datasets)
  p_svar <- ggplot2::ggplot(df_svar, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Scaled variance", x = "scaled variance", y = "density") +
    theme_human()

  # 7. Fraction zero cell
  zc_r <- colMeans(ref_mat == 0)
  df_zeroc_list <- list(data.frame(val = zc_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_zeroc_list[[sn]] <- data.frame(val = colMeans(sims[[sn]] == 0), Dataset = sn)
  }
  df_zeroc <- do.call(rbind, df_zeroc_list)
  df_zeroc$Dataset <- factor(df_zeroc$Dataset, levels = all_datasets)
  p_zeroc <- ggplot2::ggplot(df_zeroc, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Fraction zero cell", x = "fraction zeros per cell", y = "density") +
    theme_human()

  # 8. Fraction zero gene
  zg_r <- rowMeans(ref_mat == 0)
  df_zerog_list <- list(data.frame(val = zg_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_zerog_list[[sn]] <- data.frame(val = rowMeans(sims[[sn]] == 0), Dataset = sn)
  }
  df_zerog <- do.call(rbind, df_zerog_list)
  df_zerog$Dataset <- factor(df_zerog$Dataset, levels = all_datasets)
  p_zerog <- ggplot2::ggplot(df_zerog, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Fraction zero gene", x = "fraction zeros per gene", y = "density") +
    theme_human()

  # 9. Cell correlation
  n_c <- min(ncol(ref_mat), 100)
  c_idx <- seq_len(n_c)
  cc_r <- stats::cor(ref_mat[, c_idx])[upper.tri(diag(n_c))]
  df_ccor_list <- list(data.frame(val = cc_r, Dataset = ref_name))
  for (sn in sim_names) {
    cc_s <- stats::cor(sims[[sn]][, c_idx])[upper.tri(diag(n_c))]
    df_ccor_list[[sn]] <- data.frame(val = cc_s, Dataset = sn)
  }
  df_ccor <- do.call(rbind, df_ccor_list)
  df_ccor$Dataset <- factor(df_ccor$Dataset, levels = all_datasets)
  p_ccor <- ggplot2::ggplot(df_ccor, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Cell correlation", x = "cell correlation", y = "density") +
    theme_human()

  # 10. Gene correlation
  top_g <- order(v_r, decreasing = TRUE)[seq_len(min(nrow(ref_mat), 100))]
  gc_r <- stats::cor(t(ref_mat[top_g, ]))[upper.tri(diag(length(top_g)))]
  df_gcor_list <- list(data.frame(val = gc_r, Dataset = ref_name))
  for (sn in sim_names) {
    gc_s <- stats::cor(t(sims[[sn]][top_g, ]))[upper.tri(diag(length(top_g)))]
    df_gcor_list[[sn]] <- data.frame(val = gc_s, Dataset = sn)
  }
  df_gcor <- do.call(rbind, df_gcor_list)
  df_gcor$Dataset <- factor(df_gcor$Dataset, levels = all_datasets)
  p_gcor <- ggplot2::ggplot(df_gcor, ggplot2::aes(x = .data$val, fill = .data$Dataset, color = .data$Dataset)) +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.65, adjust = 1.1) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Gene correlation", x = "gene correlation", y = "density") +
    theme_human()

  # Return density grid if layout == "density" or patchwork not available
  if (layout == "density" || !requireNamespace("patchwork", quietly = TRUE)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_density_grid <- (p_lib | p_mean | p_var | p_svar) /
                        (p_zeroc | p_zerog | p_ccor | p_gcor) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
      if (!is.null(title) && nzchar(title)) {
        p_density_grid <- p_density_grid + patchwork::plot_annotation(
          title = title,
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.15, hjust = 0.5))
        )
      }
      return(p_density_grid)
    } else {
      return(p_lib)
    }
  }

  # 11. Mean vs variance (scatter)
  df_mv_list <- list(data.frame(m = m_r, v = v_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_mv_list[[sn]] <- data.frame(m = rowMeans(sims[[sn]]), v = apply(sims[[sn]], 1, stats::var), Dataset = sn)
  }
  df_mv <- do.call(rbind, df_mv_list)
  df_mv$Dataset <- factor(df_mv$Dataset, levels = all_datasets)
  p_mv <- ggplot2::ggplot(df_mv, ggplot2::aes(x = .data$m, y = .data$v, color = .data$Dataset)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.75, stroke = 0) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Mean vs variance", x = "mean expression", y = "variance of gene expression") +
    theme_human()

  # 12. Mean vs fraction zero (scatter)
  df_mz_list <- list(data.frame(m = m_r, z = zg_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_mz_list[[sn]] <- data.frame(m = rowMeans(sims[[sn]]), z = rowMeans(sims[[sn]] == 0), Dataset = sn)
  }
  df_mz <- do.call(rbind, df_mz_list)
  df_mz$Dataset <- factor(df_mz$Dataset, levels = all_datasets)
  p_mz <- ggplot2::ggplot(df_mz, ggplot2::aes(x = .data$m, y = .data$z, color = .data$Dataset)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.75, stroke = 0) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Mean vs fraction zero", x = "mean expression", y = "fraction zero per gene") +
    theme_human()

  # 13. Library size vs fraction zero (scatter)
  df_lz_list <- list(data.frame(lib = lib_r, z = zc_r, Dataset = ref_name))
  for (sn in sim_names) {
    df_lz_list[[sn]] <- data.frame(lib = colSums(sims[[sn]]), z = colMeans(sims[[sn]] == 0), Dataset = sn)
  }
  df_lz <- do.call(rbind, df_lz_list)
  df_lz$Dataset <- factor(df_lz$Dataset, levels = all_datasets)
  p_lz <- ggplot2::ggplot(df_lz, ggplot2::aes(x = .data$lib, y = .data$z, color = .data$Dataset)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.75, stroke = 0) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title = "Library size vs fraction zero", x = "library size", y = "fraction zero per cell") +
    theme_human()

  # 14. Biological signals barplot
  calc_signals <- function(mat, groups = NULL) {
    if (is.null(groups) || length(unique(groups)) < 2) {
      lib <- colSums(mat)
      groups <- factor(ifelse(lib >= stats::median(lib), "High", "Low"))
    }
    lvls <- levels(factor(groups))
    g1 <- which(groups == lvls[1])
    g2 <- which(groups == lvls[2])
    if (length(g1) < 3 || length(g2) < 3) return(c(DE = 0.1, DV = 0.1, DD = 0.1, DP = 0.1, BD = 0.1))

    de_cnt <- 0; dv_cnt <- 0; dd_cnt <- 0; dp_cnt <- 0; bd_cnt <- 0
    n_genes <- min(nrow(mat), 300)
    for (i in seq_len(n_genes)) {
      x <- mat[i, g1]; y <- mat[i, g2]
      p_mean <- tryCatch(stats::t.test(x, y)$p.value, error = function(e) 1)
      p_var  <- tryCatch(stats::var.test(x, y)$p.value, error = function(e) 1)
      p_ks   <- tryCatch(suppressWarnings(stats::ks.test(x, y)$p.value), error = function(e) 1)
      z_tab  <- matrix(c(sum(x == 0), sum(x > 0), sum(y == 0), sum(y > 0)), nrow = 2)
      p_dp   <- tryCatch(suppressWarnings(stats::chisq.test(z_tab)$p.value), error = function(e) 1)

      is_mean <- !is.na(p_mean) && p_mean < 0.05
      is_var  <- !is.na(p_var)  && p_var < 0.05
      is_dp   <- !is.na(p_dp)   && p_dp < 0.05
      is_ks   <- !is.na(p_ks)   && p_ks < 0.05

      if (is_mean && is_var) {
        bd_cnt <- bd_cnt + 1
      } else if (is_mean) {
        de_cnt <- de_cnt + 1
      } else if (is_var) {
        dv_cnt <- dv_cnt + 1
      } else if (is_dp) {
        dp_cnt <- dp_cnt + 1
      } else if (is_ks) {
        dd_cnt <- dd_cnt + 1
      }
    }
    c(DE = de_cnt, DV = dv_cnt, DD = dd_cnt, DP = dp_cnt, BD = bd_cnt) / n_genes
  }

  sig_r <- calc_signals(ref_mat, cell_types)
  types <- c("DE", "DV", "DD", "DP", "BD")
  df_bio_list <- list(data.frame(Signal = factor(types, levels = types), Proportion = sig_r, Dataset = ref_name))
  for (sn in sim_names) {
    sig_s <- calc_signals(sims[[sn]], cell_types)
    df_bio_list[[sn]] <- data.frame(Signal = factor(types, levels = types), Proportion = sig_s, Dataset = sn)
  }
  df_bio <- do.call(rbind, df_bio_list)
  df_bio$Dataset <- factor(df_bio$Dataset, levels = all_datasets)
  p_bio <- ggplot2::ggplot(df_bio, ggplot2::aes(x = .data$Signal, y = .data$Proportion, fill = .data$Dataset)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(title = "Biological signals", x = "types of differentially present genes", y = "proportion") +
    theme_human() +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))

  # 15. Clean Standalone Dataset Legend
  df_leg <- data.frame(
    Dataset = factor(all_datasets, levels = all_datasets),
    y = rev(seq_along(all_datasets))
  )
  p_legend <- ggplot2::ggplot(df_leg, ggplot2::aes(y = .data$y)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0.15, xmax = 0.50, ymin = .data$y - 0.28, ymax = .data$y + 0.28, fill = .data$Dataset), color = NA) +
    ggplot2::geom_text(ggplot2::aes(x = 0.65, label = .data$Dataset), hjust = 0, size = base_size * 0.38, color = "black") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.3, length(all_datasets) + 0.9), expand = c(0, 0)) +
    ggplot2::labs(title = "Dataset") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 0.95, hjust = 0.12, margin = ggplot2::margin(b = 6)),
      plot.margin = ggplot2::margin(l = 15, t = 8)
    )

  # Assemble 4x4 layout matching canonical single-cell benchmark publication
  composite <- p_lib + p_tmm + p_eff + p_mean +
    p_var + p_svar + p_zeroc + p_zerog +
    p_ccor + p_gcor + p_mv + p_mz +
    p_lz + p_bio + p_legend +
    patchwork::plot_layout(design = "
ABCD
EFGH
IJKL
MNO#
")

  if (!is.null(title) && nzchar(title)) {
    composite <- composite + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.15, hjust = 0.5, margin = ggplot2::margin(b = 6)))
    )
  }

  composite
}


# ============================================================================
# 2. plot_metric_boxplots()
# ============================================================================

#' Plot Comparison Boxplots Across Simulators
#' Plot Benchmark Metric Score Distributions Across Simulators
#'
#' Generates 600 DPI boxplots with jittered points showing the distribution
#' of evaluation metric scores across simulated single-cell and multiomics datasets.
#' Supports standardized direction-aligned fidelity scores (where higher is universally
#' superior, resolving metric polarity differences across distance and correlation measures)
#' or original unnormalized raw values with optional direction separation or per-metric faceting.
#'
#' \if{html}{\figure{metric_boxplots.png}{options: width="800" alt="Multi-Metric Distribution Boxplots"}}
#'
#' @param benchmark_data A data frame (such as \code{demo$benchmark_summary_table} or output from
#'   \code{\link{evaluate_simulation_accuracy}()}) or a named list of benchmark result tables.
#' @param categories Optional character vector of categories to include. Default \code{NULL} (all 8 categories).
#' @param metrics Optional character vector of specific metrics to include. Default \code{NULL} (all).
#' @param score_type Character. Score formulation: \code{"normalized"} (default; direction-inverted
#'   standardized fidelity scores in [0, 1] where higher is universally superior across all measures)
#'   or \code{"raw"} (unnormalized original metric values).
#' @param separate_direction Logical. When \code{score_type = "raw"}, whether to separate metrics by
#'   optimization direction (lower vs higher is better) into distinct sub-panels. Default \code{FALSE}.
#' @param facet_by Character. Facet layout: \code{"category"} (default; 8 canonical category panels)
#'   or \code{"metric"} (independent panels for each individual metric).
#' @param palette Optional named character vector of simulator colors.
#' @param ncol Integer. Number of columns in facet wrap. Default \code{4}.
#' @param base_size Numeric. Base font size. Default \code{11}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_metric_boxplots <- function(
  benchmark_data,
  categories         = NULL,
  metrics            = NULL,
  score_type         = c("normalized", "raw"),
  separate_direction = FALSE,
  facet_by           = c("category", "metric"),
  palette            = NULL,
  ncol               = 4,
  base_size          = 11
) {
  score_type <- match.arg(score_type)
  facet_by   <- match.arg(facet_by)
  df <- .ingest_bubble_data(benchmark_data)

  if (!is.null(categories)) {
    df <- df[df$Category %in% categories, , drop = FALSE]
  }
  if (!is.null(metrics)) {
    df <- df[df$Metric %in% metrics, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No data remaining after applying category/metric filters.")

  df <- .normalize_bubble_scores(df)

  # Map category names to compact two-line labels to ensure complete visibility
  cat_chr <- as.character(df$Category)
  df$Category <- ifelse(cat_chr %in% names(.CATEGORY_TWO_LINE_MAP),
                        .CATEGORY_TWO_LINE_MAP[cat_chr], cat_chr)
  ordered_cats <- unname(.CATEGORY_TWO_LINE_MAP[.CANONICAL_CATEGORY_ORDER])
  ordered_cats <- c(ordered_cats[ordered_cats %in% unique(df$Category)],
                    setdiff(unique(df$Category), ordered_cats))
  df$Category <- factor(df$Category, levels = ordered_cats)

  is_hib <- df$Metric %in% .HIGHER_IS_BETTER_METRICS
  df$Direction <- ifelse(is_hib, "\u2191 Higher is better", "\u2193 Lower is better")
  disp_m <- ifelse(df$Metric %in% names(.METRIC_DISPLAY_LABELS), .METRIC_DISPLAY_LABELS[df$Metric], df$Metric)
  df$Display_Metric <- paste0(disp_m, ifelse(is_hib, " (\u2191)", " (\u2193)"))

  default_sim_cols <- c(
    "scDesign3" = "#1E8449",
    "Splatter"  = "#2E86AB",
    "SCRIP"     = "#E67E22",
    "SymSim"    = "#8E44AD",
    "dyngen"    = "#D4AC0D",
    "simATAC"   = "#C0392B"
  )
  sim_methods <- unique(df$Method)
  if (is.null(palette)) {
    active_pal <- default_sim_cols[sim_methods]
    missing_sims <- sim_methods[is.na(active_pal)]
    if (length(missing_sims) > 0) {
      extra_cols <- grDevices::hcl.colors(length(missing_sims), palette = "Dark 3")
      names(extra_cols) <- missing_sims
      active_pal[missing_sims] <- extra_cols
    }
  } else {
    active_pal <- palette
  }

  y_var <- if (score_type == "normalized") "Score_Norm" else "Score_Raw"
  y_lab <- if (score_type == "normalized") {
    "Standardized Fidelity Score [0, 1] (\u2191 Better)"
  } else {
    "Original Raw Metric Score"
  }
  sub_txt <- if (score_type == "normalized") {
    "Per-metric min-max normalized fidelity scores in [0, 1] (distance/error inverted; higher is always better)"
  } else {
    "Original raw evaluation metric values (unnormalized; points denote individual evaluated measures)"
  }
  cap_txt <- if (score_type == "normalized") {
    "Boxplot indicates median and IQR across metrics per category; points denote individual evaluated measures.\nMetrics are direction-standardized to [0, 1] across simulators (1.0 = best observed performance; distance/error metrics inverted)."
  } else {
    "Boxplot indicates median and IQR across metrics per category; points denote individual evaluated measures in original units."
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Method, y = .data[[y_var]], fill = .data$Method)) +
    ggplot2::geom_boxplot(alpha = 0.72, outlier.shape = NA, width = 0.65, color = "#2C3E50", linewidth = 0.45) +
    ggplot2::geom_jitter(width = 0.18, height = 0, size = 1.6, alpha = 0.65, shape = 21, color = "#1A252F") +
    ggplot2::scale_fill_manual(values = active_pal, name = "Simulator") +
    ggplot2::labs(
      title    = "Benchmarking Metric Score Distributions Across Simulators",
      subtitle = sub_txt,
      x        = "Simulator Framework",
      y        = y_lab,
      caption  = cap_txt
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1, face = "bold", size = base_size * 0.9),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.15, hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.90, hjust = 0.5, color = "#566573"),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size * 0.82, lineheight = 0.92, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "#F2F4F4", color = "#D5D8DC", linewidth = 0.5)
    )

  if (score_type == "normalized") {
    p <- p + ggplot2::scale_y_continuous(limits = c(-0.02, 1.05), breaks = seq(0, 1, 0.25))
  }

  if (facet_by == "metric") {
    p <- p + ggplot2::facet_wrap(~ .data$Display_Metric, scales = "free_y", ncol = ncol)
  } else if (score_type == "raw" && isTRUE(separate_direction)) {
    p <- p + ggplot2::facet_grid(Category ~ Direction, scales = "free_y")
  } else if (length(unique(df$Category)) > 1) {
    p <- p + ggplot2::facet_wrap(~ .data$Category, scales = "free_y", ncol = ncol)
  }

  p
}


# ============================================================================
# 2a. plot_individual_metric_bar()
# ============================================================================

#' Plot Individual Metric Performance Barplot
#'
#' Generates a publication-quality barplot comparing simulator methods on a single
#' evaluation metric, with exact score labels, direction awareness, and performance ranking.
#'
#' @param benchmark_data A benchmark data.frame or named list containing benchmark tables.
#' @param metric Character. The name of the metric to plot (e.g. \code{"KS"}, \code{"Wasserstein"}, \code{"ARI"}).
#' @param score_type Character. Either \code{"normalized"} (standardized [0, 1] fidelity score; default)
#'   or \code{"raw"} (original unnormalized metric value).
#' @param palette Optional named character vector of simulator colors.
#' @param base_size Numeric. Base font size. Default \code{12}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_individual_metric_bar <- function(
  benchmark_data,
  metric,
  score_type = c("normalized", "raw"),
  palette    = NULL,
  base_size  = 12
) {
  score_type <- match.arg(score_type)
  df <- .ingest_bubble_data(benchmark_data)

  if ("Category" %in% colnames(df)) {
    cat_chr <- as.character(df$Category)
    df$Category <- ifelse(cat_chr %in% names(.LEGACY_CATEGORY_MAP),
                          .LEGACY_CATEGORY_MAP[cat_chr], cat_chr)
  }

  avail_m <- unique(df$Metric)
  if (!metric %in% avail_m) {
    for (m_key in names(.METRIC_DISPLAY_LABELS)) {
      if (identical(.METRIC_DISPLAY_LABELS[[m_key]], metric) && m_key %in% avail_m) {
        metric <- m_key
        break
      }
    }
  }

  df_sub <- df[df$Metric == metric, , drop = FALSE]
  if (nrow(df_sub) == 0) {
    matched <- avail_m[tolower(avail_m) == tolower(metric)]
    if (length(matched) > 0) {
      metric <- matched[1]
      df_sub <- df[df$Metric == metric, , drop = FALSE]
    } else {
      stop(sprintf("Metric '%s' not found in benchmark data.", metric))
    }
  }

  df_all_norm <- .normalize_bubble_scores(df)
  df_sub <- df_all_norm[df_all_norm$Metric == metric, , drop = FALSE]

  if (any(duplicated(df_sub$Method))) {
    df_sub <- stats::aggregate(
      cbind(Score_Norm, Score_Raw) ~ Method + Metric + Category,
      data = df_sub, FUN = mean, na.rm = TRUE
    )
  }

  disp_m <- if (metric %in% names(.METRIC_DISPLAY_LABELS)) .METRIC_DISPLAY_LABELS[metric] else metric
  is_hib <- metric %in% .HIGHER_IS_BETTER_METRICS
  dir_txt <- if (is_hib) "\u2191 Higher is better" else "\u2193 Lower is better"
  cat_name <- if (nrow(df_sub) > 0 && !is.na(df_sub$Category[1])) as.character(df_sub$Category[1]) else "Benchmark"

  default_sim_cols <- c(
    "scDesign3" = "#1E8449",
    "Splatter"  = "#2E86AB",
    "SCRIP"     = "#E67E22",
    "SymSim"    = "#8E44AD",
    "dyngen"    = "#D4AC0D",
    "simATAC"   = "#C0392B"
  )
  sim_methods <- unique(df_sub$Method)
  if (is.null(palette)) {
    active_pal <- default_sim_cols[sim_methods]
    missing_sims <- sim_methods[is.na(active_pal)]
    if (length(missing_sims) > 0) {
      extra_cols <- grDevices::hcl.colors(length(missing_sims), palette = "Dark 3")
      names(extra_cols) <- missing_sims
      active_pal[missing_sims] <- extra_cols
    }
  } else {
    active_pal <- palette
  }

  # Order simulators along the X-axis by performance (best performer first on the left)
  if (score_type == "normalized") {
    df_sub <- df_sub[order(-df_sub$Score_Norm), ]
  } else {
    if (is_hib) {
      df_sub <- df_sub[order(-df_sub$Score_Raw), ]
    } else {
      df_sub <- df_sub[order(df_sub$Score_Raw), ]
    }
  }
  df_sub$Method <- factor(df_sub$Method, levels = df_sub$Method)

  y_var <- if (score_type == "normalized") "Score_Norm" else "Score_Raw"
  y_lab <- if (score_type == "normalized") {
    "Standardized Fidelity Score [0, 1] (\u2191 Higher is Better)"
  } else {
    sprintf("Raw Score Value: %s (%s)", disp_m, dir_txt)
  }

  title_txt <- sprintf("Benchmark Barplot: %s", disp_m)
  sub_txt <- sprintf("Category: %s   |   Optimization: %s", cat_name, dir_txt)

  cap_txt <- if (score_type == "normalized") {
    "Direction-aware standardized score in [0, 1] (1.0 = best agreement with reference; distance/error inverted).\nSimulators ranked by fidelity from best performer (left) to lowest (right)."
  } else {
    sprintf("Unnormalized raw metric values in original units (%s).\nSimulators ranked from best performer (left) to lowest (right).", dir_txt)
  }

  df_sub$Value_Label <- if (score_type == "normalized") {
    sprintf("%.3f", df_sub$Score_Norm)
  } else {
    ifelse(abs(df_sub$Score_Raw) >= 100, sprintf("%.1f", df_sub$Score_Raw),
    ifelse(abs(df_sub$Score_Raw) >= 1,   sprintf("%.3f", df_sub$Score_Raw),
    sprintf("%.4f", df_sub$Score_Raw)))
  }

  p <- ggplot2::ggplot(df_sub, ggplot2::aes(x = Method, y = .data[[y_var]], fill = Method)) +
    ggplot2::geom_col(width = 0.58, color = "#1E293B", linewidth = 0.5, alpha = 0.90) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["Value_Label"]]),
      vjust = ifelse(df_sub[[y_var]] >= 0, -0.45, 1.2),
      fontface = "bold", size = base_size * 0.32, color = "#0F172A"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
    ggplot2::labs(
      title    = title_txt,
      subtitle = sub_txt,
      x        = "Simulator Framework",
      y        = y_lab,
      caption  = cap_txt
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(face = "bold", size = base_size * 0.92, color = "#0F172A"),
      axis.text.y      = ggplot2::element_text(face = "bold", size = base_size * 0.88),
      axis.title.x     = ggplot2::element_text(face = "bold", size = base_size * 0.95, margin = ggplot2::margin(t = 10)),
      axis.title.y     = ggplot2::element_text(face = "bold", size = base_size * 0.95, margin = ggplot2::margin(r = 10)),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.25, hjust = 0.5, color = "#0F172A"),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.92, hjust = 0.5, color = "#475569", margin = ggplot2::margin(b = 10)),
      plot.caption     = ggplot2::element_text(size = base_size * 0.78, color = "#64748B", hjust = 0, margin = ggplot2::margin(t = 8)),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (score_type == "normalized") {
    p <- p +
      ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed", color = "#10B981", linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", color = "#94A3B8", linewidth = 0.5) +
      ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0, 0.05)))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15)))
  }

  p
}


# ============================================================================
# 2ab. plot_category_metric_bars()
# ============================================================================

#' Plot Benchmark Metric Performance Barplots by Category
#'
#' Generates publication-ready faceted barplots comparing simulator methods across all individual
#' evaluation metrics belonging to a specified benchmark category. Each metric is presented in its own
#' sub-panel with exact numeric score labels, direction-awareness indicators (indicating whether higher
#' or lower values are optimal), and standardized or original raw score formulation.
#'
#' @param benchmark_data A benchmark data.frame (e.g. \code{demo$benchmark_summary_table} or output from
#'   \code{\link{evaluate_simulation_accuracy}()}) or a named list of benchmark result tables.
#' @param category Character. The benchmark category to visualize (e.g., \code{"(I) Distributional Properties"}
#'   or \code{"Distributional Properties"}).
#' @param score_type Character. Either \code{"normalized"} (standardized [0, 1] fidelity score where 1.0 is optimal; default)
#'   or \code{"raw"} (unnormalized original metric values in native measurement units).
#' @param palette Optional named character vector of simulator colors.
#' @param ncol Integer. Optional number of columns in facet layout. Defaults to dynamic layout based on metric count.
#' @param base_size Numeric. Base font size. Default \code{11}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_category_metric_bars <- function(
  benchmark_data,
  category,
  score_type = c("normalized", "raw"),
  palette    = NULL,
  ncol       = NULL,
  base_size  = 11
) {
  score_type <- match.arg(score_type)
  df <- .ingest_bubble_data(benchmark_data)

  if ("Category" %in% colnames(df)) {
    cat_chr <- as.character(df$Category)
    df$Category <- ifelse(cat_chr %in% names(.LEGACY_CATEGORY_MAP),
                          .LEGACY_CATEGORY_MAP[cat_chr], cat_chr)
  }

  if (!missing(category) && !is.null(category) && !identical(category, "all")) {
    if (category %in% names(.LEGACY_CATEGORY_MAP)) {
      category <- .LEGACY_CATEGORY_MAP[[category]]
    }
    df <- df[df$Category == category, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop(sprintf("No data found for category '%s'.", category))
  }

  df <- .normalize_bubble_scores(df)

  if (any(duplicated(df[, c("Method", "Metric")]))) {
    df <- stats::aggregate(
      cbind(Score_Norm, Score_Raw) ~ Method + Metric + Category,
      data = df, FUN = mean, na.rm = TRUE
    )
  }

  is_hib <- df$Metric %in% .HIGHER_IS_BETTER_METRICS
  disp_m <- ifelse(df$Metric %in% names(.METRIC_DISPLAY_LABELS), .METRIC_DISPLAY_LABELS[df$Metric], df$Metric)
  dir_symbol <- ifelse(is_hib, "\u2191", "\u2193")
  df$Display_Metric <- paste0(disp_m, " (", dir_symbol, ")")

  unique_metrics <- unique(df$Metric)
  canonical_order <- names(.METRIC_DISPLAY_LABELS)
  ordered_m <- intersect(canonical_order, unique_metrics)
  ordered_m <- c(ordered_m, setdiff(unique_metrics, ordered_m))

  disp_levels <- character(length(ordered_m))
  for (i in seq_along(ordered_m)) {
    m_i <- ordered_m[i]
    dm <- if (m_i %in% names(.METRIC_DISPLAY_LABELS)) .METRIC_DISPLAY_LABELS[[m_i]] else m_i
    ds <- if (m_i %in% .HIGHER_IS_BETTER_METRICS) "\u2191" else "\u2193"
    disp_levels[i] <- paste0(dm, " (", ds, ")")
  }
  df$Display_Metric <- factor(df$Display_Metric, levels = disp_levels)

  default_sim_cols <- c(
    "scDesign3" = "#1E8449",
    "Splatter"  = "#2E86AB",
    "SCRIP"     = "#E67E22",
    "SymSim"    = "#8E44AD",
    "dyngen"    = "#D4AC0D",
    "simATAC"   = "#C0392B"
  )
  sim_methods <- unique(df$Method)
  if (is.null(palette)) {
    active_pal <- default_sim_cols[sim_methods]
    missing_sims <- sim_methods[is.na(active_pal)]
    if (length(missing_sims) > 0) {
      extra_cols <- grDevices::hcl.colors(length(missing_sims), palette = "Dark 3")
      names(extra_cols) <- missing_sims
      active_pal[missing_sims] <- extra_cols
    }
  } else {
    active_pal <- palette
  }

  y_var <- if (score_type == "normalized") "Score_Norm" else "Score_Raw"
  df$Value_Label <- if (score_type == "normalized") {
    sprintf("%.2f", df$Score_Norm)
  } else {
    ifelse(abs(df$Score_Raw) >= 100, sprintf("%.1f", df$Score_Raw),
    ifelse(abs(df$Score_Raw) >= 1,   sprintf("%.2f", df$Score_Raw),
    sprintf("%.3f", df$Score_Raw)))
  }

  n_m <- length(unique(df$Metric))
  ncol_use <- if (!is.null(ncol)) ncol else if (n_m <= 4) n_m else 4

  title_txt <- sprintf("Category Metric Benchmark: %s", category)
  sub_txt <- if (score_type == "normalized") {
    "Standardized Fidelity Scores in [0, 1] per simulator (\u2191 Higher is Better; distance/error inverted)"
  } else {
    "Original unnormalized raw metric scores (\u2191 = Higher is better, \u2193 = Lower is better)"
  }
  cap_txt <- if (score_type == "normalized") {
    "Bars represent standardized fidelity scores [0, 1] across evaluated metrics (1.0 = best observed agreement with reference).\nArrows in panel headers indicate optimization polarity: (\u2191) higher is better, (\u2193) lower is better."
  } else {
    "Bars represent unnormalized raw metric scores in native measurement units.\nArrows in panel headers indicate optimization polarity: (\u2191) higher is better, (\u2193) lower is better."
  }
  y_lab <- if (score_type == "normalized") "Fidelity Score [0, 1] (\u2191 Better)" else "Raw Metric Score"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Method, y = .data[[y_var]], fill = .data$Method)) +
    ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.88) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Value_Label),
      vjust = ifelse(df[[y_var]] >= 0, -0.35, 1.2),
      size = base_size * 0.26, fontface = "bold", color = "#1E293B"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, name = "Simulator") +
    ggplot2::facet_wrap(~ Display_Metric, scales = "free_y", ncol = ncol_use) +
    ggplot2::labs(
      title    = title_txt,
      subtitle = sub_txt,
      x        = "Simulator Framework",
      y        = y_lab,
      caption  = cap_txt
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1, face = "bold", size = base_size * 0.82),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.15, hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.88, hjust = 0.5, color = "#566573"),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size * 0.78, lineheight = 0.92, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "#F2F4F4", color = "#D5D8DC", linewidth = 0.5),
      panel.spacing    = grid::unit(0.7, "lines")
    )

  if (score_type == "normalized") {
    p <- p +
      ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed", color = "#10B981", linewidth = 0.5, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", color = "#94A3B8", linewidth = 0.5, alpha = 0.7) +
      ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22)))
  }

  p
}


# ============================================================================
# 2b. plot_scalability_benchmark()
# ============================================================================

#' Plot Multi-Dimensional Computational Scalability Benchmark
#'
#' Evaluates computational resource efficiency across single-cell simulators using raw
#' benchmarking measurements (elapsed real wall-clock time, total CPU time, peak resident
#' Plot Multi-Dimensional Computational Scalability Benchmark
#'
#' Evaluates computational resource efficiency across single-cell simulators using raw
#' benchmarking measurements (elapsed real wall-clock time in seconds and peak resident
#' memory consumption in MiB).
#'
#' \if{html}{\figure{scalability_benchmark.png}{options: width="800" alt="Computational Scalability and Resource Footprint"}}
#'
#' @param benchmark_data Benchmark summary table or named list containing Category VIII metrics.
#' @param type Character. Visualization type: \code{"composite"} (multi-panel dashboard),
#'   \code{"runtime"} (elapsed wall-clock time barplot), \code{"memory"} (peak RAM barplot),
#'   \code{"tradeoff"} (runtime vs memory scatter/line curve), \code{"cost"} (combined resource footprint index),
#'   \code{"throughput"} (cells/sec), or \code{"efficiency"} (CPU/elapsed concurrency ratio if CPU time available).
#'   Default \code{"composite"}.
#' @param layout Character. For composite dashboards: \code{"4panel"} (runtime, memory, tradeoff, cost footprint)
#'   or \code{"6panel"} (adds throughput). Default \code{"4panel"}.
#' @param cell_count Integer. Number of cells processed in the benchmark, used to calculate
#'   simulation throughput (\code{cells / second}). Default \code{1000}.
#' @param palette Optional named character vector of simulator colors.
#' @param base_size Numeric. Base font size. Default \code{11}.
#'
#' @return A \code{ggplot} or \code{patchwork} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_scalability_benchmark <- function(
  benchmark_data,
  type       = c("composite", "runtime", "memory", "tradeoff", "cost", "throughput", "efficiency"),
  layout     = c("4panel", "6panel"),
  cell_count = 1000,
  palette    = NULL,
  base_size  = 11
) {
  type   <- match.arg(type)
  layout <- match.arg(layout)

  df <- .ingest_bubble_data(benchmark_data)
  df$Metric[df$Metric == "elapsed_time"] <- "elapsed_time_seconds"
  df$Metric[df$Metric == "memory_mb"] <- "peak_memory_mb"
  sc_mets <- c("elapsed_time_seconds", "peak_memory_mb")
  sc_df <- df[df$Metric %in% sc_mets, , drop = FALSE]
  if (nrow(sc_df) == 0) stop("No computational scalability metrics (Category VIII) found in benchmark_data.")

  wide_sc <- stats::reshape(
    sc_df[, c("Method", "Metric", "Score_Raw")],
    timevar   = "Metric",
    idvar     = "Method",
    direction = "wide"
  )
  colnames(wide_sc) <- gsub("^Score_Raw\\.", "", colnames(wide_sc))

  for (m_req in c("elapsed_time_seconds", "peak_memory_mb")) {
    if (!m_req %in% colnames(wide_sc)) wide_sc[[m_req]] <- NA_real_
  }

  wide_sc$resource_cost        <- (wide_sc$elapsed_time_seconds * wide_sc$peak_memory_mb) / 1000
  wide_sc$throughput_cells_sec <- cell_count / pmax(0.1, wide_sc$elapsed_time_seconds)
  if ("cpu_time_seconds" %in% colnames(wide_sc)) {
    wide_sc$cpu_efficiency     <- wide_sc$cpu_time_seconds / pmax(0.1, wide_sc$elapsed_time_seconds)
  } else {
    wide_sc$cpu_efficiency     <- NA_real_
  }

  default_sim_cols <- c(
    "scDesign3" = "#1E8449",
    "Splatter"  = "#2E86AB",
    "SCRIP"     = "#E67E22",
    "SymSim"    = "#8E44AD",
    "dyngen"    = "#D4AC0D",
    "simATAC"   = "#C0392B"
  )
  sim_methods <- unique(wide_sc$Method)
  if (is.null(palette)) {
    active_pal <- default_sim_cols[sim_methods]
    missing_sims <- sim_methods[is.na(active_pal)]
    if (length(missing_sims) > 0) {
      extra_cols <- grDevices::hcl.colors(length(missing_sims), palette = "Dark 3")
      names(extra_cols) <- missing_sims
      active_pal[missing_sims] <- extra_cols
    }
  } else {
    active_pal <- palette
  }

  # Order methods by elapsed time (fastest first)
  wide_sc <- wide_sc[order(wide_sc$elapsed_time_seconds, decreasing = FALSE), ]
  wide_sc$Method <- factor(wide_sc$Method, levels = rev(wide_sc$Method))

  # 1. Panel A: Elapsed Wall-Clock Runtime Barplot
  p_runtime <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$elapsed_time_seconds, y = .data$Method, fill = .data$Method)) +
    ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f s", .data$elapsed_time_seconds)),
      hjust = -0.15, size = base_size * 0.28, fontface = "bold", color = "#1E293B"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
    ggplot2::labs(
      title    = "A. Elapsed Execution Time",
      subtitle = "Wall-clock simulation elapsed time in seconds (lower is faster)",
      x        = "Elapsed Time (Seconds, Lower = Faster)",
      y        = "Simulator"
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))

  if (type == "runtime") return(p_runtime)

  # 2. Panel B: Peak Memory Barplot
  p_memory <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$peak_memory_mb, y = .data$Method, fill = .data$Method)) +
    ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f MiB", .data$peak_memory_mb)),
      hjust = -0.15, size = base_size * 0.28, fontface = "bold", color = "#1E293B"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
    ggplot2::labs(
      title    = "B. Peak Memory Usage",
      subtitle = "Maximum resident memory consumption in MiB (lower is better)",
      x        = "Peak RAM (MiB, Lower = More Efficient)",
      y        = NULL
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))

  if (type == "memory") return(p_memory)

  # 3. Panel C: Runtime vs Memory Tradeoff
  med_x <- stats::median(wide_sc$elapsed_time_seconds, na.rm = TRUE)
  med_y <- stats::median(wide_sc$peak_memory_mb, na.rm = TRUE)

  p_tradeoff <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$elapsed_time_seconds, y = .data$peak_memory_mb, color = .data$Method)) +
    ggplot2::geom_vline(xintercept = med_x, linetype = "dashed", color = "#CBD5E1", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = med_y, linetype = "dashed", color = "#CBD5E1", linewidth = 0.5) +
    ggplot2::geom_point(size = 4.5, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Method), vjust = -1.1, fontface = "bold", size = base_size * 0.30, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = active_pal, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.18))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.18))) +
    ggplot2::labs(
      title    = "C. Runtime vs. Memory Trade-Off",
      subtitle = "Bivariate resource consumption frontier (dashed lines = medians)",
      x        = "Elapsed Time (Seconds, Lower = Faster)",
      y        = "Peak RAM (MiB, Lower = Less Memory)"
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))

  if (type == "tradeoff") return(p_tradeoff)

  # 4. Panel D: Resource Cost Footprint
  p_cost <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$resource_cost, y = .data$Method, fill = .data$Method)) +
    ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f kMiB-s", .data$resource_cost)),
      hjust = -0.15, size = base_size * 0.28, fontface = "bold", color = "#1E293B"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(
      title    = "D. Resource Cost Footprint",
      subtitle = "Combined time-memory footprint: (Elapsed Time * Peak RAM) / 1000",
      x        = "Resource Footprint (kilo-MiB-seconds, Lower = Cheaper)",
      y        = NULL
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))

  if (type == "cost") return(p_cost)

  # 5. Panel E: Throughput
  p_thru <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$throughput_cells_sec, y = .data$Method, fill = .data$Method)) +
    ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f cells/s", .data$throughput_cells_sec)),
      hjust = -0.15, size = base_size * 0.28, fontface = "bold", color = "#1E293B"
    ) +
    ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(
      title    = "Simulation Throughput",
      subtitle = sprintf("Cells generated per elapsed wall-clock second (benchmark: %d cells)", cell_count),
      x        = "Throughput (Cells / Second, Higher = Faster)",
      y        = NULL
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))

  if (type == "throughput") return(p_thru)

  # Optional Efficiency panel (if CPU time is provided)
  if (type == "efficiency") {
    if (all(is.na(wide_sc$cpu_efficiency))) {
      message("Note: cpu_efficiency requires CPU execution time; displaying simulation throughput instead.")
      return(p_thru)
    }
    p_eff <- ggplot2::ggplot(wide_sc, ggplot2::aes(x = .data$cpu_efficiency, y = .data$Method, fill = .data$Method)) +
      ggplot2::geom_vline(xintercept = 1.0, linetype = "dashed", color = "#64748B", linewidth = 0.6) +
      ggplot2::geom_col(width = 0.65, color = "#1E293B", linewidth = 0.4, alpha = 0.9) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2fx", .data$cpu_efficiency)),
        hjust = -0.15, size = base_size * 0.28, fontface = "bold", color = "#1E293B"
      ) +
      ggplot2::scale_fill_manual(values = active_pal, guide = "none") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
      ggplot2::labs(
        title    = "CPU Multithreading Efficiency",
        subtitle = "Ratio of CPU Time to Elapsed Time (> 1.0x indicates multicore concurrency)",
        x        = "Efficiency Ratio (CPU / Elapsed)",
        y        = NULL
      ) +
      .pub_theme(base_size = base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0), plot.subtitle = ggplot2::element_text(hjust = 0))
    return(p_eff)
  }

  # Composite Assembly
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(p_runtime)
  }

  if (layout == "6panel") {
    composite <- (p_runtime | p_memory) / (p_tradeoff | p_cost) / (p_thru | patchwork::plot_spacer()) +
      patchwork::plot_annotation(
        title    = "Computational Scalability & Resource Utilization Benchmark",
        subtitle = "Multi-panel profiling of execution runtime, peak memory footprint, trade-off frontier, and computational cost",
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = base_size * 1.35, hjust = 0.5, color = "#0F172A"),
          plot.subtitle = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5, color = "#475569", margin = ggplot2::margin(b = 10))
        )
      )
  } else {
    composite <- (p_runtime | p_memory) / (p_tradeoff | p_cost) +
      patchwork::plot_annotation(
        title    = "Computational Scalability & Resource Utilization Benchmark",
        subtitle = "Multi-panel profiling of elapsed execution time, peak resident RAM, and resource tradeoff frontier",
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0.5, color = "#0F172A"),
          plot.subtitle = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5, color = "#475569", margin = ggplot2::margin(b = 10))
        )
      )
  }

  composite
}


# ============================================================================
# 3. plot_benchmark_summary()
# ============================================================================

#' Plot Benchmark Summary Metrics (Ranked Bar Chart)
#'
#' Displays a ranked horizontal bar chart of evaluated benchmarking discrepancy
#' metrics across biological properties, colored by canonical evaluation category.
#'
#' @param benchmark_results A data frame from \code{evaluate_multiomics_accuracy()$benchmark_summary_table}
#'   or \code{evaluate_simulation_accuracy()$metrics_summary_table}, or the full result list.
#' @param category Optional character string or vector to filter by evaluation category.
#' @param modality Optional character string to filter by omics modality.
#' @param metric_type Character. Distance metric to plot. Default \code{"KS"}.
#' @param top_n Integer. Show only top N properties. Default \code{NULL} (all).
#' @param palette Optional named color vector. If \code{NULL}, uses canonical palette.
#'
#' @return A \code{ggplot} object of sorted discrepancy scores.
#' @keywords internal
#' @noRd
#' @examples
#' data(example_multiomics)
#' res <- evaluate_multiomics_accuracy(
#'   example_multiomics$ref_multi,
#'   example_multiomics$sim_multi,
#'   verbose = FALSE
#' )
#' p <- plot_benchmark_summary(res$benchmark_summary_table, metric_type = "KS")
#' if (requireNamespace("ggplot2", quietly = TRUE)) print(p)
plot_benchmark_summary <- function(
  benchmark_results,
  category    = NULL,
  modality    = NULL,
  metric_type = "KS",
  top_n       = NULL,
  palette     = NULL
) {
  df <- if (is.data.frame(benchmark_results)) {
    benchmark_results
  } else if (is.list(benchmark_results) && !is.null(benchmark_results$benchmark_summary_table)) {
    benchmark_results$benchmark_summary_table
  } else if (is.list(benchmark_results) && !is.null(benchmark_results$metrics_summary_table)) {
    benchmark_results$metrics_summary_table
  } else {
    stop("benchmark_results must be a data.frame or benchmark output list.")
  }

  if ("Category" %in% colnames(df)) {
    df$Category <- ifelse(df$Category %in% names(.LEGACY_CATEGORY_MAP),
                          .LEGACY_CATEGORY_MAP[df$Category], df$Category)
  }
  if (!is.null(category)) df <- df[df$Category %in% category, , drop = FALSE]
  if (!is.null(modality) && "Modality" %in% colnames(df))
    df <- df[df$Modality %in% modality, , drop = FALSE]

  sub_df <- df[df$Metric == metric_type, , drop = FALSE]
  if (nrow(sub_df) == 0) { sub_df <- head(df, 30); metric_type <- "All Metrics" }

  sub_df$Property_Clean <- gsub("_", " ", sub_df$Property)
  sub_df <- sub_df[order(abs(sub_df$Value), decreasing = TRUE), ]
  if (!is.null(top_n)) sub_df <- head(sub_df, top_n)
  sub_df$Property_Clean <- factor(sub_df$Property_Clean,
                                   levels = rev(unique(sub_df$Property_Clean)))

  fill_col <- if ("Category" %in% colnames(sub_df)) "Category" else
              if ("Modality" %in% colnames(sub_df)) "Modality" else "Property_Clean"

  active_palette <- if (!is.null(palette)) palette else
                    if (fill_col == "Category") .CANONICAL_CATEGORY_COLORS else NULL

  p <- ggplot2::ggplot(sub_df,
                       ggplot2::aes(x = .data$Value, y = .data$Property_Clean,
                                    fill = .data[[fill_col]])) +
    ggplot2::geom_col(width = 0.72, alpha = 0.88, color = "#2C3E50", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$Value)),
                       hjust = -0.12, size = 3.0, color = "#2C3E50", fontface = "bold") +
    ggplot2::geom_vline(xintercept = 0, color = "#2C3E50", linewidth = 0.5, linetype = "dashed") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title    = paste0("Benchmark Accuracy Scores: ", metric_type, " Distance"),
      subtitle = "Ranked discrepancy between real and simulated distributions (lower = better fidelity)",
      x        = paste0("Discrepancy (", metric_type, ")"),
      y        = "Evaluated Biological Property",
      fill     = fill_col,
      caption  = "Bars sorted by absolute discrepancy magnitude."
    ) +
    .pub_theme(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold", size = 9),
                   legend.position = if (length(unique(sub_df[[fill_col]])) > 1) "right" else "none")

  if (!is.null(active_palette)) {
    p <- p + ggplot2::scale_fill_manual(values = active_palette)
  } else if (fill_col != "Category") {
    p <- p + ggplot2::scale_fill_brewer(palette = "Set2")
  }
  p
}






# ============================================================================
# 6. plot_metric_heatmap()
# ============================================================================

#' Plot Multi-Simulator Comparative Metric Heatmap Across Canonical Categories
#'
#' Produces a 600 DPI comparative heatmap displaying simulators on the x-axis
#' and evaluated biological/computational measures on the y-axis, grouped by the 8 canonical
#' evaluation categories. Displays the exact original raw evaluation score in text inside every cell,
#' with cell fill colors scaled by relative fidelity (direction-aware, ensuring balanced visual contrast
#' across all 62 measures without scale distortion from high-magnitude metrics like RAM or runtime).
#'
#' \if{html}{\figure{metric_heatmap.png}{options: width="800" alt="Cross-Metric Z-Scored Heatmap"}}
#'
#' @param benchmark_data A tidy benchmark summary table or a named list of benchmark tables.
#' @param category Optional character vector to filter by canonical evaluation category. Default \code{NULL} (all 8 categories).
#' @param metrics Optional character vector of specific metrics to include. Default \code{NULL} (all).
#' @param scale_fill Character. Cell color fill scaling: \code{"relative"} (default; direction-aware
#'   min-max normalized fidelity in [0, 1] across methods per metric, ensuring balanced color contrast
#'   across all 62 measures while cell text displays exact raw numbers) or \code{"raw"} (direct numeric values).
#' @param top_n_properties Optional integer. If specified, restricts to top N most variable metrics. Default \code{NULL} (all metrics).
#' @param facet_by_category Logical. Whether to group rows into vertical category banner panels across all 8 categories. Default \code{TRUE}.
#' @param cluster_rows Logical. Hierarchically cluster measures along y-axis. Default \code{FALSE}.
#' @param cluster_cols Logical. Hierarchically cluster simulators along x-axis. Default \code{FALSE}.
#' @param base_size Numeric. Base font size. Default \code{10}.
#' @param metric_name Optional legacy alias for \code{metrics}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_metric_heatmap <- function(
  benchmark_data,
  category          = NULL,
  metrics           = NULL,
  scale_fill        = c("relative", "raw"),
  top_n_properties  = NULL,
  facet_by_category = TRUE,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  base_size         = 10,
  metric_name       = NULL
) {
  if (!is.null(metric_name) && is.null(metrics)) metrics <- metric_name
  scale_fill <- match.arg(scale_fill)
  df <- .ingest_bubble_data(benchmark_data)

  if (!is.null(category)) {
    df <- df[df$Category %in% category, , drop = FALSE]
  }
  if (!is.null(metrics)) {
    df <- df[df$Metric %in% metrics, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No data remaining after applying category/metric filters.")

  # Clean measure display names with polarity indicators
  df$Display_Metric <- if ("Metric" %in% colnames(df)) {
    ifelse(df$Metric %in% names(.METRIC_DISPLAY_LABELS), .METRIC_DISPLAY_LABELS[df$Metric], df$Metric)
  } else df$Metric

  is_hib <- df$Metric %in% .HIGHER_IS_BETTER_METRICS
  df$Display_Metric <- paste0(df$Display_Metric, ifelse(is_hib, " (+)", " (-)"))

  # Normalize per metric for relative color scaling
  df <- .normalize_bubble_scores(df)

  # Filter top discriminating properties only if explicitly requested
  all_props <- unique(df$Display_Metric)
  if (!is.null(top_n_properties) && length(all_props) > top_n_properties) {
    var_by_prop <- tapply(df$Score_Raw, df$Display_Metric, stats::var, na.rm = TRUE)
    top_props   <- names(sort(var_by_prop, decreasing = TRUE))[seq_len(top_n_properties)]
    df          <- df[df$Display_Metric %in% top_props, , drop = FALSE]
  }

  # Text color based on lightness & exact raw number formatting
  df$Score_Label <- ifelse(is.na(df$Score_Raw), "-",
                           ifelse(abs(df$Score_Raw) >= 100, sprintf("%.0f", df$Score_Raw),
                           ifelse(abs(df$Score_Raw) >= 10,  sprintf("%.1f", df$Score_Raw),
                           sprintf("%.2f", df$Score_Raw))))

  row_order <- unique(df$Display_Metric)
  col_order <- unique(df$Method)

  # Reshape wide for optional clustering
  if ((cluster_rows || cluster_cols) && length(row_order) > 1 && length(col_order) > 1) {
    wide_m <- stats::reshape(
      df[, c("Method", "Display_Metric", "Score_Raw")],
      timevar = "Method", idvar = "Display_Metric", direction = "wide"
    )
    rnames <- wide_m$Display_Metric
    mat <- as.matrix(wide_m[, -1, drop = FALSE])
    rownames(mat) <- rnames
    colnames(mat) <- gsub("^Score_Raw\\.", "", colnames(mat))
    mat[is.na(mat)] <- 0

    if (cluster_rows && nrow(mat) > 1) {
      hc_r <- stats::hclust(stats::dist(mat))
      row_order <- rnames[hc_r$order]
    }
    if (cluster_cols && ncol(mat) > 1) {
      hc_c <- stats::hclust(stats::dist(t(mat)))
      col_order <- colnames(mat)[hc_c$order]
    }
  }

  df$Display_Metric <- factor(df$Display_Metric, levels = rev(row_order))
  df$Method         <- factor(df$Method, levels = col_order)

  # Map category names to compact two-line labels to reduce right-side panel width
  cat_chr <- as.character(df$Category)
  df$Category <- ifelse(cat_chr %in% names(.CATEGORY_TWO_LINE_MAP),
                        .CATEGORY_TWO_LINE_MAP[cat_chr], cat_chr)
  ordered_cats <- unname(.CATEGORY_TWO_LINE_MAP[.CANONICAL_CATEGORY_ORDER])
  ordered_cats <- c(ordered_cats[ordered_cats %in% unique(df$Category)],
                    setdiff(unique(df$Category), ordered_cats))
  df$Category <- factor(df$Category, levels = ordered_cats)

  fill_var <- if (scale_fill == "relative") "Score_Norm" else "Score_Raw"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Method, y = .data$Display_Metric, fill = .data[[fill_var]])) +
    ggplot2::geom_tile(color = "#FFFFFF", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Score_Label),
      size = base_size * 0.26, fontface = "bold", color = "#1E293B"
    )

  if (scale_fill == "relative") {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#EBF5FB", high = "#2E86AB",
      na.value = "#F2F4F4", name = "Relative Fidelity [0, 1]\n(Higher = Better)",
      limits = c(0, 1),
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(7, "cm"), barheight = ggplot2::unit(0.35, "cm"),
        title.position = "top", title.hjust = 0.5
      )
    )
  } else {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#EBF5FB", mid = "#AED6F1", high = "#2E86AB",
      na.value = "#F2F4F4", name = "Raw Score",
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(7, "cm"), barheight = ggplot2::unit(0.35, "cm"),
        title.position = "top", title.hjust = 0.5
      )
    )
  }

  p <- p +
    ggplot2::labs(
      title    = "Multi-Simulator Benchmark Heatmap Across Canonical Categories",
      subtitle = "Cell text displays exact unnormalized raw scores; cell fill color indicates relative performance",
      x        = "Simulator Framework",
      y        = "Evaluated Biological & Computational Metric",
      caption  = "Values displayed inside cells are original unnormalized scores. Fill indicates standardized fidelity [0, 1] (1.0 = best). (+) Higher is better, (-) Lower is better."
    ) +
    .pub_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1, face = "bold", size = base_size * 0.90),
      axis.text.y      = ggplot2::element_text(size = base_size * 0.75),
      panel.grid       = ggplot2::element_blank(),
      strip.text.y     = ggplot2::element_text(angle = 0, face = "bold", size = base_size * 0.70, lineheight = 0.90, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "#F4F6F7", color = "#BDC3C7", linewidth = 0.4),
      legend.position  = "bottom"
    )

  if (isTRUE(facet_by_category) && length(unique(df$Category)) > 1) {
    p <- p + ggplot2::facet_grid(Category ~ ., scales = "free_y", space = "free_y")
  }

  p
}


# ============================================================================
# 6b. plot_metric_mds()
# ============================================================================

#' Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or Simulators
#'
#' Projects evaluation metric profiles or simulator performances into a 2D MDS space
#' for benchmarking single-cell simulation methods.
#'
#' When \code{ordination_by = "summaries"}, each point represents an evaluated property/metric,
#' colored by its biological level: \code{"gene"} (red), \code{"cell"} (blue), or \code{"global"} (green).
#' If \code{by_category = TRUE}, a multi-panel dashboard is produced showing one MDS ordination per
#' category (by default across the 6 canonical categories with >= 3 measures, excluding Trajectory and
#' Scalability which each contain only 2 measures).
#'
#' @param benchmark_data Benchmark summary data.frame, matrix, or list.
#' @param ordination_by Character. Either \code{"summaries"} (default; ordinated points represent
#'   biological and computational summaries) or \code{"simulators"} (ordinated points represent
#'   candidate simulators).
#' @param by_category Logical. If \code{TRUE} and \code{ordination_by = "summaries"}, generates
#'   a multi-panel grid displaying an MDS ordination for each evaluation category separately
#'   (skipping categories with fewer than 3 measures). Default is \code{FALSE}.
#' @param as_list Logical. If \code{TRUE} and \code{by_category = TRUE}, returns a named list of ggplot objects. Default \code{FALSE}.
#' @param category Optional character vector to filter by evaluation category (e.g. \code{"(I) Distributional Properties"}).
#' @param exclude_categories Character vector of categories to exclude from the ordination.
#'   Default is \code{c("(VI) Trajectory & Lineage Dynamics", "(VIII) Computational Scalability")}
#'   because these two categories contain only two measures each, leaving the 6 canonical categories.
#' @param metrics Optional character vector of specific metrics to include.
#' @param normalize_scores Logical. If \code{TRUE} (default), scores are direction-aware normalized to [0, 1]
#'   so that metrics on disparate scales (e.g. MB vs distance) do not artificially distort Euclidean distance.
#' @param palette Optional named character vector of colors. For \code{ordination_by = "summaries"},
#'   maps \code{"gene"}, \code{"cell"}, and \code{"global"}.
#' @param title Optional plot title.
#' @param base_size Numeric. Base font size. Default \code{11}.
#' @param ncol Integer. Number of columns when \code{by_category = TRUE}. Default \code{2}.
#'
#' @return A \code{ggplot} or \code{patchwork} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_metric_mds <- function(
  benchmark_data,
  ordination_by      = c("summaries", "simulators"),
  by_category        = FALSE,
  as_list            = FALSE,
  category           = NULL,
  exclude_categories = c("(VI) Trajectory & Lineage Dynamics", "(VIII) Computational Scalability"),
  metrics            = NULL,
  normalize_scores   = TRUE,
  palette            = NULL,
  title              = NULL,
  base_size          = 11,
  ncol               = 2
) {
  ordination_by <- match.arg(ordination_by)
  df <- .ingest_bubble_data(benchmark_data)

  if (!is.null(category) && !identical(category, "all")) {
    df <- df[df$Category %in% category, , drop = FALSE]
  } else if (!is.null(exclude_categories)) {
    ex_pattern <- paste(exclude_categories, collapse = "|")
    df <- df[!grepl(ex_pattern, df$Category, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(metrics)) {
    df <- df[df$Metric %in% metrics, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No data remaining after applying category/metric filters.")

  if (normalize_scores) {
    if ("Score_Norm" %in% colnames(df)) {
      score_col <- "Score_Norm"
    } else {
      df <- .normalize_bubble_scores(df)
      score_col <- "Score_Norm"
    }
  } else {
    score_col <- if ("Score_Raw" %in% colnames(df)) "Score_Raw" else "Score"
  }

  type_colors <- c(
    "gene"   = "#E04038",  # red (gene-level)
    "cell"   = "#3876E0",  # blue (cell-level)
    "global" = "#28A745"   # green (global structure / mixing)
  )
  if (!is.null(palette) && ordination_by == "summaries") {
    for (nm in names(palette)) {
      if (nm %in% names(type_colors)) type_colors[nm] <- palette[nm]
    }
  }

  if (ordination_by == "summaries") {
    # Helper to build a single summary ordination panel
    .build_summary_panel <- function(sub_df, panel_title = NULL, panel_tag = NULL, show_legend = TRUE, p_base_size = base_size) {
      wide_m <- stats::reshape(
        sub_df[, c("Method", "Metric", score_col)],
        timevar   = "Method",
        idvar     = "Metric",
        direction = "wide"
      )
      rnames <- wide_m$Metric
      mat <- as.matrix(wide_m[, -1, drop = FALSE])
      rownames(mat) <- rnames
      mat[is.na(mat)] <- 0

      if (nrow(mat) < 3) return(NULL)

      if (!normalize_scores) {
        mat <- t(scale(t(mat)))
        mat[is.na(mat)] <- 0
      }

      dist_mat <- stats::dist(mat)
      mds      <- stats::cmdscale(dist_mat, k = 2, eig = TRUE)

      mds_df <- data.frame(
        Metric = rownames(mat),
        MDS1   = mds$points[, 1],
        MDS2   = mds$points[, 2],
        stringsAsFactors = FALSE
      )

      # Classify into gene, cell, global
      mds_df$Type <- .SUMMARY_TYPE_MAP[mds_df$Metric]
      for (i in seq_len(nrow(mds_df))) {
        if (is.na(mds_df$Type[i])) {
          m <- tolower(mds_df$Metric[i])
          if (grepl("gene|feat|deg|cpm|var|disp|fc", m)) {
            mds_df$Type[i] <- "gene"
          } else if (grepl("cell|lib|ldf|knn", m)) {
            mds_df$Type[i] <- "cell"
          } else {
            mds_df$Type[i] <- "global"
          }
        }
      }
      present_types <- intersect(c("gene", "cell", "global"), unique(as.character(mds_df$Type)))
      mds_df$Type   <- factor(mds_df$Type, levels = present_types)

      # Display labels
      mds_df$Display_Label <- .SUMMARY_DISPLAY_NAME_MAP[mds_df$Metric]
      idx_na <- is.na(mds_df$Display_Label)
      if (any(idx_na)) {
        mds_df$Display_Label[idx_na] <- gsub("_", " ", mds_df$Metric[idx_na])
      }

      p_panel <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$MDS1, y = .data$MDS2)) +
        ggplot2::geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
        ggplot2::geom_point(
          ggplot2::aes(fill = .data$Type, color = .data$Type),
          shape = 21, size = ifelse(by_category, 3.2, 3.8), stroke = 0.75, alpha = 0.88
        )

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p_panel <- p_panel + ggrepel::geom_text_repel(
          ggplot2::aes(label = .data$Display_Label, color = .data$Type),
          fontface = "bold", size = p_base_size * 0.28,
          segment.color = "#999999", segment.size = 0.32,
          box.padding = 0.28, point.padding = 0.25,
          max.overlaps = 60, show.legend = FALSE
        )
      } else {
        p_panel <- p_panel + ggplot2::geom_text(
          ggplot2::aes(label = .data$Display_Label, color = .data$Type),
          vjust = -0.8, fontface = "bold", size = p_base_size * 0.28, show.legend = FALSE
        )
      }

      p_panel <- p_panel +
        ggplot2::scale_color_manual(values = type_colors[present_types], name = NULL, drop = TRUE) +
        ggplot2::scale_fill_manual(values = type_colors[present_types], name = NULL, drop = TRUE) +
        ggplot2::labs(
          title = panel_title,
          tag   = panel_tag,
          x     = "MDS dim. 1",
          y     = "MDS dim. 2"
        ) +
        ggplot2::theme_bw(base_size = p_base_size) +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_line(color = "#EFEFEF", linewidth = 0.5),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.title       = ggplot2::element_text(face = "bold", size = p_base_size * 0.95),
          plot.tag         = ggplot2::element_text(face = "bold", size = p_base_size * 1.15),
          axis.title       = ggplot2::element_text(face = "bold", size = p_base_size * 0.90),
          axis.text        = ggplot2::element_text(size = p_base_size * 0.80, color = "black"),
          legend.position  = if (show_legend) "bottom" else "none",
          legend.direction = "horizontal",
          legend.text      = ggplot2::element_text(face = "bold", size = p_base_size * 0.85),
          legend.key       = ggplot2::element_blank(),
          plot.margin      = ggplot2::margin(8, 8, 8, 8)
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.8, shape = 21, stroke = 0.8)))

      p_panel
    }

    if (by_category) {
      ordered_cats <- .CANONICAL_CATEGORY_ORDER[.CANONICAL_CATEGORY_ORDER %in% unique(df$Category)]
      if (length(ordered_cats) == 0) ordered_cats <- unique(df$Category)

      plot_list <- list()
      tags <- letters[seq_along(ordered_cats)]
      for (idx in seq_along(ordered_cats)) {
        c_nm <- ordered_cats[idx]
        sub_c <- df[df$Category == c_nm, , drop = FALSE]
        if (length(unique(sub_c$Metric)) >= 3) {
          p_c <- .build_summary_panel(
            sub_c,
            panel_title = if (as_list) sprintf("Metric MDS Plot - %s", c_nm) else c_nm,
            panel_tag   = if (as_list) NULL else tags[idx],
            show_legend = if (as_list) TRUE else (length(plot_list) == 0),
            p_base_size = base_size * 0.9
          )
          if (!is.null(p_c)) {
            plot_list[[length(plot_list) + 1]] <- p_c
          }
        }
      }

      if (length(plot_list) == 0) {
        stop("None of the evaluated categories have >= 3 distinct metrics for MDS ordination.")
      }

      valid_cats <- character(0)
      for (idx in seq_along(ordered_cats)) {
        c_nm <- ordered_cats[idx]
        sub_c <- df[df$Category == c_nm, , drop = FALSE]
        if (length(unique(sub_c$Metric)) >= 3) valid_cats <- c(valid_cats, c_nm)
      }
      names(plot_list) <- valid_cats

      if (as_list) {
        return(plot_list)
      }

      if (requireNamespace("patchwork", quietly = TRUE)) {
        p_final <- patchwork::wrap_plots(plot_list, ncol = ncol, guides = "collect") &
          ggplot2::theme(legend.position = "bottom")
        if (!is.null(title)) {
          p_final <- p_final + patchwork::plot_annotation(
            title = title,
            theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.15))
          )
        }
        return(p_final)
      } else {
        return(plot_list[[1]])
      }
    } else {
      default_title <- if (!is.null(title)) {
        title
      } else if (!is.null(category) && length(category) == 1) {
        sprintf("Metric MDS Plot - %s", category)
      } else if (length(unique(df$Category)) == 1) {
        sprintf("Metric MDS Plot - %s", unique(df$Category))
      } else {
        "Metric MDS Plot"
      }
      p <- .build_summary_panel(
        df,
        panel_title = default_title,
        panel_tag   = NULL,
        show_legend = TRUE,
        p_base_size = base_size
      )
      if (is.null(p)) stop("MDS ordination of summaries requires at least 3 distinct metrics.")
      return(p)
    }

  } else {
    # -------------------------------------------------------------------------
    # Simulators Ordination (Backwards compatible)
    # -------------------------------------------------------------------------
    wide_m <- stats::reshape(
      df[, c("Method", "Metric", score_col)],
      timevar   = "Metric",
      idvar     = "Method",
      direction = "wide"
    )
    rnames <- wide_m$Method
    mat <- as.matrix(wide_m[, -1, drop = FALSE])
    rownames(mat) <- rnames
    colnames(mat) <- gsub(paste0("^", score_col, "\\."), "", colnames(mat))
    mat[is.na(mat)] <- 0

    if (nrow(mat) < 3) stop("MDS ordination requires at least 3 distinct simulators.")

    dist_mat <- stats::dist(mat)
    mds      <- stats::cmdscale(dist_mat, k = 2, eig = TRUE)
    var_exp  <- round(100 * (mds$eig[seq_len(2)] / sum(abs(mds$eig))), 1)

    mds_df <- data.frame(
      Simulator = rownames(mat),
      MDS1      = mds$points[, 1],
      MDS2      = mds$points[, 2],
      stringsAsFactors = FALSE
    )

    default_sim_cols <- c(
      "scDesign3" = "#1E8449",
      "Splatter"  = "#2E86AB",
      "SCRIP"     = "#E67E22",
      "SymSim"    = "#8E44AD",
      "dyngen"    = "#D4AC0D",
      "simATAC"   = "#C0392B"
    )
    sim_methods <- unique(mds_df$Simulator)
    if (is.null(palette)) {
      active_pal <- default_sim_cols[sim_methods]
      missing_sims <- sim_methods[is.na(active_pal)]
      if (length(missing_sims) > 0) {
        extra_cols <- grDevices::hcl.colors(length(missing_sims), palette = "Dark 3")
        names(extra_cols) <- missing_sims
        active_pal[missing_sims] <- extra_cols
      }
    } else {
      active_pal <- palette
    }

    p <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$MDS1, y = .data$MDS2, color = .data$Simulator, label = .data$Simulator)) +
      ggplot2::geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
      ggplot2::geom_point(size = 4.8, alpha = 0.9)

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        fill = "white", fontface = "bold", size = base_size * 0.30,
        label.padding = ggplot2::unit(0.22, "lines"), label.r = ggplot2::unit(0.35, "lines"),
        label.size = 0.7, box.padding = 0.4, show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(vjust = -1.0, fontface = "bold", size = base_size * 0.32, show.legend = FALSE)
    }

    p <- p +
      ggplot2::scale_color_manual(values = active_pal, guide = "none") +
      ggplot2::labs(
        title = title,
        x     = paste0("MDS dim. 1 (", var_exp[1], "% var)"),
        y     = paste0("MDS dim. 2 (", var_exp[2], "% var)")
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "#EFEFEF", linewidth = 0.5),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title       = ggplot2::element_text(face = "bold", size = base_size * 0.95),
        axis.text        = ggplot2::element_text(size = base_size * 0.85, color = "black"),
        plot.margin      = ggplot2::margin(10, 10, 10, 10)
      )
    return(p)
  }
}

# 6c. plot_metric_pca()
# ============================================================================

#' Principal Component Analysis (PCA) Dashboard of Benchmark Metrics and Methods
#'
#' Generates a dual-panel PCA dashboard for evaluating single-cell simulation benchmark metrics
#' and methods.
#'
#' All metrics present in \code{benchmark_data} (e.g. all 62 canonical metrics across categories)
#' are mathematically considered in the PCA decomposition. Panel \code{"both"} displays:
#' \itemize{
#'   \item \strong{Panel a (Methods Ordination)}: Simulators mapped into PC1 vs PC2 coordinates with
#'     solid black crosshairs, centroid markers, and white rounded badge labels with colored borders.
#'   \item \strong{Panel b (PC Loadings Vectors)}: Directional loading vectors (arrows) radiating from
#'     the origin \code{(0, 0)} to loading coordinates \code{(w1, w2)}, colored by summary type
#'     (\code{"gene"} = red, \code{"cell"} = blue, \code{"global"} = green) with rounded badge labels.
#'     Pass \code{top_n_loadings = NULL} or \code{62} to render all metric loading vectors.
#' }
#'
#' @param benchmark_data Benchmark summary data.frame, matrix, or list.
#' @param panel Character. Which panel(s) to render: \code{"both"} (default; 2-panel stacked layout),
#'   \code{"methods"} (panel a only), or \code{"loadings"} (panel b only).
#' @param by_category Logical. If \code{TRUE}, perform separate PCA ordinations per evaluation category. Default \code{FALSE}.
#' @param as_list Logical. If \code{TRUE} and \code{by_category = TRUE}, returns a named list of ggplot objects. Default \code{FALSE}.
#' @param category Optional character vector to filter by evaluation category.
#' @param exclude_categories Optional character vector of categories to exclude (e.g. scalability).
#' @param metrics Optional character vector of specific metrics to include.
#' @param palette Optional named character vector of colors for simulators.
#' @param top_n_loadings Integer or \code{NULL}. Maximum number of loading vectors to plot in panel b
#'   (ordered by vector magnitude in PC1-PC2 space). Default is \code{14} to prevent visual clutter;
#'   set to \code{NULL} or \code{62} to display all evaluated metrics.
#' @param base_size Numeric. Base font size. Default \code{11}.
#'
#' @return A \code{ggplot} or \code{patchwork} object representing the ordination figure.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_metric_pca <- function(
  benchmark_data,
  panel              = c("both", "methods", "loadings"),
  by_category        = FALSE,
  as_list            = FALSE,
  category           = NULL,
  exclude_categories = NULL,
  metrics            = NULL,
  palette            = NULL,
  top_n_loadings     = 14,
  base_size          = 11
) {
  panel <- match.arg(panel)
  df <- .ingest_bubble_data(benchmark_data)

  if (by_category) {
    ordered_cats <- .CANONICAL_CATEGORY_ORDER[.CANONICAL_CATEGORY_ORDER %in% unique(df$Category)]
    if (length(ordered_cats) == 0) ordered_cats <- unique(df$Category)
    if (!is.null(exclude_categories)) {
      ex_pattern <- paste(exclude_categories, collapse = "|")
      ordered_cats <- ordered_cats[!grepl(ex_pattern, ordered_cats, ignore.case = TRUE)]
    } else {
      ordered_cats <- ordered_cats[!grepl("Trajectory|Scalability", ordered_cats, ignore.case = TRUE)]
    }

    plot_list <- list()
    for (cn in ordered_cats) {
      sub_c <- df[df$Category == cn, , drop = FALSE]
      if (length(unique(sub_c$Metric)) >= 2) {
        p_c <- plot_metric_pca(
          benchmark_data     = sub_c,
          panel              = panel,
          by_category        = FALSE,
          as_list            = FALSE,
          category           = NULL,
          exclude_categories = NULL,
          metrics            = NULL,
          palette            = palette,
          top_n_loadings     = if (is.null(top_n_loadings)) NULL else top_n_loadings,
          base_size          = base_size
        )
        plot_list[[cn]] <- p_c
      }
    }
    return(plot_list)
  }

  if (!is.null(category) && !identical(category, "all")) {
    df <- df[df$Category %in% category, , drop = FALSE]
  } else if (!is.null(exclude_categories)) {
    ex_pattern <- paste(exclude_categories, collapse = "|")
    df <- df[!grepl(ex_pattern, df$Category, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(metrics)) {
    df <- df[df$Metric %in% metrics, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No data remaining after applying category/metric filters.")

  score_col <- if ("Score_Norm" %in% colnames(df)) "Score_Norm" else
               if ("Score_Raw"  %in% colnames(df)) "Score_Raw"  else "Score"

  # Reshape: Simulators (rows) x Metrics (columns)
  wide_m <- stats::reshape(
    df[, c("Method", "Metric", score_col)],
    timevar   = "Metric",
    idvar     = "Method",
    direction = "wide"
  )
  rnames <- wide_m$Method
  mat <- as.matrix(wide_m[, -1, drop = FALSE])
  rownames(mat) <- rnames
  colnames(mat) <- gsub(paste0("^", score_col, "\\."), "", colnames(mat))
  mat[is.na(mat)] <- 0

  # Remove zero-variance columns for PCA
  col_vars <- apply(mat, 2, stats::var)
  mat <- mat[, col_vars > 1e-12, drop = FALSE]

  if (nrow(mat) < 3 || ncol(mat) < 2) stop("PCA requires at least 3 simulators and 2 non-constant metrics.")

  pca     <- stats::prcomp(mat, scale. = TRUE, center = TRUE)
  pca_var <- round(100 * summary(pca)$importance[2, seq_len(2)], 1)

  # Scores data frame (Methods)
  scores_df <- data.frame(
    Method = rownames(pca$x),
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  method_names <- unique(scores_df$Method)
  n_m <- length(method_names)
  base_pal <- c("#3B82F6", "#10B981", "#EF4444", "#F59E0B", "#8B5CF6", "#EC4899", "#14B8A6", "#6366F1", "#84CC16", "#06B6D4")
  if (!is.null(palette)) {
    method_cols <- palette
  } else if (n_m <= length(base_pal)) {
    method_cols <- stats::setNames(base_pal[seq_len(n_m)], method_names)
  } else {
    method_cols <- stats::setNames(grDevices::hcl.colors(n_m, "Dark 3"), method_names)
  }

  type_colors <- c(
    "gene"   = "#E04038",  # red
    "cell"   = "#3876E0",  # blue
    "global" = "#28A745"   # green
  )

  # ---------------------------------------------------------------------------
  # Panel a: Methods Ordination
  # ---------------------------------------------------------------------------
  p_a <- ggplot2::ggplot(scores_df, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Method, fill = .data$Method),
      shape = 21, size = 4.8, stroke = 0.9, alpha = 0.92
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_a <- p_a + ggrepel::geom_label_repel(
      ggplot2::aes(label = .data$Method, color = .data$Method),
      fill = "white", fontface = "bold", size = base_size * 0.30,
      label.padding = ggplot2::unit(0.22, "lines"), label.r = ggplot2::unit(0.35, "lines"),
      label.size = 0.7, box.padding = 0.4, point.padding = 0.3,
      max.overlaps = 50, show.legend = FALSE
    )
  } else {
    p_a <- p_a + ggplot2::geom_text(
      ggplot2::aes(label = .data$Method, color = .data$Method),
      vjust = -0.8, fontface = "bold", size = base_size * 0.30, show.legend = FALSE
    )
  }

  p_a <- p_a +
    ggplot2::scale_color_manual(values = method_cols, name = "method") +
    ggplot2::scale_fill_manual(values = method_cols, name = "method") +
    ggplot2::labs(
      title    = "Simulator Performance Biplot (PC1 vs. PC2)",
      subtitle = sprintf("Variance Explained: PC1 = %s%%, PC2 = %s%% (Cumulative = %s%%)",
                         pca_var[1], pca_var[2], round(pca_var[1] + pca_var[2], 1)),
      x = sprintf("PC1 (%s%%)", pca_var[1]),
      y = sprintf("PC2 (%s%%)", pca_var[2])
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#EFEFEF", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.05),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.85, color = "#475569"),
      axis.title       = ggplot2::element_text(face = "bold", size = base_size * 0.95),
      axis.text        = ggplot2::element_text(size = base_size * 0.85, color = "black"),
      legend.title     = ggplot2::element_text(face = "bold", size = base_size * 0.90),
      legend.text      = ggplot2::element_text(size = base_size * 0.82),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 4, shape = 21, fill = method_cols)),
      fill  = "none"
    )

  if (panel == "methods") return(p_a)

  # ---------------------------------------------------------------------------
  # Panel b: PC Loadings Vectors
  # ---------------------------------------------------------------------------
  loadings_df <- as.data.frame(pca$rotation[, seq_len(2), drop = FALSE])
  loadings_df$Metric <- rownames(loadings_df)
  loadings_df$mag <- sqrt(loadings_df$PC1^2 + loadings_df$PC2^2)

  if (!is.null(top_n_loadings) && is.finite(top_n_loadings) && top_n_loadings < nrow(loadings_df)) {
    loadings_df <- head(loadings_df[order(loadings_df$mag, decreasing = TRUE), ], top_n_loadings)
  }

  # Assign type: gene, cell, global
  loadings_df$Type <- .SUMMARY_TYPE_MAP[loadings_df$Metric]
  for (i in seq_len(nrow(loadings_df))) {
    if (is.na(loadings_df$Type[i])) {
      m <- tolower(loadings_df$Metric[i])
      if (grepl("gene|feat|deg|cpm|var|disp|fc", m)) {
        loadings_df$Type[i] <- "gene"
      } else if (grepl("cell|lib|ldf|knn", m)) {
        loadings_df$Type[i] <- "cell"
      } else {
        loadings_df$Type[i] <- "global"
      }
    }
  }
  present_types <- intersect(c("gene", "cell", "global"), unique(as.character(loadings_df$Type)))
  loadings_df$Type <- factor(loadings_df$Type, levels = present_types)

  loadings_df$Display_Label <- .SUMMARY_DISPLAY_NAME_MAP[loadings_df$Metric]
  idx_na <- is.na(loadings_df$Display_Label)
  if (any(idx_na)) {
    loadings_df$Display_Label[idx_na] <- gsub("_", " ", loadings_df$Metric[idx_na])
  }

  p_b <- ggplot2::ggplot(loadings_df) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2, color = .data$Type),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.24, "cm"), type = "closed"),
      linewidth = 0.72, alpha = 0.85
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_b <- p_b + ggrepel::geom_label_repel(
      ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$Display_Label, color = .data$Type),
      fill = "white", fontface = "bold", size = base_size * 0.25,
      label.padding = ggplot2::unit(0.18, "lines"), label.r = ggplot2::unit(0.30, "lines"),
      label.size = 0.65, box.padding = 0.35, point.padding = 0.25,
      max.overlaps = 70, show.legend = FALSE
    )
  } else {
    p_b <- p_b + ggplot2::geom_text(
      ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$Display_Label, color = .data$Type),
      vjust = -0.6, fontface = "bold", size = base_size * 0.25, show.legend = FALSE
    )
  }

  p_b <- p_b +
    ggplot2::scale_color_manual(values = type_colors[present_types], name = "summary type", drop = TRUE) +
    ggplot2::labs(
      title    = sprintf("Metric Loadings Vectors (Top %d Discriminating Metrics)", nrow(loadings_df)),
      subtitle = "Vector length and direction indicate metric influence on PC1 and PC2 separation",
      x = sprintf("PC1 (%s%%)", pca_var[1]),
      y = sprintf("PC2 (%s%%)", pca_var[2])
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#EFEFEF", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.05),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.85, color = "#475569"),
      axis.title       = ggplot2::element_text(face = "bold", size = base_size * 0.95),
      axis.text        = ggplot2::element_text(size = base_size * 0.85, color = "black"),
      legend.title     = ggplot2::element_text(face = "bold", size = base_size * 0.90),
      legend.text      = ggplot2::element_text(size = base_size * 0.82),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.0, shape = NA))
    )

  if (panel == "loadings") return(p_b)

  # ---------------------------------------------------------------------------
  # Stacked layout (both panel a and panel b)
  # ---------------------------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_final <- patchwork::wrap_plots(p_a, p_b, ncol = 1) +
      patchwork::plot_annotation(tag_levels = "a") &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = base_size * 1.1))
    return(p_final)
  } else {
    return(p_a)
  }
}

# ============================================================================
# 8. plot_benchmark_bubble_matrix()
# ============================================================================

#' Plot Multi-Dimensional Benchmarking Bubble Matrix
#'
#' Produces a multi-dimensional benchmarking bubble matrix comparing simulation methods.
#' Methods appear as rows; evaluation metrics appear as columns grouped under the eight
#' canonical evaluation categories (I-VIII) displayed as colored header strips at the top.
#'
#' \if{html}{\figure{benchmark_bubble_matrix.png}{options: width="800" alt="Benchmarking Bubble Matrix"}}
#'
#' \strong{Bubble encoding:}
#' \itemize{
#'   \item \strong{Size}: Normalized fidelity score [0, 1] -- larger = better performance.
#'   \item \strong{Fill color}: Evaluation category (one of the eight canonical groups).
#'   \item \strong{Circle (shape 21)}: Regular performance (score < 0.96).
#'   \item \strong{Square (shape 22)}: Top performer (score >= 0.96).
#'   \item \strong{Small grey dot}: Metric not available/computed for that method.
#' }
#'
#' \strong{Score normalization} is direction-aware: distance/error metrics are inverted
#' so that 1 always means best performance, while concordance/correlation metrics are used
#' directly (higher = better).
#'
#' @param data Input data. One of:
#'   \itemize{
#'     \item An \code{scSimEval_consolidated} object from \code{\link{evaluate_multiple_datasets}()}.
#'     \item A named list of benchmark outputs from \code{\link{evaluate_simulation_accuracy}()}.
#'     \item A tidy \code{data.frame} with columns: \code{Method}, \code{Metric}, \code{Category},
#'           and \code{Score} (or \code{Value}).
#'   }
#' @param method_classes Optional named list mapping method names to row group labels.
#'   E.g. \code{list("scRNA-seq" = c("Splatter", "SymSim"), "Multiomics" = c("scDesign3", "dyngen"))}.
#' @param category_colors Optional named character vector to override the canonical palette.
#' @param metrics_order Optional character vector for explicit metric ordering (left to right).
#' @param methods_order Optional character vector for explicit method ordering (top to bottom).
#' @param rank_methods Logical. Automatically rank methods along the vertical axis by mean overall fidelity score. Default \code{TRUE}.
#' @param compact_strips Logical. Use compact two-line facet strip titles for clean horizontal layout. Default \code{TRUE}.
#' @param normalize_scores Logical. Normalize raw distances into [0, 1] fidelity scores. Default \code{TRUE}.
#' @param show_missing_dots Logical. Show grey dots for missing metric combinations. Default \code{TRUE}.
#' @param bubble_size_range Numeric vector of length 2: min and max bubble size. Default \code{c(2.0, 7.8)}.
#' @param title Character. Main plot title. If \code{NULL}, dynamically derived from metric count.
#' @param subtitle Character. Plot subtitle.
#' @param base_size Numeric. Base font size. Default \code{11}.
#' @param show_score_labels Logical. Annotate bubbles with numeric score. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object rendering the multi-dimensional bubble matrix.
#' @export
#' @examples
#' # Tidy benchmark data frame of standardized scores
#' bm_df <- data.frame(
#'   Method   = rep(c("Splat", "scDesign3", "SymSim", "dyngen"), each = 4),
#'   Category = rep(c("Accuracy", "Accuracy", "Cellular Structure", "Scalability"), 4),
#'   Metric   = rep(c("KS Distance", "Wasserstein", "Silhouette ASW", "CPU Time"), 4),
#'   Score    = c(0.8, 0.7, 0.9, 0.6, 0.5, 0.8, 0.6, 0.5, 0.8, 0.7, 0.6, 0.5, 0.7, 0.6, 0.9, 0.5)
#' )
#' p <- plot_benchmark_bubble_matrix(bm_df)
#' p
#' \donttest{
#' # Load benchmark summary across simulators
#' demo_file <- system.file("shiny/scSimEvalApp/data/demo_benchmark_data.rds", package = "scSimEval")
#' if (file.exists(demo_file)) {
#'   demo <- readRDS(demo_file)
#'   p <- plot_benchmark_bubble_matrix(
#'     data = demo$benchmark_summary_table,
#'     title = "Single-Cell Multiomics Simulation Benchmark",
#'     method_classes = list(
#'       "scRNA-seq"  = c("Splatter", "SymSim"),
#'       "scATAC-seq" = c("simATAC", "SCRIP"),
#'       "Multiomics" = c("scDesign3", "dyngen")
#'     )
#'   )
#'   print(p)
#' }
#' }
plot_benchmark_bubble_matrix <- function(
  data,
  method_classes    = NULL,
  category_colors   = NULL,
  metrics_order     = NULL,
  methods_order     = NULL,
  rank_methods      = TRUE,
  compact_strips    = TRUE,
  normalize_scores  = TRUE,
  show_missing_dots = TRUE,
  bubble_size_range = c(2.0, 7.8),
  title             = NULL,
  subtitle          = "Standardized Direction-Aware Fidelity Scores [0, 1] Across 8 Canonical Evaluation Categories",
  base_size         = 11,
  show_score_labels = FALSE
) {
  # --- 1. Ingest ---------------------------------------------------------
  df <- .ingest_bubble_data(data)

  # --- 2. Category remapping ---------------------------------------------
  if ("Category" %in% colnames(df)) {
    df$Display_Category <- ifelse(df$Category %in% names(.LEGACY_CATEGORY_MAP),
                                   .LEGACY_CATEGORY_MAP[df$Category], df$Category)
  } else {
    df$Display_Category <- "(I) Distributional Properties"
  }

  # Metric-level category override (highest priority)
  if ("Metric" %in% colnames(df)) {
    mc_override <- .METRIC_CATEGORY_MAP[df$Metric]
    idx         <- !is.na(mc_override)
    df$Display_Category[idx] <- mc_override[idx]
  }

  # --- 3. Display metric labels ------------------------------------------
  if ("Metric" %in% colnames(df)) {
    lbl_map <- .METRIC_DISPLAY_LABELS[df$Metric]
    df$Display_Metric <- ifelse(!is.na(lbl_map), lbl_map, df$Metric)
  } else {
    df$Display_Metric <- as.character(df$Metric)
  }

  # --- 4. Score normalization --------------------------------------------
  if (normalize_scores && "Score_Raw" %in% colnames(df)) {
    df <- .normalize_bubble_scores(df)
  } else if (!"Score_Norm" %in% colnames(df)) {
    df$Score_Norm <- pmax(0, pmin(1, df$Score_Raw))
  }

  # --- 5. Method class ---------------------------------------------------
  if (!is.null(method_classes)) {
    df$Method_Class <- NA_character_
    for (cls in names(method_classes))
      df$Method_Class[df$Method %in% method_classes[[cls]]] <- cls
    df$Method_Class[is.na(df$Method_Class)] <- "Other Methods"
  } else {
    df$Method_Class <- "All Methods"
  }

  # --- 6. Factor ordering ------------------------------------------------
  all_methods <- unique(df$Method)
  if (!is.null(methods_order)) {
    ordered_m <- c(methods_order[methods_order %in% all_methods], setdiff(all_methods, methods_order))
  } else if (rank_methods && "Score_Norm" %in% colnames(df)) {
    ranks_df <- stats::aggregate(Score_Norm ~ Method, data = df, FUN = mean, na.rm = TRUE)
    ordered_m <- ranks_df$Method[order(ranks_df$Score_Norm, decreasing = TRUE)]
    ordered_m <- c(ordered_m, setdiff(all_methods, ordered_m))
  } else {
    ordered_m <- all_methods
  }
  # rev(ordered_m) places the top-performing / first simulator at the top of the Y-axis
  df$Method <- factor(df$Method, levels = rev(ordered_m))

  all_cats_present <- unique(as.character(df$Display_Category))
  if (compact_strips) {
    strip_map <- .COMPACT_CATEGORY_STRIPS
    df$Display_Category <- as.character(df$Display_Category)
    df$Display_Category <- ifelse(df$Display_Category %in% names(strip_map),
                                  strip_map[df$Display_Category], df$Display_Category)
    cat_order <- unname(strip_map[.CANONICAL_CATEGORY_ORDER[.CANONICAL_CATEGORY_ORDER %in% names(strip_map)]])
    cat_order <- c(cat_order, setdiff(unique(df$Display_Category), cat_order))
    df$Display_Category <- factor(df$Display_Category, levels = cat_order)
  } else {
    cat_order <- c(.CANONICAL_CATEGORY_ORDER[.CANONICAL_CATEGORY_ORDER %in% all_cats_present],
                   setdiff(all_cats_present, .CANONICAL_CATEGORY_ORDER))
    df$Display_Category <- factor(df$Display_Category, levels = cat_order)
  }

  # Metric ordering: custom or alphabetical within category
  if (!is.null(metrics_order)) {
    all_mets <- unique(df$Display_Metric)
    ord_mets <- c(metrics_order[metrics_order %in% all_mets], setdiff(all_mets, metrics_order))
    df$Display_Metric <- factor(df$Display_Metric, levels = ord_mets)
  } else {
    category_metric_order <- character(0)
    for (cat in levels(df$Display_Category)) {
      mets_in <- sort(unique(df$Display_Metric[df$Display_Category == cat]))
      category_metric_order <- c(category_metric_order, mets_in)
    }
    df$Display_Metric <- factor(df$Display_Metric, levels = category_metric_order)
  }

  if (!is.null(method_classes)) {
    cls_lvls <- names(method_classes)
    has_other <- "Other Methods" %in% df$Method_Class
    df$Method_Class <- factor(df$Method_Class,
                               levels = c(cls_lvls[cls_lvls %in% df$Method_Class],
                                          if (has_other) "Other Methods"))
  }

  df$Is_Top <- df$Score_Norm >= 0.96

  # --- 7. Category palette -----------------------------------------------
  base_pal <- .CANONICAL_CATEGORY_COLORS
  if (compact_strips) {
    active_palette <- setNames(
      base_pal[names(.COMPACT_CATEGORY_STRIPS)],
      unname(.COMPACT_CATEGORY_STRIPS)
    )
    fill_labels <- names(.COMPACT_CATEGORY_STRIPS)
    names(fill_labels) <- unname(.COMPACT_CATEGORY_STRIPS)
  } else {
    active_palette <- base_pal
    fill_labels <- ggplot2::waiver()
  }

  if (!is.null(category_colors))
    for (k in names(category_colors)) active_palette[k] <- category_colors[[k]]
  cats_present <- levels(df$Display_Category)
  active_palette <- active_palette[names(active_palette) %in% cats_present]
  unmapped <- setdiff(cats_present, names(active_palette))
  if (length(unmapped) > 0) {
    greys <- grDevices::gray.colors(length(unmapped), start = 0.45, end = 0.75)
    names(greys) <- unmapped
    active_palette <- c(active_palette, greys)
  }

  # --- 8. Build ggplot ---------------------------------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Display_Metric, y = .data$Method))

  # Subtle horizontal row guides
  p <- p + ggplot2::geom_hline(
    yintercept = seq_along(levels(df$Method)),
    color = "#F1F5F9", linewidth = 0.5
  )

  # Background grey grid dots for missing metric combinations
  if (show_missing_dots) {
    grid_df <- expand.grid(Method = levels(df$Method), Display_Metric = levels(df$Display_Metric),
                            stringsAsFactors = FALSE)
    mc_map  <- unique(df[, c("Display_Metric", "Display_Category")])
    grid_df <- merge(grid_df, mc_map, by = "Display_Metric", all.x = TRUE)
    p <- p + ggplot2::geom_point(
      data = grid_df, ggplot2::aes(x = .data$Display_Metric, y = .data$Method),
      shape = 16, size = 0.85, color = "#D5D8DC", alpha = 0.60, inherit.aes = FALSE
    )
  }

  # Regular bubbles (circles)
  df_reg <- df[!df$Is_Top, ]
  if (nrow(df_reg) > 0)
    p <- p + ggplot2::geom_point(
      data = df_reg,
      ggplot2::aes(size = .data$Score_Norm, fill = .data$Display_Category),
      shape = 21, color = "#1E293B", stroke = 0.5, alpha = 0.92
    )

  # Top performer bubbles (squares)
  df_top <- df[df$Is_Top, ]
  if (nrow(df_top) > 0)
    p <- p + ggplot2::geom_point(
      data = df_top,
      ggplot2::aes(size = .data$Score_Norm, fill = .data$Display_Category),
      shape = 22, color = "#0F172A", stroke = 0.90, alpha = 0.98
    )

  if (show_score_labels)
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$Score_Norm)),
      size = 2.1, color = "#1A252F", fontface = "bold"
    )

  # --- 9. Scales ---------------------------------------------------------
  p <- p +
    ggplot2::scale_fill_manual(
      values = active_palette,
      labels = fill_labels,
      name   = "Evaluation Category (Full Provenance Mapping)",
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 5.0, shape = 21), nrow = 2,
        title.position = "top", order = 1
      )
    ) +
    ggplot2::scale_size_continuous(
      range  = bubble_size_range, limits = c(0, 1),
      breaks = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0),
      labels = c("0.10", "0.30", "0.50", "0.70", "0.90", "1.00 (Best)"),
      name   = "Fidelity Score (Larger Bubble = Higher Fidelity)",
      guide  = ggplot2::guide_legend(
        override.aes = list(fill = "#94A3B8", shape = 21, color = "#1E293B"),
        nrow = 1, title.position = "top", order = 2
      )
    )

  # --- 10. Faceting ------------------------------------------------------
  has_multi_class <- length(unique(df$Method_Class)) > 1
  has_cats        <- length(levels(df$Display_Category)) >= 1

  if (has_multi_class && has_cats) {
    p <- p + ggplot2::facet_grid(Method_Class ~ Display_Category,
                                  scales = "free", space = "free")
  } else if (has_cats) {
    p <- p + ggplot2::facet_grid(~ Display_Category,
                                  scales = "free_x", space = "free_x")
  } else if (has_multi_class) {
    p <- p + ggplot2::facet_grid(Method_Class ~ .,
                                  scales = "free_y", space = "free_y")
  }

  # --- 11. Publication theme ---------------------------------------------
  if (is.null(title)) {
    n_metrics <- length(unique(df$Display_Metric))
    title <- sprintf("Benchmarking Single-Cell Simulation Fidelity Across %d Canonical Measures", n_metrics)
  }
  y_title <- if (rank_methods && is.null(methods_order)) "Ranked Simulators" else NULL

  p <- p +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = NULL,
      y        = y_title,
      caption  = paste0(
        "Circle: standard performance (< 0.96)  |  Square: top performer (>= 0.96)",
        if (show_missing_dots) "  |  Grey dot: metric not computed." else ".\n",
        if (show_missing_dots) "\n" else "",
        "All ", length(unique(df$Display_Metric)), " metrics direction-normalized: for error/distance metrics, scores are inverted as 1 - norm(x) so 1.0 always indicates closest agreement to empirical reference."
      )
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1,
                                                 face = "bold", size = base_size * 0.77,
                                                 color = "#1E293B"),
      axis.text.y        = ggplot2::element_text(face = "bold", size = base_size * 1.0,
                                                 color = "#0F172A"),
      axis.title.y       = ggplot2::element_text(face = "bold", size = base_size * 1.10,
                                                 color = "#0F172A",
                                                 margin = ggplot2::margin(r = 12)),
      axis.ticks         = ggplot2::element_line(color = "#CBD5E1", linewidth = 0.4),
      panel.grid.major.x = ggplot2::element_line(color = "#F8FAFC", linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "#F1F5F9", linewidth = 0.5),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(color = "#94A3B8", fill = NA, linewidth = 0.6),
      panel.spacing.x    = grid::unit(0.35, "lines"),
      panel.spacing.y    = grid::unit(0.28, "lines"),
      strip.background   = ggplot2::element_rect(fill = "#F8FAFC", color = "#94A3B8",
                                                 linewidth = 0.6),
      strip.text.x       = ggplot2::element_text(face = "bold", size = base_size * 0.68,
                                                 lineheight = 0.88, color = "#0F172A",
                                                 margin = ggplot2::margin(t = 4, b = 4, l = 1, r = 1)),
      strip.text.y       = ggplot2::element_text(angle = 0, face = "bold",
                                                 size = base_size * 0.83, color = "#0F172A"),
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.margin      = ggplot2::margin(t = 6, b = 2),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.86),
      legend.text        = ggplot2::element_text(size = base_size * 0.77),
      plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.36,
                                                 color = "#0F172A", hjust = 0.5,
                                                 margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size * 0.95, color = "#475569",
                                                 hjust = 0.5,
                                                 margin = ggplot2::margin(b = 10)),
      plot.caption       = ggplot2::element_text(size = base_size * 0.78, color = "#64748B",
                                                 hjust = 0,
                                                 margin = ggplot2::margin(t = 8)),
      plot.margin        = ggplot2::margin(t = 12, r = 14, b = 10, l = 12)
    )

  p
}

#' @rdname plot_benchmark_bubble_matrix
#' @export
plot_bubble_matrix <- plot_benchmark_bubble_matrix


# ============================================================================
# INTERNAL HELPERS
# ============================================================================

# Ingest bubble matrix input into a tidy data.frame
# @keywords internal
.ingest_bubble_data <- function(data) {
  # scSimEval_consolidated object
  if (inherits(data, "scSimEval_consolidated") ||
      (is.list(data) && "consolidated_summary_table" %in% names(data))) {
    tbl <- data$consolidated_summary_table
    return(data.frame(
      Method    = tbl$Data_Name,
      Category  = tbl$Category,
      Metric    = tbl$Metric,
      Score_Raw = if (!is.null(tbl$Score)) tbl$Score else tbl$Value,
      stringsAsFactors = FALSE
    ))
  }

  # Direct evaluation result object with benchmark_summary_table or metrics_summary_table
  if (is.list(data) && !is.data.frame(data)) {
    if ("benchmark_summary_table" %in% names(data)) {
      tbl <- data$benchmark_summary_table
      if (is.data.frame(tbl)) {
        if (!"Method" %in% colnames(tbl)) tbl$Method <- "Simulation"
        return(.ingest_bubble_data(tbl))
      }
    }
    if ("metrics_summary_table" %in% names(data)) {
      tbl <- data$metrics_summary_table
      if (is.data.frame(tbl)) {
        if (!"Method" %in% colnames(tbl)) tbl$Method <- "Simulation"
        return(.ingest_bubble_data(tbl))
      }
    }
    if ("summary_table" %in% names(data)) {
      tbl <- data$summary_table
      if (is.data.frame(tbl)) {
        if (!"Method" %in% colnames(tbl)) tbl$Method <- "Simulation"
        return(.ingest_bubble_data(tbl))
      }
    }
  }

  # Named list of individual benchmark outputs
  if (is.list(data) && !is.data.frame(data) && !is.null(names(data))) {
    rows <- list()
    for (nm in names(data)) {
      item <- data[[nm]]
      tbl  <- NULL
      if (is.list(item) && !is.null(item$benchmark_summary_table)) tbl <- item$benchmark_summary_table
      if (is.list(item) && !is.null(item$metrics_summary_table))   tbl <- item$metrics_summary_table
      if (is.list(item) && !is.null(item$summary_table))           tbl <- item$summary_table
      if (is.data.frame(item)) tbl <- item
      if (!is.null(tbl)) {
        val <- if (!is.null(tbl$Score)) tbl$Score else if (!is.null(tbl$Value)) tbl$Value else NA_real_
        mth <- if (!is.null(tbl$Method)) tbl$Method else nm
        cat_col <- if (!is.null(tbl$Category)) tbl$Category else "Uncategorized"
        rows[[nm]] <- data.frame(Method = mth, Category = cat_col,
                                  Metric = tbl$Metric, Score_Raw = val,
                                  stringsAsFactors = FALSE)
      }
    }
    if (length(rows) == 0) stop("No valid benchmark tables found in the named list.")
    return(do.call(rbind, rows))
  }

  # Pre-formed tidy data.frame
  if (is.data.frame(data)) {
    df <- data
    if (!"Method" %in% colnames(df)) {
      if ("Simulator" %in% colnames(df)) {
        df$Method <- df$Simulator
      } else if ("Data_Name" %in% colnames(df)) {
        df$Method <- df$Data_Name
      } else {
        df$Method <- "Simulation"
      }
    }
    if (!"Category" %in% colnames(df)) df$Category <- "Uncategorized"
    if (!"Score_Raw" %in% colnames(df) && "Score" %in% colnames(df)) df$Score_Raw <- df$Score
    if (!"Score_Raw" %in% colnames(df) && "Value" %in% colnames(df)) df$Score_Raw <- df$Value
    if (!"Metric" %in% colnames(df) && "Display_Metric" %in% colnames(df)) df$Metric <- df$Display_Metric
    missing <- setdiff(c("Method", "Metric", "Score_Raw"), colnames(df))
    if (length(missing) > 0)
      stop(paste("data.frame missing required columns:", paste(missing, collapse = ", ")))
    return(df)
  }

  stop("Unsupported `data` format. See ?plot_benchmark_bubble_matrix for accepted inputs.")
}


# Normalize bubble scores (direction-aware, per-metric across methods)
# @keywords internal
.normalize_bubble_scores <- function(df) {
  df$Score_Norm <- NA_real_
  for (met in unique(df$Metric)) {
    idx     <- df$Metric == met
    raw_val <- df$Score_Raw[idx]
    non_na  <- raw_val[!is.na(raw_val)]
    if (length(non_na) == 0) {
      df$Score_Norm[idx] <- NA_real_
      next
    }
    is_hib  <- met %in% .HIGHER_IS_BETTER_METRICS
    min_v   <- min(non_na)
    max_v   <- max(non_na)
    rng     <- max_v - min_v
    if (rng < 1e-10) {
      normed <- rep(1.0, sum(idx))
    } else if (is_hib) {
      normed <- (raw_val - min_v) / rng
    } else {
      normed <- 1 - (raw_val - min_v) / rng   # invert: lower raw = better = higher score
    }
    df$Score_Norm[idx] <- pmax(0, pmin(1, normed))
  }
  still_na <- is.na(df$Score_Norm)
  if (any(still_na))
    df$Score_Norm[still_na] <- pmax(0, pmin(1, df$Score_Raw[still_na]))
  df
}



# -------------------------------------------------------------------------
# Evaluation Summary Horizontal Bar Matrix
# -------------------------------------------------------------------------

#' Plot Single-Cell Simulator Evaluation Summary
#'
#' Produces a 600 DPI horizontal bar matrix ranking single-cell simulators
#' across the 8 canonical evaluation categories and overall composite performance.
#'
#' \if{html}{\figure{evaluation_summary_bars.png}{options: width="800" alt="Single-Cell Simulator Evaluation Summary"}}
#'
#' @param data A benchmark summary data frame, a list containing \code{benchmark_summary_table},
#'   or the output of \code{\link{evaluate_simulation_accuracy}}.
#' @param method_classes Optional named vector mapping method names to class/category labels
#'   (e.g., \code{c("Splatter" = "Class 1", "scDesign3" = "Class 2")}). If provided, rows are
#'   grouped into class panels.
#' @param category_colors Optional named vector specifying custom colors for category panels and overall score.
#' @param normalize_scores Logical indicating whether to normalize scores to [0, 1] (direction-aware where larger is better). Default is \code{TRUE}.
#' @param simulators_order Optional vector specifying a custom ordering of simulators. Default is ranking by descending Overall Score.
#' @param show_labels Logical indicating whether to display numerical score labels on/beside each bar. Default is \code{TRUE}.
#' @param title Character plot title. Default is \code{"Single-Cell Simulator Evaluation Summary"}.
#' @param subtitle Character plot subtitle.
#' @param base_size Numeric base font size (default 11).
#'
#' @return A \code{ggplot} object representing the multi-panel evaluation summary.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
#' @aliases plot_benchmark_summary_bars plot_summary_bars
#' @examples
#' \donttest{
#' demo <- readRDS(system.file("shiny/scSimEvalApp/data/demo_benchmark_data.rds", package = "scSimEval"))
#' p <- plot_evaluation_summary(demo$benchmark_summary_table)
#' print(p)
#' }
plot_evaluation_summary <- function(
  data,
  method_classes    = NULL,
  category_colors   = NULL,
  normalize_scores  = TRUE,
  simulators_order  = NULL,
  show_labels       = TRUE,
  title             = "Single-Cell Simulator Evaluation Summary",
  subtitle          = "Average fidelity scores across 8 canonical evaluation categories and overall benchmark performance",
  base_size         = 11
) {
  # 1. Ingest data
  df <- .ingest_bubble_data(data)

  # 2. Canonical category mapping
  if ("Category" %in% colnames(df)) {
    df$Display_Category <- ifelse(df$Category %in% names(.LEGACY_CATEGORY_MAP),
                                  .LEGACY_CATEGORY_MAP[df$Category], df$Category)
  } else {
    df$Display_Category <- "(I) Distributional Properties"
  }

  if ("Metric" %in% colnames(df)) {
    mc_override <- .METRIC_CATEGORY_MAP[df$Metric]
    idx         <- !is.na(mc_override)
    df$Display_Category[idx] <- mc_override[idx]
  }

  # 3. Score normalization
  if (normalize_scores && "Score_Raw" %in% colnames(df)) {
    df <- .normalize_bubble_scores(df)
  } else if (!"Score_Norm" %in% colnames(df)) {
    df$Score_Norm <- pmax(0, pmin(1, df$Score_Raw))
  }

  # 4. Method classes (optional)
  has_classes <- FALSE
  if (!is.null(method_classes)) {
    has_classes <- TRUE
    df$Method_Class <- as.character(method_classes[df$Method])
    df$Method_Class[is.na(df$Method_Class)] <- "Other Simulators"
  } else if ("Method_Class" %in% colnames(df) || "Class" %in% colnames(df)) {
    cls_col <- if ("Method_Class" %in% colnames(df)) "Method_Class" else "Class"
    df$Method_Class <- as.character(df[[cls_col]])
    has_classes <- TRUE
  }

  strip_map <- c(
    "(I) Distributional Properties"              = "(I)\nDistribution",
    "(II) Correlations & Zero-Inflation"         = "(II)\nCorrelation",
    "(III) Cellular Structure & Concordance"     = "(III)\nCell Structure",
    "(IV) Batch Effects & Confounder Mixing"     = "(IV)\nBatch Mixing",
    "(V) Biological Signal & Downstream Fidelity"= "(V)\nBio-Signal & DE",
    "(VI) Trajectory & Lineage Dynamics"         = "(VI)\nTraj.",
    "(VII) Cross-Modal Coupling & Modularity"    = "(VII)\nCross-Modal",
    "(VIII) Computational Scalability"           = "(VIII)\nScalability"
  )

  df$Display_Category <- as.character(df$Display_Category)
  df$Short_Category <- strip_map[df$Display_Category]
  df$Short_Category[is.na(df$Short_Category)] <- "(I)\nDistribution"

  # 5. Aggregate category means
  group_cols <- c("Method", "Short_Category")
  if (has_classes) group_cols <- c(group_cols, "Method_Class")
  cat_means <- stats::aggregate(
    stats::as.formula(paste("Score_Norm ~", paste(group_cols, collapse = " + "))),
    data = df, FUN = mean
  )
  colnames(cat_means)[colnames(cat_means) == "Short_Category"] <- "Panel"

  # 6. Aggregate overall means
  ov_group_cols <- "Method"
  if (has_classes) ov_group_cols <- c(ov_group_cols, "Method_Class")
  overall_means <- stats::aggregate(
    stats::as.formula(paste("Score_Norm ~", paste(ov_group_cols, collapse = " + "))),
    data = df, FUN = mean
  )
  overall_means$Panel <- "Overall\nScore"

  combined <- rbind(cat_means, overall_means)

  # 7. Simulator ordering
  if (!is.null(simulators_order)) {
    ranked_methods <- simulators_order
  } else if (has_classes) {
    cls_ord <- order(overall_means$Method_Class, -overall_means$Score_Norm)
    ranked_methods <- unique(overall_means$Method[cls_ord])
  } else {
    ord <- order(overall_means$Score_Norm, decreasing = TRUE)
    ranked_methods <- overall_means$Method[ord]
  }

  combined$Method <- factor(combined$Method, levels = rev(ranked_methods))

  panel_levels <- c(unname(strip_map[.CANONICAL_CATEGORY_ORDER]), "Overall\nScore")
  combined$Panel <- factor(combined$Panel, levels = panel_levels)

  if (has_classes) {
    combined$Method_Class <- factor(combined$Method_Class, levels = unique(df$Method_Class))
  }

  # 8. Colors (earlier pleasant tinted palette)
  default_colors <- c(
    "(I)\nDistribution"   = "#60A5FA",  # Sky Blue
    "(II)\nCorrelation"    = "#34D399",  # Emerald Green
    "(III)\nCell Structure"= "#A78BFA",  # Purple / Lavender
    "(IV)\nBatch Mixing"   = "#FBBF24",  # Amber / Warm Yellow
    "(V)\nBio-Signal & DE" = "#FB923C",  # Coral / Orange
    "(VI)\nTraj."          = "#2DD4BF",  # Teal / Mint
    "(VII)\nCross-Modal"   = "#818CF8",  # Indigo
    "(VIII)\nScalability"  = "#4ADE80",  # Light Green
    "Overall\nScore"       = "#94A3B8"   # Slate Grey
  )

  if (!is.null(category_colors)) {
    for (nm in names(category_colors)) {
      if (nm %in% names(default_colors)) default_colors[nm] <- category_colors[nm]
    }
  }

  combined$Label <- sprintf("%.2f", combined$Score_Norm)

  # Zebra striping
  n_methods <- length(ranked_methods)
  zebra_df <- data.frame(
    ymin = seq(0.5, n_methods - 0.5, by = 1),
    ymax = seq(1.5, n_methods + 0.5, by = 1),
    is_even = rep(c(TRUE, FALSE), length.out = n_methods)
  )
  zebra_df <- zebra_df[zebra_df$is_even, ]

  p <- ggplot2::ggplot(combined, ggplot2::aes(x = Score_Norm, y = Method))

  if (!has_classes) {
    p <- p + ggplot2::geom_rect(
      data = zebra_df,
      ggplot2::aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      fill = "#F8FAFC", inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = seq(0.5, n_methods + 0.5, by = 1), color = "#E2E8F0", linewidth = 0.5)
  }

  p <- p +
    ggplot2::geom_vline(xintercept = 0, color = "#475569", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted", color = "#CBD5E1", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "dashed", color = "#64748B", linewidth = 0.5) +
    ggplot2::geom_col(
      ggplot2::aes(fill = Panel),
      width     = 0.68,
      color     = "#1E293B",
      linewidth = 0.45,
      alpha     = 0.95
    )

  if (show_labels) {
    p <- p +
      ggplot2::geom_text(
        data = combined[combined$Score_Norm < 0.28, ],
        ggplot2::aes(x = Score_Norm + 0.02, label = Label),
        hjust = 0, size = base_size * 0.29, fontface = "bold", color = "#0F172A"
      ) +
      ggplot2::geom_text(
        data = combined[combined$Score_Norm >= 0.28, ],
        ggplot2::aes(x = Score_Norm - 0.03, label = Label),
        hjust = 1, size = base_size * 0.29, fontface = "bold", color = "#0F172A"
      )
  }

  if (has_classes) {
    p <- p + ggplot2::facet_grid(Method_Class ~ Panel, scales = "free_y", space = "free_y")
  } else {
    p <- p + ggplot2::facet_grid(. ~ Panel, scales = "fixed")
  }

  p <- p +
    ggplot2::scale_fill_manual(values = default_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.5, 1.0),
      labels = c("0", "0.5", "1.0"),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Average Evaluation Score (0 - 1, Higher = Better)",
      y        = "Ranked Simulators"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.35, hjust = 0.5, color = "#0F172A", margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5, color = "#475569", margin = ggplot2::margin(b = 14)),
      axis.title.x       = ggplot2::element_text(face = "bold", size = base_size * 0.95, color = "#1E293B", margin = ggplot2::margin(t = 8)),
      axis.title.y       = ggplot2::element_text(face = "bold", size = base_size * 1.0, color = "#0F172A", margin = ggplot2::margin(r = 10)),
      axis.text.y        = ggplot2::element_text(face = "bold", size = base_size * 0.95, color = "#0F172A"),
      axis.text.x        = ggplot2::element_text(size = base_size * 0.77, color = "#64748B"),
      strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.83, color = "#0F172A", lineheight = 1.15),
      strip.background   = ggplot2::element_rect(fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.7),
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(color = "#CBD5E1", fill = NA, linewidth = 0.7),
      panel.spacing      = grid::unit(4, "pt"),
      plot.margin        = ggplot2::margin(12, 14, 12, 12)
    )

  p
}

#' @rdname plot_evaluation_summary
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_benchmark_summary_bars <- plot_evaluation_summary

#' @rdname plot_evaluation_summary
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
plot_summary_bars <- plot_evaluation_summary


# ==============================================================================
# DIMENSION REDUCTION & EMBEDDING VISUALIZATION (UMAP / t-SNE / PCA)
# ==============================================================================

#' Compute Low-Dimensional Embeddings for Reference and Simulated Single-Cell Datasets
#'
#' Projects empirical biological reference and simulated count matrices into low-dimensional
#' representations (UMAP, t-SNE, or PCA) for side-by-side and multi-method topological benchmarking.
#' Automatically performs library-size normalization, variance-based feature selection, principal
#' component analysis, and optional unsupervised clustering.
#'
#' @param reference Biological reference count matrix (features x cells), dgCMatrix, SingleCellExperiment, or Seurat object.
#' @param simulated A single simulated count matrix or a named list of simulated count matrices (e.g. \code{list("Splatter" = sim1, "scDesign3" = sim2)}).
#' @param reduction Character string specifying the dimensionality reduction method: \code{"umap"} (default), \code{"tsne"}, or \code{"pca"}.
#' @param n_pcs Integer specifying the number of principal components to calculate (default: 30).
#' @param perplexity Numeric perplexity for t-SNE (default: 30; automatically adapted for small sample sizes).
#' @param n_neighbors Integer number of nearest neighbors for UMAP (default: 15).
#' @param min_dist Numeric minimum distance parameter for UMAP (default: 0.3).
#' @param seed Random seed for reproducibility (default: 42).
#' @param cell_types Optional factor or character vector of cell type annotations for reference cells (or named list if per-dataset).
#' @param batch Optional factor or character vector of batch annotations.
#'
#' @return A tidy \code{data.frame} containing cell coordinates (\code{Dim1}, \code{Dim2}),
#'   \code{Dataset} name, \code{Dataset_Type} ("Reference" vs "Simulated"), \code{Cell_Type},
#'   \code{Cluster}, \code{Library_Size}, \code{Detected_Features}, and \code{Batch}.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' emb <- compute_dataset_embeddings(
#'   reference = example_scrna$ref,
#'   simulated = list("Splatter" = example_scrna$sim),
#'   reduction = "umap",
#'   cell_types = example_scrna$cell_types
#' )
#' head(emb)
#' @export
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
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set_seed_fn <- get("set.seed", asNamespace("base"))
    set_seed_fn(seed)
  }
  
  extract_counts <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "SingleCellExperiment")) {
      SingleCellExperiment::counts(obj)
    } else if (inherits(obj, "Seurat")) {
      obj[["RNA"]]$counts
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
    libs <- if (inherits(mat, "dgCMatrix")) Matrix::colSums(mat) else colSums(mat)
    det_feats <- if (inherits(mat, "dgCMatrix")) Matrix::colSums(mat > 0) else colSums(mat > 0)
    
    # 2. Library size scaling and log-transformation
    scale_factor <- stats::median(libs[libs > 0])
    if (is.na(scale_factor) || scale_factor == 0) scale_factor <- 1e4
    
    mat_dense <- as.matrix(mat)
    norm_mat <- log2(sweep(mat_dense, 2, pmax(libs, 1), "/") * scale_factor + 1)
    
    # 3. Variance-based feature selection (top 2,000 highly variable features)
    vars <- apply(norm_mat, 1, stats::var)
    vars[is.na(vars)] <- 0
    top_n <- min(2000, n_feats)
    top_idx <- order(vars, decreasing = TRUE)[seq_len(top_n)]
    sub_mat <- norm_mat[top_idx, , drop = FALSE]
    
    # 4. Principal Component Analysis (PCA)
    k_pc <- min(n_pcs, n_cells - 1, top_n - 1)
    if (k_pc < 2) k_pc <- 2
    
    pca_res <- if (requireNamespace("irlba", quietly = TRUE) && k_pc < (n_cells - 2) && k_pc < (top_n - 2)) {
      tryCatch(
        irlba::prcomp_irlba(t(sub_mat), n = k_pc, center = TRUE, scale. = FALSE),
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
        paste0("Cluster_", km_fit$cluster)
      }
    } else {
      paste0("Cluster_", km_fit$cluster)
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


#' Plot High-Dimensional Dataset Embeddings for Reference and Simulated Single-Cell Data
#'
#' Generates publication-ready comparative visualizations of empirical reference and simulated single-cell
#' datasets across UMAP, t-SNE, or PCA coordinate spaces. Supports side-by-side comparison, comprehensive
#' multi-method facet grids, and co-embedded overlays, with customizable color mappings and themes.
#'
#' @param embedding_data A \code{data.frame} produced by \code{\link{compute_dataset_embeddings}}.
#' @param reduction Character string specifying reduction used: \code{"umap"}, \code{"tsne"}, or \code{"pca"}.
#' @param layout Character string specifying the comparison layout:
#'   \describe{
#'     \item{\code{"facet"}}{Faceted grid showing Reference alongside all selected simulated datasets simultaneously (default).}
#'     \item{\code{"side_by_side"}}{Direct side-by-side comparison of Reference versus a single chosen simulator.}
#'     \item{\code{"overlay"}}{Overlaid single coordinate space showing Reference and Simulated cells together.}
#'   }
#' @param color_by Feature used for coloring cells: \code{"cell_type"} (default), \code{"cluster"}, \code{"library_size"}, \code{"dataset"}, or \code{"batch"}.
#' @param selected_methods Optional character vector of simulator names to display. If \code{NULL}, displays all available simulators.
#' @param pt_size Numeric point size (default: 0.8).
#' @param alpha Numeric transparency in [0, 1] (default: 0.75).
#' @param palette Character string specifying discrete color palette from RColorBrewer (default: \code{"Set1"}).
#' @param title Optional title character string. If \code{NULL}, an informative default title is generated.
#' @param base_size Base font size for ggplot2 rendering (default: 12).
#'
#' @return A \code{ggplot2} object.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' emb <- compute_dataset_embeddings(
#'   reference = example_scrna$ref,
#'   simulated = list("Splatter" = example_scrna$sim),
#'   reduction = "umap",
#'   cell_types = example_scrna$cell_types
#' )
#' plot_dataset_embeddings(emb, reduction = "umap", layout = "facet")
#' @export
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
  
  # Strip labels with sample sizes
  tally <- table(df$Dataset)
  levels_with_n <- paste0(names(tally), " (n = ", formatC(as.integer(tally), format = "d", big.mark = ","), ")")
  names(levels_with_n) <- names(tally)
  df$Facet_Label <- factor(levels_with_n[as.character(df$Dataset)], levels = levels_with_n)
  
  color_col <- switch(color_by,
    "cell_type" = "Cell_Type",
    "cluster" = "Cluster",
    "library_size" = "Library_Size",
    "dataset" = "Dataset",
    "batch" = "Batch",
    "Cell_Type"
  )
  
  is_numeric_color <- is.numeric(df[[color_col]])
  legend_name <- gsub("_", " ", color_col)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2))
  
  if (is_numeric_color) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_col]]), size = pt_size, alpha = alpha) +
      ggplot2::scale_color_viridis_c(name = legend_name, option = "viridis")
  } else {
    n_cats <- length(unique(stats::na.omit(df[[color_col]])))
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_col]]), size = pt_size, alpha = alpha)
    
    if (n_cats <= 9 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + ggplot2::scale_color_brewer(palette = palette, name = legend_name)
    } else {
      cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_cats)), palette))(n_cats)
      p <- p + ggplot2::scale_color_manual(values = cols, name = legend_name)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  }
  
  # Dynamic title and subtitle
  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste0(red_label, " Projection: Biological Reference vs. Simulators")
  }
  
  plot_sub <- paste0("Color: ", legend_name, " | Layout: ", layout, " | Total Cells: ", formatC(nrow(df), format = "d", big.mark = ","))
  
  p <- p +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_sub,
      x = paste0(red_label, " 1"),
      y = paste0(red_label, " 2")
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


#' Compute Quantitative Quality Metrics for Low-Dimensional Cell Embeddings
#'
#' Computes quantitative cluster separability (Silhouette width), unsupervised cluster recovery
#' (Adjusted Rand Index), and distribution preservation metrics (Mean UMI, Mean detected features,
#' sparsity percentage, and relative percentage error) between reference and simulated datasets.
#'
#' @param embedding_data A \code{data.frame} produced by \code{\link{compute_dataset_embeddings}}.
#' @param metric_coords Character vector of the 2 coordinate columns used to compute distances (default: \code{c("Dim1", "Dim2")}).
#'
#' @return A \code{data.frame} summarizing quantitative embedding quality metrics for each dataset.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' emb <- compute_dataset_embeddings(
#'   reference = example_scrna$ref,
#'   simulated = list("Splatter" = example_scrna$sim),
#'   reduction = "umap",
#'   cell_types = example_scrna$cell_types
#' )
#' compute_embedding_quality_metrics(emb)
#' @export
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
    
    rows[[d]] <- data.frame(
      Dataset = d,
      Role = unique(sub$Dataset_Type),
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


# ============================================================================
# 11. compute_method_leaderboard()
# ============================================================================

#' Compute Simulator Method Performance Leaderboard
#'
#' Computes the overall performance leaderboard ranking simulation methods across
#' all benchmark evaluation metrics. Metrics are first direction-inverted and min-max
#' standardized into fidelity scores in [0, 1] (where 1.0 represents best observed performance).
#' Simulators are then rank-ordered by mean overall fidelity score.
#'
#' @param benchmark_data A benchmark data.frame (such as \code{demo$benchmark_summary_table} or output from
#'   \code{\link{evaluate_simulation_accuracy}()}) or a named list of benchmark result tables.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{Overall_Rank}: Integer rank (1 = top-performing simulator).
#'   \item \code{Method}: Simulator method name.
#'   \item \code{Average_Fidelity}: Formatted percentage string (e.g., \code{"58.5\%"}).
#'   \item \code{Fidelity_Score}: Numeric composite fidelity score in [0, 1] rounded to 4 decimals.
#' }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' lb <- compute_method_leaderboard(res)
#' @export
compute_method_leaderboard <- function(benchmark_data) {
  df <- .ingest_bubble_data(benchmark_data)
  df_norm <- .normalize_bubble_scores(df)

  leaderboard <- stats::aggregate(
    Score_Norm ~ Method,
    data = df_norm,
    FUN  = mean,
    na.rm = TRUE
  )
  colnames(leaderboard)[2] <- "Score"

  leaderboard$Overall_Rank     <- rank(-leaderboard$Score, ties.method = "min")
  leaderboard                  <- leaderboard[order(leaderboard$Overall_Rank), ]
  leaderboard$Average_Fidelity <- sprintf("%.1f%%", leaderboard$Score * 100)
  leaderboard$Fidelity_Score   <- round(leaderboard$Score, 4)

  leaderboard[, c("Overall_Rank", "Method", "Average_Fidelity", "Fidelity_Score"), drop = FALSE]
}


