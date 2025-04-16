#' Plot Comparative Climate Correlation Heatmap
#'
#' Generates a heatmap comparing climate correlation results (e.g., growth vs.
#' aggregated climate variable) across multiple climate datasets. Assumes input
#' data is the combined output from running `clim_corr` for several datasets.
#'
#' @param combined_correlation_data A data frame containing combined correlation
#'   results. Must include columns: 'month_abbr' (or 'end_month'), 'correlation',
#'   'p_value', and 'dataset' (character identifier for each dataset).
#' @param dataset_display_names Optional named character vector or list to map
#'   dataset identifiers (from the 'dataset' column) to nicer names for the
#'   plot axis. E.g., `c(era5 = "ERA5", cru = "CRU TS 4.0")`. If NULL, original
#'   dataset identifiers are used.
#' @param significance_level Numeric. The threshold for statistical significance
#'   (default: 0.05). P-values >= this level will be marked.
#' @param order_by_month Character. The month abbreviation (e.g., "Jun", "Jul")
#'   used to order the datasets vertically on the heatmap based on descending
#'   correlation in that month. If NULL, datasets are ordered alphabetically.
#'   Default: "Jun".
#' @param mark_non_significant Character. The symbol to append to correlation values
#'   where p >= significance_level (default: "*"). Set to "" to disable marking.
#' @param plot_title Character string. The main title for the plot.
#' @param plot_subtitle Character string or NULL. Subtitle for the plot. Can be
#'   automatically generated if NULL to show period and significance.
#' @param axis_label_x Character string. Label for the X-axis (months).
#' @param axis_label_y Character string. Label for the Y-axis (datasets).
#' @param legend_title Character string. Title for the color scale legend.
#' @param color_palette_diverging Character vector. Color palette from
#'   `RColorBrewer::brewer.pal` for diverging scales (used if negative
#'   correlations are present). Default: `rev(RColorBrewer::brewer.pal(11, "RdBu"))`.
#' @param color_palette_sequential Character vector. Color palette for sequential
#'   scales (used if correlations are mostly non-negative). Starts with "white".
#'   Default: `c("white", RColorBrewer::brewer.pal(9, "YlOrRd"))`.
#' @param base_font_size Numeric. Base font size for the plot theme (default: 11).
#'
#' @return A `ggplot` object representing the heatmap. The user can then print,
#'   save (using `ggplot2::ggsave`), or modify it.
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats complete.cases
#' @importFrom utils head tail
#'
#' @examples
#' # --- Create dummy combined correlation data for example ---
#' datasets <- c("Dataset_A", "Dataset_B", "Dataset_C")
#' months <- month.abb
#' dummy_data_list <- list()
#' for (d in datasets) {
#'   set.seed(which(datasets == d)) # for reproducibility
#'   dummy_data_list[[d]] <- data.frame(
#'      month_abbr = factor(months, levels = month.abb),
#'      correlation = runif(12, -0.3, 0.6),
#'      p_value = runif(12, 0, 0.2),
#'      dataset = d,
#'      n_overlap = sample(30:50, 12, replace = TRUE) # Example n
#'   )
#' }
#' combined_dummy_data <- dplyr::bind_rows(dummy_data_list)
#'
#' # --- Basic Plot ---
#' p1 <- plot_corr_db(combined_dummy_data)
#' # print(p1) # Suppressing print in examples
#'
#' # --- Plot with Options ---
#' display_names <- c(Dataset_A = "Alpha Data", Dataset_B = "Beta Climate",
#'                    Dataset_C = "Gamma Source")
#' p2 <- plot_corr_db(
#'   combined_dummy_data,
#'   dataset_display_names = display_names,
#'   order_by_month = "Jul",
#'   significance_level = 0.10,
#'   plot_title = "Example Comparative Heatmap",
#'   axis_label_x = "Final Month of Window"
#' )
#' # print(p2)
#'
plot_corr_db <- function(combined_correlation_data,
                         dataset_display_names = NULL,
                         significance_level = 0.05,
                         order_by_month = "Jun",
                         mark_non_significant = "*",
                         plot_title = "Climate Correlation Comparison",
                         plot_subtitle = NULL,
                         axis_label_x = "End Month of Window",
                         axis_label_y = "Dataset",
                         legend_title = "  r",
                         color_palette_diverging = NULL, # Set defaults later
                         color_palette_sequential = NULL, # Set defaults later
                         base_font_size = 11
) {

  # --- Input Validation ---
  req_cols <- c("correlation", "p_value", "dataset")
  # Check if either month_abbr or end_month exists
  if (!("month_abbr" %in% colnames(combined_correlation_data) || "end_month" %in% colnames(combined_correlation_data))) {
    stop("Input data must contain either 'month_abbr' or 'end_month' column.", call. = FALSE)
  }
  if (!all(req_cols %in% colnames(combined_correlation_data))) {
    stop("Input data must contain columns: ", paste(req_cols, collapse=", "), call. = FALSE)
  }
  if (!is.numeric(combined_correlation_data$correlation) || !is.numeric(combined_correlation_data$p_value)){
    stop("'correlation' and 'p_value' columns must be numeric.", call. = FALSE)
  }

  # --- Set Color Palette Defaults ---
  if (is.null(color_palette_diverging)) {
    color_palette_diverging <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  }
  if (is.null(color_palette_sequential)) {
    color_palette_sequential <- c("white", RColorBrewer::brewer.pal(9, "YlOrRd"))
  }


  # --- Data Preparation ---
  plot_data <- combined_correlation_data %>%
    # Ensure month_abbr exists and is a factor
    { if ("month_abbr" %in% colnames(.)) . else mutate(., month_abbr = month.abb[as.integer(.data$end_month)]) } %>%
    mutate(month_abbr = factor(.data$month_abbr, levels = month.abb)) %>%
    # Ensure dataset is character for reliable mapping later
    mutate(dataset = as.character(.data$dataset))

  # --- Dataset Ordering ---
  available_datasets <- unique(plot_data$dataset)
  if (!is.null(order_by_month) && order_by_month %in% month.abb) {
    message("Ordering datasets by correlation in: ", order_by_month)
    order_df <- plot_data %>%
      filter(.data$month_abbr == order_by_month & !is.na(.data$correlation)) %>%
      arrange(desc(.data$correlation))
    ordered_datasets <- unique(order_df$dataset) # Get ordered datasets with data
    # Add any remaining datasets (e.g., those with NA in order_by_month) alphabetically
    ordered_datasets_complete <- c(ordered_datasets, sort(setdiff(available_datasets, ordered_datasets)))
  } else {
    if (!is.null(order_by_month)) warning("`order_by_month` was '", order_by_month, "' but not found in month.abb. Using alphabetical order.", call. = FALSE)
    message("Ordering datasets alphabetically.")
    ordered_datasets_complete <- sort(available_datasets)
  }

  # --- Apply Display Names and Factorize ---
  if (!is.null(dataset_display_names) && (is.list(dataset_display_names) || is.vector(dataset_display_names))) {
    # Ensure it's a named vector/list
    if(is.null(names(dataset_display_names))) {
      warning("'dataset_display_names' should be a named vector/list (e.g., c(id='Display')). Using original ids.", call. = FALSE)
      display_map <- setNames(available_datasets, available_datasets)
    } else {
      display_map <- dataset_display_names
      # Add any missing datasets using their original id as the display name
      missing_in_map <- setdiff(available_datasets, names(display_map))
      if(length(missing_in_map) > 0) {
        names(missing_in_map) <- missing_in_map
        display_map <- c(display_map, missing_in_map)
      }
    }
  } else {
    # If no map provided, use original ids
    display_map <- setNames(available_datasets, available_datasets)
  }
  # Ensure display_map covers all datasets to avoid errors
  plot_data <- plot_data %>%
    mutate(display_name = unlist(display_map[.data$dataset])) # Apply mapping

  # Factorize based on calculated order and display names
  ordered_display_names <- unlist(display_map[ordered_datasets_complete])
  plot_data <- plot_data %>%
    mutate(display_name_f = factor(.data$display_name, levels = rev(ordered_display_names))) # Use rev() for ggplot y-axis

  # --- Create Formatted Labels ---
  plot_data <- plot_data %>%
    mutate(
      label_formatted = case_when(
        is.na(.data$correlation) ~ "", # No label if correlation is NA
        # Mark if p-value is NA or >= significance level
        is.na(.data$p_value) | .data$p_value >= significance_level ~ paste0(sprintf("%.3f", .data$correlation), mark_non_significant),
        TRUE ~ sprintf("%.3f", .data$correlation) # Otherwise, just the number
      )
    )

  # --- Determine Color Scale ---
  corr_range <- range(plot_data$correlation, na.rm = TRUE)
  if (all(is.na(corr_range))) {
    warning("All correlations are NA. Heatmap may be empty or grey.", call. = FALSE)
    color_limit_min <- 0
    color_limit_max <- 0
    heatmap_colors <- "grey90"
    color_limits <- c(0,0)
  } else {
    color_limit_min <- min(c(0, corr_range), na.rm = TRUE)
    color_limit_max <- max(c(0, corr_range), na.rm = TRUE)
    # Use diverging palette if significant negative values exist
    if (color_limit_min < -0.01) {
      message("Using diverging color palette (RdBu).")
      heatmap_colors <- color_palette_diverging
      abs_max <- max(abs(c(color_limit_min, color_limit_max)), na.rm = TRUE)
      color_limits <- c(-abs_max, abs_max)
    } else {
      message("Using sequential color palette (YlOrRd based).")
      heatmap_colors <- color_palette_sequential
      # Ensure lower limit is at least 0 for sequential scale
      color_limits <- c(max(0, color_limit_min), color_limit_max)
    }
  }

  # --- Generate Plot Subtitle if NULL ---
  if (is.null(plot_subtitle)) {
    plot_subtitle <- paste0("Mark '", mark_non_significant, "' indicates p >= ", significance_level)
  }


  # --- Build ggplot ---
  message("Generating ggplot object...")
  # Define text size for cells (relative to base size)
  text_in_cell_pt <- base_font_size * 0.6 # Adjust multiplier as needed
  text_in_cell_geom_size <- text_in_cell_pt / ggplot2::.pt # Convert to ggplot units

  heatmap_plot <- ggplot(plot_data, aes(x = .data$month_abbr, y = .data$display_name_f, fill = .data$correlation)) +
    geom_tile(color = "grey80", linewidth = 0.2) +
    geom_text(
      aes(label = .data$label_formatted),
      size = text_in_cell_geom_size,
      color = "black" # Consider adaptive color (white/black) based on fill later if needed
    ) +
    scale_fill_gradientn(
      colors = heatmap_colors,
      limits = color_limits,
      oob = scales::squish, # How to handle values outside limits
      na.value = "grey90", # Color for NA correlation values
      name = legend_title
    ) +
    scale_y_discrete(drop = FALSE, name = axis_label_y) + # Use factor for correct order
    scale_x_discrete(drop = FALSE, position = "bottom", name = axis_label_x) +
    coord_fixed(ratio = 1) + # Ensure square cells
    labs(
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    theme_minimal(base_size = base_font_size) +
    theme(
      axis.text = element_text(size = rel(0.8)), # Relative sizing
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.7)),
      legend.title = element_text(size = rel(0.8)),
      plot.title = element_text(hjust = 0.5, size = rel(1.1), face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.85)),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10))
    )

  message("Plot generation complete.")
  return(heatmap_plot)
}
