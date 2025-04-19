#' Calculate Monthly Climate Correlation Heatmap Data & Save Significant Results
#'
#' Calculates the correlation between a tree-growth variable and a climate
#' variable over a range of moving window widths. Returns results formatted
#' for heatmap analysis (wide matrices) and optionally saves a tidy data frame
#' of statistically significant correlations (p < sig) to a file.
#'
#' @param growth_data A data frame with tree growth data ('year', `growth_var`).
#' @param climate_data A data frame with climate data (output of `read_climate_data`).
#' @param growth_var Character string. Growth index column name (default: "res").
#' @param climate_var_name Character string. Base name for climate variable (default: "climate").
#'   Used for internal messages and default output filename.
#' @param window_widths Integer vector. Sequence of moving window widths in months
#'   to test (e.g., `3:24`). Default: `3:24`.
#' @param start_year Optional integer. First year for analysis.
#' @param end_year Optional integer. Last year for analysis.
#' @param cor_method Character string. Correlation method (default: "pearson").
#' @param daily_agg_fun Function or Character string. Function for daily aggregation
#'   ('auto', 'sum', 'mean', or custom function).
#' @param value_col Character string. Name of value column if `climate_data` is daily. Default: "value".
#' @param sig Numeric. Threshold for p-value significance (default: 0.05).
#' @param save_significant Logical. Should the tidy data frame of significant
#'   results be saved to a text file? (Default: TRUE).
#' @param output_filename Character string or NULL. Path and filename for the saved
#'   significant results (.txt). If NULL (default), a name like
#'   "significant_correlations_[climate_var_name].txt" is generated in the
#'   working directory.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `correlations`: Wide data frame (rows=window_width, cols=month_abbr) of correlations.
#'     \item `p_values`: Wide data frame of p-values.
#'     \item `n_overlap`: Wide data frame of overlap counts (N).
#'     \item `significant_tidy`: Tidy data frame containing only rows where
#'       p_value < `sig`. Columns include window_width, end_month,
#'       month_abbr, correlation, p_value, n_overlap.
#'   }
#'   Returns NULL if calculations fail for all windows. As a side effect, may save a .txt file.
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom zoo rollsumr
#' @importFrom stats cor.test complete.cases aggregate sd
#' @importFrom lubridate year month day is.Date
#' @importFrom purrr map_dfr possibly
#' @importFrom rlang .data
#' @importFrom utils write.table
#'
#' @examples
#' # --- Generate Sample Data ---
#' growth <- data.frame(year = 1950:2000, res = rnorm(51))
#' climate_m <- data.frame(year = 1949:2001); for(m in 1:12) climate_m[[as.character(m)]]<-runif(53,10,60)
#' attr(climate_m, "frequency") <- "monthly"; attr(climate_m, "variable_type") <- "P"
#'
#' # --- Run clim_corr_m (default saving) ---
#' temp_dir <- tempdir() # Use temporary directory for example file
#' default_fname <- file.path(temp_dir, "significant_correlations_prec_m.txt")
#' results_list <- clim_corr_m(growth, climate_m,
#'                             window_widths = 6:18,
#'                             climate_var_name = "prec_m",
#'                             sig = 0.10, # Use different significance for example
#'                             output_filename = default_fname) # Save to temp dir
#'
#' # --- Inspect results ---
#' if (!is.null(results_list)) {
#'   print("Significant Tidy Results (p < 0.10):")
#'   print(results_list$significant_tidy)
#'   if(file.exists(default_fname)) {
#'      print(paste("Significant results saved to:", default_fname))
#'      # print(readLines(default_fname, n=5)) # View file content
#'      unlink(default_fname) # Clean up temp file
#'   } else {
#'      print("No significant results found to save (p < 0.10).")
#'   }
#'
#'   # Example: Run without saving
#'   # results_nosave <- clim_corr_m(growth, climate_m, save_significant = FALSE)
#'   # print(results_nosave$significant_tidy) # Still available in the list
#' }
#'
clim_corr_m <- function(growth_data,
                        climate_data,
                        growth_var = "res",
                        climate_var_name = "climate",
                        window_widths = 1:24,
                        start_year = NULL,
                        end_year = NULL,
                        cor_method = "pearson",
                        daily_agg_fun = "auto",
                        value_col = "value",
                        sig = 0.05,
                        save_significant = TRUE,
                        output_filename = NULL
) {

  # --- Input Validation ---
  if (!is.data.frame(growth_data) || !all(c("year", growth_var) %in% colnames(growth_data))) stop("`growth_data` must be a data frame with 'year' and '", growth_var, "' columns.", call. = FALSE)
  if (!is.data.frame(climate_data)) stop("`climate_data` must be a data frame.", call. = FALSE)
  if (!is.numeric(growth_data$year)) stop("'year' column in `growth_data` must be numeric.", call.=FALSE)
  if (!is.numeric(growth_data[[growth_var]])) stop("Growth variable column '", growth_var, "' must be numeric.", call.=FALSE)
  if (!is.numeric(window_widths) || any(window_widths <= 0) || any(window_widths != floor(window_widths))) stop("`window_widths` must be a vector or sequence of positive integers.", call. = FALSE)
  window_widths <- sort(unique(as.integer(window_widths)))
  if(length(window_widths) == 0) stop("`window_widths` resulted in an empty sequence.", call.=FALSE)
  if (!is.numeric(sig) || sig <= 0 || sig >= 1) stop("`sig` (significance level) must be between 0 and 1.", call.=FALSE) # Check renamed parameter
  if (!is.logical(save_significant) || length(save_significant) != 1) stop("`save_significant` must be TRUE or FALSE.", call.=FALSE)
  if (!is.null(output_filename) && (!is.character(output_filename) || length(output_filename) != 1 || output_filename == "")) stop("`output_filename` must be NULL or a single character string.", call.=FALSE)

  message("Calculating correlations for window widths: ", min(window_widths), " to ", max(window_widths), " months.")
  message("Significance level (alpha / sig): ", sig) # Message updated

  # --- Detect Climate Data Frequency & Variable Type ---
  frequency <- attr(climate_data, "frequency")
  variable_type <- attr(climate_data, "variable_type")
  is_daily <- FALSE
  climate_monthly_long <- NULL

  if (!is.null(frequency) && frequency == "daily") is_daily <- TRUE
  else if (!is.null(frequency) && frequency == "monthly") is_daily <- FALSE
  else {
    message("Frequency attribute missing, guessing based on columns...")
    if (all(c("date", value_col) %in% colnames(climate_data)) && lubridate::is.Date(climate_data$date)) { is_daily <- TRUE; message("Found 'date','value', assuming daily.") }
    else if (all(c("year", "month", "day", value_col) %in% colnames(climate_data))) { is_daily <- TRUE; message("Found 'year','month','day','value', assuming daily.") }
    else if (all(c("year", as.character(1:12)) %in% colnames(climate_data))) { is_daily <- FALSE; message("Found 'year','1'..'12', assuming monthly.") }
    else stop("Could not determine climate data frequency (daily/monthly).", call. = FALSE)
  }

  # --- Determine Aggregation Function for Daily Data ---
  agg_fun_to_use <- NULL
  if (is_daily) {
    if (is.character(daily_agg_fun) && daily_agg_fun == "auto") {
      if (is.null(variable_type) || !(variable_type %in% c("P", "T"))) {
        warning("Cannot auto-determine aggregation function: variable_type attribute is '",
                variable_type, "' or missing. Defaulting to 'sum'. Specify 'daily_agg_fun'.", call. = FALSE)
        agg_fun_to_use <- sum; agg_fun_name <- "sum"
      } else if (variable_type == "P") {
        agg_fun_to_use <- sum; agg_fun_name <- "sum"; message("Variable type is 'P', using 'sum' for daily aggregation.")
      } else if (variable_type == "T") {
        agg_fun_to_use <- mean; agg_fun_name <- "mean"; message("Variable type is 'T', using 'mean' for daily aggregation.")
      }
    } else if (is.function(daily_agg_fun)) {
      agg_fun_to_use <- daily_agg_fun; agg_fun_name <- deparse(substitute(daily_agg_fun)); message("Using custom function '", agg_fun_name, "' for daily aggregation.")
    } else if (is.character(daily_agg_fun) && daily_agg_fun == "sum") {
      agg_fun_to_use <- sum; agg_fun_name <- "sum"; message("Using 'sum' (specified by user) for daily aggregation.")
    } else if (is.character(daily_agg_fun) && daily_agg_fun == "mean") {
      agg_fun_to_use <- mean; agg_fun_name <- "mean"; message("Using 'mean' (specified by user) for daily aggregation.")
    } else {
      stop("Invalid 'daily_agg_fun'. Should be 'auto', 'sum', 'mean', or a function.", call. = FALSE)
    }
  }

  # --- Aggregate Daily to Monthly OR Pivot Monthly Wide to Long ---
  if (is_daily) {
    message("Aggregating daily data to monthly using '", agg_fun_name, "'.")
    value_col_present <- value_col %in% colnames(climate_data)
    if(!value_col_present) stop("Value column '", value_col,"' not found in daily climate_data.", call.=FALSE)
    if(!is.numeric(climate_data[[value_col]])) stop("Value column '", value_col,"' must be numeric.", call.=FALSE)
    if (all(c("date", value_col) %in% colnames(climate_data))) {
      climate_monthly_long <- climate_data %>%
        filter(!is.na(.data$date) & !is.na(.data[[value_col]])) %>%
        mutate(year = lubridate::year(.data$date), month = lubridate::month(.data$date)) %>%
        group_by(.data$year, .data$month) %>%
        summarise(monthly_value = agg_fun_to_use(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
    } else if (all(c("year", "month", "day", value_col) %in% colnames(climate_data))) {
      climate_monthly_long <- climate_data %>%
        filter(!is.na(.data$year) & !is.na(.data$month) & !is.na(.data$day) & !is.na(.data[[value_col]])) %>%
        group_by(.data$year, .data$month) %>%
        summarise(monthly_value = agg_fun_to_use(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
    } else stop("Daily data detected, but required columns not found.", call.=FALSE)
    message("Aggregation complete.")
  } else { # Is monthly
    message("Pivoting monthly data to long format...")
    if(!all(as.character(1:12) %in% colnames(climate_data))) stop("Monthly data must have columns '1' through '12'.", call.=FALSE)
    month_cols <- climate_data %>% select(all_of(as.character(1:12)))
    if(!all(sapply(month_cols, is.numeric))) stop("Monthly columns '1' through '12' must be numeric.", call.=FALSE)
    if(!("year" %in% colnames(climate_data)) || !is.numeric(climate_data$year)) stop("Monthly data must have a numeric 'year' column.", call.=FALSE)
    climate_monthly_long <- climate_data %>%
      select(year, all_of(as.character(1:12))) %>%
      pivot_longer(cols = all_of(as.character(1:12)), names_to = "month", values_to = "monthly_value") %>%
      mutate(month = as.integer(.data$month)) %>% filter(!is.na(monthly_value))
  }

  # --- Filter Years (on standardized monthly long data) ---
  climate_proc <- climate_monthly_long
  growth_filt <- growth_data
  climate_proc$year <- as.integer(climate_proc$year)
  growth_filt$year <- as.integer(growth_filt$year)
  if (!is.null(start_year)) { climate_proc <- climate_proc %>% filter(.data$year >= start_year); growth_filt <- growth_filt %>% filter(.data$year >= start_year) }
  if (!is.null(end_year)) { climate_proc <- climate_proc %>% filter(.data$year <= end_year); growth_filt <- growth_filt %>% filter(.data$year <= end_year) }

  # --- Define Function to Calculate Correlation for ONE Window Width ---
  calculate_for_single_width <- function(w_width, climate_monthly_base, growth_base, growth_var_name, climate_var_base_name, corr_method_used) {

    # message("Processing window width: ", w_width, " months...") # Reduced verbosity
    aggregated_climate_col_name <- paste0(climate_var_base_name, "_", w_width, "m")

    # Calculate Rolling Sum for this specific width 'w_width'
    climate_rolling <- climate_monthly_base %>%
      arrange(.data$year, .data$month) %>%
      mutate(
        "{aggregated_climate_col_name}" := zoo::rollsumr(.data$monthly_value, k = w_width, fill = NA)
      ) %>%
      filter(!is.na(.data[[aggregated_climate_col_name]])) # Remove NAs generated by rollsum

    if (nrow(climate_rolling) == 0) {
      # warning("Window width ", w_width, "m: No climate data remaining after calculating rolling sum.", call.=FALSE) # Reduced verbosity
      return(NULL) # Return NULL for this width if no data
    }

    # Find Common Years & Merge for this specific width
    common_years <- intersect(growth_base$year, climate_rolling$year)
    if (length(common_years) < 3) {
      # warning("Window width ", w_width, "m: Fewer than 3 overlapping years (", length(common_years), ").", call. = FALSE) # Reduced verbosity
      return(NULL) # Return NULL if not enough overlap
    }
    # message("Window width ", w_width, "m: Found ", length(common_years), " overlapping years (", min(common_years), "-", max(common_years), ").") # Reduced verbosity

    growth_final <- growth_base %>% filter(.data$year %in% common_years) %>% select(year, all_of(growth_var_name))
    climate_final <- climate_rolling %>% filter(.data$year %in% common_years)
    # Use inner_join to ensure only common years are kept and data is aligned
    merged_data <- inner_join(climate_final, growth_final, by = "year")

    # Calculate Correlations by Month for this width
    correlation_results <- merged_data %>%
      group_by(end_month = .data$month) %>%
      summarise(
        test_result = list(tryCatch(
          stats::cor.test(x = .data[[aggregated_climate_col_name]], y = .data[[growth_var_name]], method = corr_method_used, use = "complete.obs"),
          error = function(e) {
            # warning("Window width ", w_width, "m, Month ", unique(.data$end_month), ": Correlation failed: ", e$message, call. = FALSE) # Reduced verbosity
            return(NULL) # Return NULL for this specific month if cor.test fails
          }
        )),
        .groups = 'drop'
      ) %>%
      mutate(
        correlation = sapply(.data$test_result, function(res) if(!is.null(res)) res$estimate else NA_real_),
        p_value = sapply(.data$test_result, function(res) if(!is.null(res)) res$p.value else NA_real_),
        # Recalculate n_overlap safely after potential NA removal by complete.obs
        n_overlap = sapply(.data$end_month, function(m) {
          sub_data <- merged_data %>% filter(.data$month == m)
          sum(stats::complete.cases(sub_data[[aggregated_climate_col_name]], sub_data[[growth_var_name]]))
        })
      ) %>%
      filter(!is.na(correlation)) %>% # Keep only if correlation could be calculated
      mutate(
        # Ensure month_abbr is a factor with correct levels for pivoting later
        month_abbr = factor(month.abb[.data$end_month], levels = month.abb),
        climate_var = aggregated_climate_col_name, # Keep for potential reference if needed later? Not strictly necessary for output.
        window_width = w_width
      ) %>%
      # Select columns needed for the final wide format + identifiers
      select(window_width, end_month, month_abbr, correlation, p_value, n_overlap) %>%
      arrange(end_month)

    if(nrow(correlation_results) == 0) {
      # warning("Window width ", w_width, "m: Correlation calculation successful for 0 months.", call. = FALSE) # Reduced verbosity
      return(NULL) # Return NULL if no months were successful for this width
    }

    # message("Window width ", w_width, "m: Correlation calculation complete.") # Reduced verbosity
    return(correlation_results)
  } # End of helper function definition

  # --- Loop Over Window Widths and Combine Tidy Results ---
  message("Calculating tidy results for all specified window widths...")
  safe_calculate <- purrr::possibly(calculate_for_single_width, otherwise = NULL, quiet = FALSE)

  all_results_tidy <- purrr::map_dfr(
    window_widths,
    ~ safe_calculate(
      w_width = .x,
      climate_monthly_base = climate_proc,
      growth_base = growth_filt,
      growth_var_name = growth_var,
      climate_var_base_name = climate_var_name,
      corr_method_used = cor_method
    )
  )

  # --- Check if any results were produced AT ALL ---
  if (is.null(all_results_tidy) || nrow(all_results_tidy) == 0) {
    warning("Correlation calculation failed or produced no results for all specified window widths.", call. = FALSE)
    return(NULL)
  }
  message("Tidy results calculation complete.")

  # --- Filter Significant Results & Optionally Save ---
  significant_results_tidy <- all_results_tidy %>%
    filter(.data$p_value < sig) %>% # Use renamed parameter `sig` for filtering
    arrange(.data$window_width, .data$end_month) # Order nicely

  if (save_significant) {
    if (nrow(significant_results_tidy) > 0) {
      # Determine filename
      if (is.null(output_filename) || output_filename == "") {
        # Sanitize climate_var_name for use in filename (basic)
        safe_clim_name <- gsub("[^[:alnum:]_.-]", "_", climate_var_name)
        output_filename_used <- paste0("significant_correlations_", safe_clim_name, ".txt")
        message("`output_filename` not specified, using default: ", output_filename_used)
      } else {
        output_filename_used <- output_filename
      }

      # Try to save the file
      message("Attempting to save significant results (p < ", sig, ") to: ", output_filename_used) # Use sig in message
      tryCatch({
        utils::write.table(significant_results_tidy,
                           file = output_filename_used,
                           sep = "\t",         # Tab separated
                           row.names = FALSE, # No row names
                           col.names = TRUE,  # Include column headers
                           quote = FALSE,     # No quotes around strings
                           na = "NA")         # How to represent NA
        message("Significant results successfully saved.")
      }, error = function(e) {
        warning("Failed to save significant results to '", output_filename_used, "'.\nError: ", e$message, call. = FALSE)
      })
    } else {
      message("No significant results found (p < ", sig, "). Nothing to save.") # Use sig in message
    }
  } else {
    message("Saving significant results skipped (`save_significant` is FALSE).")
  }

  # --- Reshape ALL Tidy Data to Wide Format ---
  message("Reshaping full results to wide format...")
  # Ensure month_abbr is a factor with correct levels if it's not already
  if(!is.factor(all_results_tidy$month_abbr)) {
    all_results_tidy$month_abbr <- factor(all_results_tidy$month_abbr, levels = month.abb)
  } else {
    # Ensure levels are correct even if it's already a factor
    all_results_tidy$month_abbr <- factor(all_results_tidy$month_abbr, levels = month.abb)
  }

  # Pivot Correlations
  correlations_wide <- all_results_tidy %>%
    select(window_width, month_abbr, correlation) %>%
    tidyr::pivot_wider(
      id_cols = window_width,         # Rows will be window widths
      names_from = month_abbr,        # Columns will be month abbreviations
      values_from = correlation,      # Cells contain correlation values
      names_sort = TRUE               # Ensure columns are ordered correctly (Jan, Feb,...)
    ) %>%
    arrange(.data$window_width) # Sort rows by window width

  # Pivot P-values
  p_values_wide <- all_results_tidy %>%
    select(window_width, month_abbr, p_value) %>%
    tidyr::pivot_wider(
      id_cols = window_width,
      names_from = month_abbr,
      values_from = p_value,
      names_sort = TRUE
    ) %>%
    arrange(.data$window_width)

  # Pivot N Overlap
  n_overlap_wide <- all_results_tidy %>%
    select(window_width, month_abbr, n_overlap) %>%
    tidyr::pivot_wider(
      id_cols = window_width,
      names_from = month_abbr,
      values_from = n_overlap,
      names_sort = TRUE
    ) %>%
    arrange(.data$window_width)
  message("Reshaping complete.")

  # --- Return List (including significant tidy data) ---
  results_list <- list(
    correlations = correlations_wide,
    p_values = p_values_wide,
    n_overlap = n_overlap_wide,
    significant_tidy = significant_results_tidy # Add the filtered tidy data frame
  )

  return(results_list)
}
