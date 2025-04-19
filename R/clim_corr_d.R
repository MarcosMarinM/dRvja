#' Calculate Daily Climate Correlation (Day-Resolution Output)
#'
#' Calculates the correlation between an annual tree-growth variable and a DAILY
#' climate variable. Computes correlations for various rolling window widths
#' (in days) ending on specific days relative to the growth year.
#' Day 1 is Jan 1 of the previous year, Day 730 is approx Dec 31 of the growth year.
#' Returns results formatted for heatmap analysis (wide matrices: width vs end_day)
#' and optionally saves significant correlations. Requires daily climate input.
#' VERY COMPUTATIONALLY INTENSIVE.
#'
#' @param growth_data A data frame with tree growth data ('year', `growth_var`).
#' @param climate_data A data frame with DAILY climate data. Must contain date and value columns.
#' @param date_col Character string. Name of the date column (must be Date class)
#'   in `climate_data`. Default: "date".
#' @param value_col Character string. Name of the numeric value column in `climate_data`.
#'   Default: "value".
#' @param growth_var Character string. Growth index column name (default: "res").
#' @param climate_var_name Character string. Base name for climate variable (default: "climate_d").
#'   Used for internal messages and default output filename.
#' @param window_widths_days Integer vector. Sequence of moving window widths in DAYS
#'   to test (e.g., `30:365`). Default: `30:365`.
#' @param end_days Integer vector. Sequence of relative end days for the window.
#'   Day 1 = Jan 1 of year BEFORE growth year; Day 730 approx Dec 31 of growth year.
#'   Default: `1:730`.
#' @param start_year Optional integer. First year for analysis.
#' @param end_year Optional integer. Last year for analysis.
#' @param cor_method Character string. Correlation method (default: "pearson").
#' @param daily_roll_fun Function or Character string. Function for daily rolling window
#'   ('sum' or 'mean'). Default: "sum".
#' @param sig Numeric. Threshold for p-value significance (default: 0.05).
#' @param save_significant Logical. Save significant results to file? (Default: TRUE).
#' @param output_filename Character string or NULL. Path/filename for significant
#'   results (.txt). If NULL, default name is generated.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `correlations`: Wide data frame (rows=window_width_days, cols=end_day) of correlations.
#'     \item `p_values`: Wide data frame of p-values.
#'     \item `n_overlap`: Wide data frame of overlap counts (N).
#'     \item `significant_tidy`: Tidy data frame (p < `sig`). Columns:
#'       window_width_days, end_day, correlation, p_value, n_overlap.
#'   }
#'   Returns NULL if calculations fail. Shows progress bar.
#' @export
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @importFrom zoo rollsumr rollmeanr
#' @importFrom stats cor.test complete.cases sd
#' @importFrom lubridate year month day is.Date days as_date
#' @importFrom purrr possibly map
#' @importFrom rlang .data :=
#' @importFrom utils write.table
#' @importFrom pbapply pblapply
#'
#' @examples
#' \dontrun{ # VERY SLOW example
#' # --- Generate Sample Data ---
#' growth <- data.frame(year = 1980:2000, res = rnorm(21)) # Shorter period
#' # Daily Climate Data
#' dates_d <- seq(as.Date("1978-01-01"), as.Date("2001-12-31"), by = "day")
#' climate_d <- data.frame(my_date = dates_d, my_value = runif(length(dates_d), 0, 10))
#' # No attributes needed now, columns specified directly
#'
#' # --- Run clim_corr_d (SMALL example range) ---
#' temp_dir <- tempdir()
#' fname <- file.path(temp_dir, "sig_corr_d_prec_d.txt")
#' results_list_d <- clim_corr_d(
#'   growth_data = growth,
#'   climate_data = climate_d,
#'   date_col = "my_date",         # Specify column names
#'   value_col = "my_value",
#'   climate_var_name = "prec_d",
#'   window_widths_days = c(30, 60, 90), # Small set of widths
#'   end_days = seq(300, 400, by=10),    # Small, sparse set of end days
#'   daily_roll_fun = "sum",
#'   output_filename = fname
#' )
#'
#' # --- Inspect results ---
#' if (!is.null(results_list_d)) {
#'   print("Correlation Matrix (Daily, Width vs End Day):")
#'   print(results_list_d$correlations)
#'   print("Significant Tidy Results (Daily, p < 0.05):")
#'   print(results_list_d$significant_tidy)
#'   if(file.exists(fname)) {
#'      print(paste("Significant results saved to:", fname))
#'      unlink(fname) # Clean up temp file
#'   }
#' }
#' } # End dontrun
#'
clim_corr_d <- function(growth_data,
                            climate_data,
                            date_col = "date",         # New param
                            value_col = "value",       # Now param (was already)
                            growth_var = "res",
                            climate_var_name = "climate_d",
                            window_widths_days = 30:365, # Default changed
                            end_days = 1:730,          # New param
                            start_year = NULL,
                            end_year = NULL,
                            cor_method = "pearson",
                            daily_roll_fun = "sum",
                            sig = 0.05,
                            save_significant = TRUE,
                            output_filename = NULL
) {

  # --- Input Validation ---
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required. Please install it.", call. = FALSE)
  if (!requireNamespace("pbapply", quietly = TRUE)) stop("Package 'pbapply' is required. Please install it.", call. = FALSE)

  if (!is.data.frame(growth_data) || !all(c("year", growth_var) %in% colnames(growth_data))) stop("`growth_data` must be a data frame with 'year' and '", growth_var, "' columns.", call. = FALSE)
  if (!is.data.frame(climate_data)) stop("`climate_data` must be a data frame.", call. = FALSE)

  # *** Validate climate data columns and types ***
  if (!date_col %in% colnames(climate_data)) stop("Specified `date_col` '", date_col, "' not found in `climate_data`.", call. = FALSE)
  if (!value_col %in% colnames(climate_data)) stop("Specified `value_col` '", value_col, "' not found in `climate_data`.", call. = FALSE)
  # Use tryCatch for as_date conversion robustness
  climate_data[[date_col]] <- tryCatch(
    lubridate::as_date(climate_data[[date_col]]),
    error = function(e) stop("Column `", date_col, "` could not be converted to Date class. Check format or content. Original error: ", e$message, call.=FALSE),
    warning = function(w) {message("Warning converting `", date_col, "` to Date: ", w$message); lubridate::as_date(climate_data[[date_col]])}
  )
  if (!lubridate::is.Date(climate_data[[date_col]])) stop("Column `", date_col, "` must be or be convertible to Date class.", call. = FALSE) # Redundant check after tryCatch
  if (!is.numeric(climate_data[[value_col]])) stop("Column `", value_col, "` must be numeric.", call. = FALSE)
  # *** End climate data validation ***

  if (!is.numeric(growth_data$year)) stop("'year' column in `growth_data` must be numeric.", call.=FALSE)
  if (!is.numeric(growth_data[[growth_var]])) stop("Growth variable column '", growth_var, "' must be numeric.", call.=FALSE)
  if (!is.numeric(window_widths_days) || any(window_widths_days <= 0) || any(window_widths_days != floor(window_widths_days))) stop("`window_widths_days` must be a vector/sequence of positive integers.", call. = FALSE)
  window_widths_days <- sort(unique(as.integer(window_widths_days)))
  if(length(window_widths_days) == 0) stop("`window_widths_days` resulted in an empty sequence.", call.=FALSE)
  if (!is.numeric(end_days) || any(end_days <= 0) || any(end_days != floor(end_days))) stop("`end_days` must be a vector/sequence of positive integers.", call. = FALSE)
  end_days <- sort(unique(as.integer(end_days)))
  if(length(end_days) == 0) stop("`end_days` resulted in an empty sequence.", call.=FALSE)
  if (!is.numeric(sig) || sig <= 0 || sig >= 1) stop("`sig` (significance level) must be between 0 and 1.", call.=FALSE)
  if (!is.logical(save_significant) || length(save_significant) != 1) stop("`save_significant` must be TRUE or FALSE.", call.=FALSE)
  if (!is.null(output_filename) && (!is.character(output_filename) || length(output_filename) != 1 || output_filename == "")) stop("`output_filename` must be NULL or a single character string.", call.=FALSE)

  # --- Validate daily_roll_fun ---
  if (is.character(daily_roll_fun)) {
    if (daily_roll_fun == "sum") {
      roll_fun_to_use <- zoo::rollsumr
      roll_fun_name <- "sum"
    } else if (daily_roll_fun == "mean") {
      roll_fun_to_use <- zoo::rollmeanr
      roll_fun_name <- "mean"
    } else {
      stop("Invalid character `daily_roll_fun`. Use 'sum' or 'mean'. For custom functions, pass the function itself.", call. = FALSE)
    }
  } else if (is.function(daily_roll_fun)) {
    roll_fun_to_use <- daily_roll_fun
    roll_fun_name <- deparse(substitute(daily_roll_fun))
  } else {
    stop("Invalid `daily_roll_fun`. Must be 'sum', 'mean', or a function.", call. = FALSE)
  }
  message("Using daily rolling function: ", roll_fun_name)
  message("Testing window widths (days): ", min(window_widths_days), " to ", max(window_widths_days))
  message("Testing end days (relative): ", min(end_days), " to ", max(end_days))
  message("Significance level (alpha / sig): ", sig)

  # --- Prepare Climate Data (using data.table) ---
  message("Preparing climate data...")
  climate_dt <- data.table::as.data.table(climate_data)[, .(date = get(date_col), value = get(value_col))] # Select and rename cols
  climate_dt <- climate_dt[!is.na(date) & !is.na(value)] # Filter NAs
  data.table::setorder(climate_dt, date)                # Ensure sorted
  climate_dt[, year := lubridate::year(date)]           # Add year column

  # --- Prepare Growth Data (using data.table) ---
  growth_dt <- data.table::as.data.table(growth_data)[, .(year = as.integer(year), growth = get(growth_var))] # Select and rename
  growth_dt <- growth_dt[!is.na(growth)]

  # --- Apply Year Filters ---
  if (!is.null(start_year)) {
    climate_dt <- climate_dt[year >= start_year]
    growth_dt <- growth_dt[year >= start_year]
  }
  if (!is.null(end_year)) {
    climate_dt <- climate_dt[year <= end_year]
    growth_dt <- growth_dt[year <= end_year]
  }

  # --- Check for Overlap and Sufficient Data ---
  valid_growth_years <- growth_dt$year
  if(length(valid_growth_years) < 3) {
    warning("Fewer than 3 growth data points remain after filtering years.", call. = FALSE); return(NULL)
  }
  # Calculate date range needed based on growth years and end_days
  min_req_date <- as.Date(paste0(min(valid_growth_years) - 1, "-01-01")) + min(end_days) - 1 - max(window_widths_days) # Earliest possible start date needed for rolling window
  max_req_date <- as.Date(paste0(max(valid_growth_years) - 1, "-01-01")) + max(end_days) - 1 # Latest possible end date needed
  climate_dt <- climate_dt[date >= min_req_date & date <= max_req_date] # Filter climate data to relevant range

  if(nrow(climate_dt) == 0) {
    warning("No climate data remains within the required date range.", call. = FALSE); return(NULL)
  }
  message("Climate data prepared. Growth years: ", length(valid_growth_years), ". Climate days: ", nrow(climate_dt))


  # --- Set Key for Fast Lookups ---
  data.table::setkey(climate_dt, date)
  data.table::setkey(growth_dt, year)


  # --- Outer Loop: Iterate over Window Widths with Progress Bar ---
  message("Calculating correlations (this may take a long time)...")

  all_results_list <- pbapply::pblapply(window_widths_days, function(w) {

    # --- 1. Calculate Rolling Column for Width 'w' ---
    rolling_col_name <- paste0("roll_", w, "d")
    # Calculate rolling value, handling potential errors within the calculation if needed
    climate_dt[, (rolling_col_name) := tryCatch(
      roll_fun_to_use(value, k = w, fill = NA, align = "right"),
      error = function(e) {
        # warning("Rolling function failed for width ", w, "d. Error: ", e$message, call.=FALSE)
        rep(NA_real_, .N) # Return NA if function fails
      }
    )]

    # If all rolling values are NA (e.g., function failed consistently), skip this width
    if (climate_dt[, all(is.na(get(rolling_col_name)))]) {
      # warning("Rolling calculation resulted in all NAs for width ", w, "d. Skipping.")
      return(NULL)
    }

    # --- 2. Inner Calculation: Iterate over End Days (using lapply) ---
    results_for_width_w <- lapply(end_days, function(d) {

      # Calculate target dates for all growth years for this specific end day 'd'
      # Base date is Jan 1 of the *previous* year
      target_dates <- lubridate::as_date(paste0(valid_growth_years - 1, "-01-01")) + lubridate::days(d - 1)

      # Fast lookup of climate values (and year) on target dates using data.table join
      # Include year from climate_dt to ensure correct join if target_dates has duplicates (unlikely here)
      climate_values_for_d <- climate_dt[J(target_dates), .(date, year, climate_roll = get(rolling_col_name)), nomatch=NA] # Keep all target_dates, fill climate with NA if date not found
      climate_values_for_d[, lookup_year := lubridate::year(date)] # Year associated with the looked-up climate value

      # Join growth data (all valid years) with the looked-up climate values
      # We need the year from the *growth* perspective to link correctly
      # Create a temporary table with growth_year and target_date
      target_lookup <- data.table(year = valid_growth_years, target_date = target_dates)
      setkey(target_lookup, target_date)
      setkey(climate_values_for_d, date) # Key for joining on date

      # Join climate values to the target lookup table
      climate_matched_to_target <- climate_values_for_d[target_lookup, on = .(date = target_date)]
      # Now join this with growth data based on the original growth year
      setkey(climate_matched_to_target, year) # Key on growth year
      merged_for_d <- growth_dt[climate_matched_to_target, on = "year", nomatch = 0] # Inner join

      # Calculate correlation if enough non-NA pairs
      valid_pairs <- merged_for_d[!is.na(climate_roll) & !is.na(growth)]
      n_complete <- nrow(valid_pairs)

      if (n_complete >= 3) {
        test <- tryCatch(
          stats::cor.test(valid_pairs$climate_roll, valid_pairs$growth, method = cor_method, use = "na.or.complete"), # use should be covered by filtering
          error = function(e) NULL)

        if (!is.null(test)) {
          return(data.table::data.table( # Return data.table for faster binding later
            window_width_days = w,
            end_day = d,
            correlation = test$estimate,
            p_value = test$p.value,
            n_overlap = n_complete
          ))
        }
      }
      return(NULL) # Return NULL if correlation fails or not enough data
    }) # End inner lapply over 'd'

    # --- 3. Combine results for this width 'w' ---
    results_dt_for_w <- data.table::rbindlist(results_for_width_w) # Fast binding

    # --- 4. Clean up the rolling column to save memory ---
    climate_dt[, (rolling_col_name) := NULL]

    if (nrow(results_dt_for_w) > 0) {
      return(results_dt_for_w)
    } else {
      return(NULL)
    }

  }) # End outer pblapply over 'w'


  # --- Combine results from all widths ---
  all_results_tidy <- data.table::rbindlist(all_results_list)

  # --- Check if any results were produced ---
  if (nrow(all_results_tidy) == 0) {
    warning("Correlation calculation failed or produced no results for all specified inputs.", call. = FALSE)
    return(NULL)
  }
  message("Tidy results calculation complete.")

  # Convert to data.frame for consistency if needed, though not strictly necessary
  all_results_tidy_df <- as.data.frame(all_results_tidy)


  # --- Filter Significant Results & Optionally Save ---
  significant_results_tidy <- all_results_tidy_df %>%
    filter(.data$p_value < sig) %>%
    arrange(.data$window_width_days, .data$end_day)

  if (save_significant) {
    if (nrow(significant_results_tidy) > 0) {
      if (is.null(output_filename) || output_filename == "") {
        safe_clim_name <- gsub("[^[:alnum:]_.-]", "_", climate_var_name)
        output_filename_used <- paste0("significant_correlations_d_", safe_clim_name, ".txt")
        message("`output_filename` not specified, using default: ", output_filename_used)
      } else {
        output_filename_used <- output_filename
      }
      message("Attempting to save significant results (p < ", sig, ") to: ", output_filename_used)
      tryCatch({
        utils::write.table(significant_results_tidy, file = output_filename_used,
                           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
        message("Significant results successfully saved.")
      }, error = function(e) {
        warning("Failed to save significant results to '", output_filename_used, "'.\nError: ", e$message, call. = FALSE)
      })
    } else {
      message("No significant results found (p < ", sig, "). Nothing to save.")
    }
  } else {
    message("Saving significant results skipped (`save_significant` is FALSE).")
  }

  # --- Reshape ALL Tidy Data to Wide Format ---
  message("Reshaping full results to wide format...")

  # Pivot requires data.frame, not data.table generally
  pivot_data <- all_results_tidy_df

  # Pivot Correlations (Rows: width, Columns: end_day)
  correlations_wide <- pivot_data %>%
    select(window_width_days, end_day, correlation) %>%
    tidyr::pivot_wider(
      id_cols = window_width_days,
      names_from = end_day,
      names_prefix = "d_", # Add prefix to ensure valid column names
      values_from = correlation,
      names_sort = TRUE
    ) %>%
    arrange(.data$window_width_days)

  # Pivot P-values
  p_values_wide <- pivot_data %>%
    select(window_width_days, end_day, p_value) %>%
    tidyr::pivot_wider(
      id_cols = window_width_days,
      names_from = end_day,
      names_prefix = "d_",
      values_from = p_value,
      names_sort = TRUE
    ) %>%
    arrange(.data$window_width_days)

  # Pivot N Overlap
  n_overlap_wide <- pivot_data %>%
    select(window_width_days, end_day, n_overlap) %>%
    tidyr::pivot_wider(
      id_cols = window_width_days,
      names_from = end_day,
      names_prefix = "d_",
      values_from = n_overlap,
      names_sort = TRUE
    ) %>%
    arrange(.data$window_width_days)

  message("Reshaping complete.")

  # --- Return List ---
  results_list <- list(
    correlations = correlations_wide,
    p_values = p_values_wide,
    n_overlap = n_overlap_wide,
    significant_tidy = significant_results_tidy
  )

  # Clean up memory potentially held by data.table internal structures
  gc()

  return(results_list)
}
