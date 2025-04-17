#' Calculate Monthly Climate Correlation with Moving Window (Handles Frequency)
#'
#' Calculates the correlation between a tree-growth variable and a climate
#' variable aggregated over a moving window (monthly analysis). Handles both daily
#' and monthly input climate data (output from `read_climate_data`). Daily data is
#' automatically aggregated to monthly sums before calculating the rolling window.
#'
#' @param growth_data A data frame with tree growth data ('year', `growth_var`).
#' @param climate_data A data frame with climate data, typically the output of
#'   `read_climate_data`. Can be daily (must contain 'date'/'value' or
#'   'year'/'month'/'day'/'value') or monthly ('year', '1'..'12'). The function
#'   checks the "frequency" attribute or column names to determine the type.
#' @param growth_var Character string. Growth index column name (default: "res").
#' @param climate_var_name Character string. Base name for climate variable (default: "climate").
#' @param window_width Integer. Moving window width in months (default: 12).
#' @param start_year Optional integer. First year for analysis.
#' @param end_year Optional integer. Last year for analysis.
#' @param cor_method Character string. Correlation method (default: "pearson").
#' @param daily_agg_fun Function or Character string. The function used to aggregate daily
#'   data to monthly. Default is "sum" (appropriate for precipitation). Can be "mean"
#'   (for temperature) or a custom function. Function checks `variable_type` attribute
#'   from `climate_data` to guess if not specified ("sum" for "P", "mean" for "T").
#' @param value_col Character string. Name of the value column in `climate_data` if
#'   it's daily (used internally if daily). Default is "value".
#'
#' @return A tidy data frame summarizing the monthly correlation results (columns:
#'   end_month, month_abbr, correlation, p_value, n_overlap, climate_var),
#'   or NULL if calculations fail.
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom zoo rollsumr
#' @importFrom stats cor.test complete.cases aggregate sd
#' @importFrom lubridate year month day is.Date
#' @importFrom utils head tail
#'
#' @examples
#' # --- Generate Sample Data ---
#' growth <- data.frame(year = 1950:2000, res = rnorm(51))
#' # Monthly Climate Data
#' climate_m <- data.frame(year = 1949:2001); for(m in 1:12) climate_m[[as.character(m)]]<-runif(53,10,60)
#' attr(climate_m, "frequency") <- "monthly"; attr(climate_m, "variable_type") <- "P"
#' # Daily Climate Data (Temp)
#' dates_d <- seq(as.Date("1949-01-01"), as.Date("2001-12-31"), by = "day")
#' climate_d_t <- data.frame(date = dates_d, value = rnorm(length(dates_d), 15, 5))
#' attr(climate_d_t, "frequency") <- "daily"; attr(climate_d_t, "variable_type") <- "T"
#'
#' # --- Run clim_corr_m with MONTHLY input ---
#' results_m <- clim_corr_m(growth, climate_m, climate_var_name = "prec_m")
#' print("Results from Monthly Input:")
#' print(results_m)
#'
#' # --- Run clim_corr_m with DAILY TEMP input (should use mean) ---
#' results_d_t <- clim_corr_m(growth, climate_d_t, climate_var_name = "temp_d")
#' print("Results from Daily Temp Input (aggregated with mean):")
#' print(results_d_t)
#'
#' # --- Run clim_corr_m forcing SUM aggregation on Temp data ---
#' # results_d_t_sum <- clim_corr_m(growth, climate_d_t, daily_agg_fun = "sum")
#' # print(results_d_t_sum) # Compare correlations
#'
clim_corr_m <- function(growth_data,
                        climate_data,
                        growth_var = "res",
                        climate_var_name = "climate",
                        window_width = 12,
                        start_year = NULL,
                        end_year = NULL,
                        cor_method = "pearson",
                        daily_agg_fun = "auto", # Default to auto-detect based on var_type
                        value_col = "value"     # Assume value column is named 'value' in daily input
) {

  # --- Input Validation ---
  if (!is.data.frame(growth_data) || !all(c("year", growth_var) %in% colnames(growth_data))) stop("`growth_data` must be a data frame with 'year' and '", growth_var, "' columns.", call. = FALSE)
  if (!is.data.frame(climate_data)) stop("`climate_data` must be a data frame.", call. = FALSE)
  if (!is.numeric(growth_data$year)) stop("'year' column in `growth_data` must be numeric.", call.=FALSE)
  if (!is.numeric(growth_data[[growth_var]])) stop("Growth variable column '", growth_var, "' must be numeric.", call.=FALSE)


  # --- Detect Climate Data Frequency & Variable Type ---
  frequency <- attr(climate_data, "frequency")
  variable_type <- attr(climate_data, "variable_type")
  is_daily <- FALSE
  climate_monthly_long <- NULL

  if (!is.null(frequency) && frequency == "daily") is_daily <- TRUE
  else if (!is.null(frequency) && frequency == "monthly") is_daily <- FALSE
  else { # Guess frequency if attribute missing
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
      # Auto-detect based on variable type attribute
      if (is.null(variable_type) || !(variable_type %in% c("P", "T"))) {
        warning("Cannot auto-determine aggregation function: variable_type attribute is '",
                variable_type, "' or missing. Defaulting to 'sum'. Specify 'daily_agg_fun'.", call. = FALSE)
        agg_fun_to_use <- sum
        agg_fun_name <- "sum"
      } else if (variable_type == "P") {
        agg_fun_to_use <- sum
        agg_fun_name <- "sum"
        message("Variable type is 'P', using 'sum' for daily aggregation.")
      } else if (variable_type == "T") {
        agg_fun_to_use <- mean
        agg_fun_name <- "mean"
        message("Variable type is 'T', using 'mean' for daily aggregation.")
      }
    } else if (is.function(daily_agg_fun)) {
      agg_fun_to_use <- daily_agg_fun
      agg_fun_name <- deparse(substitute(daily_agg_fun))
      message("Using custom function '", agg_fun_name, "' for daily aggregation.")
    } else if (is.character(daily_agg_fun) && daily_agg_fun == "sum") {
      agg_fun_to_use <- sum
      agg_fun_name <- "sum"
      message("Using 'sum' (specified by user) for daily aggregation.")
    } else if (is.character(daily_agg_fun) && daily_agg_fun == "mean") {
      agg_fun_to_use <- mean
      agg_fun_name <- "mean"
      message("Using 'mean' (specified by user) for daily aggregation.")
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
        filter(!is.na(.data$date) & !is.na(.data[[value_col]])) %>% # Filter NAs before grouping
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
      select(year, all_of(as.character(1:12))) %>% # Ensure only needed cols
      pivot_longer(cols = all_of(as.character(1:12)),
                   names_to = "month", values_to = "monthly_value") %>%
      mutate(month = as.integer(.data$month)) %>%
      filter(!is.na(monthly_value)) # Remove rows where monthly value was NA
  }


  # --- Filter Years (on standardized monthly long data) ---
  climate_proc <- climate_monthly_long
  growth_filt <- growth_data
  climate_proc$year <- as.integer(climate_proc$year)
  growth_filt$year <- as.integer(growth_filt$year)

  if (!is.null(start_year)) { climate_proc <- climate_proc %>% filter(.data$year >= start_year); growth_filt <- growth_filt %>% filter(.data$year >= start_year) }
  if (!is.null(end_year)) { climate_proc <- climate_proc %>% filter(.data$year <= end_year); growth_filt <- growth_filt %>% filter(.data$year <= end_year) }

  # --- Prepare Climate Data (Rolling Sum) ---
  message("Calculating rolling sum on monthly data...")
  aggregated_climate_col_name <- paste0(climate_var_name, "_", window_width, "m")

  climate_rolling <- climate_proc %>%
    arrange(.data$year, .data$month) %>%
    mutate(
      "{aggregated_climate_col_name}" := zoo::rollsumr(.data$monthly_value, k = window_width, fill = NA)
    ) %>%
    filter(!is.na(.data[[aggregated_climate_col_name]]))

  if (nrow(climate_rolling) == 0) { warning("No climate data remaining after calculating rolling sum.", call.=FALSE); return(NULL) }

  # --- Find Common Years & Merge ---
  common_years <- intersect(growth_filt$year, climate_rolling$year)
  if (length(common_years) < 3) { warning("Fewer than 3 overlapping years after processing.", call. = FALSE); return(NULL) }
  message("Found ", length(common_years), " overlapping years (", min(common_years), "-", max(common_years), ").")

  growth_final <- growth_filt %>% filter(.data$year %in% common_years) %>% select(year, all_of(growth_var))
  climate_final <- climate_rolling %>% filter(.data$year %in% common_years)
  merged_data <- inner_join(climate_final, growth_final, by = "year")

  # --- Calculate Correlations by Month ---
  message("Calculating correlations for each ending month...")
  correlation_results <- merged_data %>%
    group_by(end_month = .data$month) %>%
    summarise(
      test_result = list(tryCatch(
        stats::cor.test(x = .data[[aggregated_climate_col_name]], y = .data[[growth_var]], method = cor_method, use = "complete.obs"),
        error = function(e) { warning("Correlation failed for month ", unique(.data$end_month), ": ", e$message, call. = FALSE); return(NULL) }
      )),
      .groups = 'drop'
    ) %>%
    mutate(
      correlation = sapply(.data$test_result, function(res) if(!is.null(res)) res$estimate else NA_real_),
      p_value = sapply(.data$test_result, function(res) if(!is.null(res)) res$p.value else NA_real_),
      n_overlap = sapply(.data$end_month, function(m) {
        sub_data <- merged_data %>% filter(.data$month == m)
        sum(stats::complete.cases(sub_data[[aggregated_climate_col_name]], sub_data[[growth_var]]))
      })
    ) %>%
    mutate(month_abbr = month.abb[.data$end_month], climate_var = aggregated_climate_col_name) %>%
    select(end_month, month_abbr, correlation, p_value, n_overlap, climate_var) %>%
    arrange(end_month)

  if(nrow(correlation_results) == 0 || all(is.na(correlation_results$correlation))) {
    warning("Correlation calculation failed for all months.", call. = FALSE); return(NULL)
  }

  message("Correlation calculation complete.")
  return(correlation_results)
}
