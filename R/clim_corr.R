#' Calculate Climate Correlation with Moving Window
#'
#' Calculates the correlation between a tree-growth variable (e.g., residual
#' chronology) and a climate variable aggregated over a specified moving window
#' (e.g., 12-month accumulated precipitation). Correlations are calculated
#' separately for windows ending in each month of the year.
#'
#' @param growth_data A data frame containing tree growth data. Must have
#'   columns 'year' and the column specified by `growth_var`.
#' @param climate_data A data frame containing monthly climate data, with columns
#'   'year' and columns '1' through '12' (output of `read_climate_data`).
#' @param growth_var Character string. The name of the column in `growth_data`
#'   containing the growth index (default: "res").
#' @param climate_var_name Character string. A base name for the climate variable
#'   used when creating the aggregated column name (e.g., "precip"). The
#'   resulting column will be named like "precip_12m". Default: "climate".
#' @param window_width Integer. The width of the moving window in months
#'   (default: 12).
#' @param start_year Optional integer. The first year to include in the analysis.
#' @param end_year Optional integer. The last year to include in the analysis.
#' @param cor_method Character string. The correlation method to use, passed
#'   to `stats::cor.test` (default: "pearson").
#'
#' @return A tidy data frame summarizing the correlation results, with columns:
#'   \describe{
#'     \item{end_month}{Integer (1-12), the ending month of the moving window.}
#'     \item{month_abbr}{Character (Jan-Dec), abbreviation of the ending month.}
#'     \item{correlation}{Numeric, the correlation coefficient (r).}
#'     \item{p_value}{Numeric, the p-value of the correlation test.}
#'     \item{n_overlap}{Integer, the number of pairs used in the correlation.}
#'     \item{climate_var}{Character, the name of the aggregated climate variable used.}
#'   }
#'   Returns NULL if calculations cannot be performed (e.g., insufficient overlap).
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom zoo rollsumr
#' @importFrom stats cor.test complete.cases
#' @importFrom utils head tail
#'
#' @examples
#' # --- Generate Sample Data ---
#' # Growth Data (like chronology output)
#' growth <- data.frame(year = 1950:2000,
#'                      res = arima.sim(model=list(ar=0.5), n=51) + rnorm(51, 1, 0.1))
#'
#' # Climate Data (like read_climate_data output)
#' climate <- data.frame(year = 1949:2001) # Need extra year for window start
#' for(m in 1:12) { climate[[as.character(m)]] <- runif(53, 10, 60) }
#'
#' # --- Run clim_corr ---
#' correlation_results <- clim_corr(
#'   growth_data = growth,
#'   climate_data = climate,
#'   growth_var = "res",
#'   climate_var_name = "precip",
#'   window_width = 12,
#'   start_year = 1960,
#'   end_year = 1990
#' )
#'
#' print(correlation_results)
#'
#' # --- Example using dplR data (requires det_32 and chrono_txt steps first) ---
#' \dontrun{
#' if(requireNamespace("dplR", quietly=TRUE) && exists("det_32") && exists("chrono_txt")) {
#'   data(co021, package="dplR")
#'   data(ca533, package="dplR") # Monthly climate data example from dplR
#'
#'   # 1. Process growth data (example steps)
#'   rwi_co021 <- det_32(co021)
#'   chron_co021 <- dplR::chron(rwi_co021)
#'
#'   # 2. Prepare climate data (dplR ca533 has year, P*, T*)
#'   # Need to convert ca533 to the format expected by clim_corr (year, 1..12)
#'   # This requires some manipulation specific to ca533 format
#'   ca533_precip <- ca533[, grep("^P", colnames(ca533))] # Select Precip cols
#'   colnames(ca533_precip) <- 1:ncol(ca533_precip) # Rename to 1..12
#'   ca533_precip$year <- rownames(ca533)
#'   ca533_precip <- ca533_precip %>% select(year, everything()) %>% mutate(year = as.integer(year))
#'
#'   # 3. Run correlation
#'   results_dplR <- clim_corr(growth_data = chron_co021, # dplR chron has 'res'
#'                            climate_data = ca533_precip,
#'                            climate_var_name = "dplR_prec",
#'                            window_width = 12)
#'   print(results_dplR)
#' }
#' }
#'
clim_corr <- function(growth_data,
                      climate_data,
                      growth_var = "res",
                      climate_var_name = "climate",
                      window_width = 12,
                      start_year = NULL,
                      end_year = NULL,
                      cor_method = "pearson") {

  # --- Input Validation ---
  if (!is.data.frame(growth_data) || !all(c("year", growth_var) %in% colnames(growth_data))) {
    stop("`growth_data` must be a data frame with 'year' and '", growth_var, "' columns.", call. = FALSE)
  }
  if (!is.data.frame(climate_data) || !("year" %in% colnames(climate_data)) || !all(as.character(1:12) %in% colnames(climate_data))) {
    stop("`climate_data` must be a data frame with 'year' and monthly columns '1' through '12'.", call. = FALSE)
  }
  if (!is.numeric(growth_data$year) || !is.numeric(climate_data$year)) {
    stop("'year' columns must be numeric.", call. = FALSE)
  }
  if (!is.numeric(window_width) || window_width < 1) {
    stop("'window_width' must be a positive integer.", call. = FALSE)
  }

  # --- Filter Years ---
  growth_filt <- growth_data
  climate_filt <- climate_data
  if (!is.null(start_year)) {
    growth_filt <- growth_filt %>% filter(.data$year >= start_year)
    climate_filt <- climate_filt %>% filter(.data$year >= start_year)
  }
  if (!is.null(end_year)) {
    growth_filt <- growth_filt %>% filter(.data$year <= end_year)
    climate_filt <- climate_filt %>% filter(.data$year <= end_year)
  }

  # --- Prepare Climate Data (Long Format + Rolling Sum) ---
  message("Preparing climate data (long format and rolling sum)...")
  aggregated_climate_col_name <- paste0(climate_var_name, "_", window_width, "m")

  climate_long <- climate_filt %>%
    tidyr::pivot_longer(cols = all_of(as.character(1:12)),
                        names_to = "month",
                        values_to = "monthly_value") %>%
    mutate(month = as.integer(.data$month)) %>%
    arrange(.data$year, .data$month) %>%
    # Calculate rolling sum - using rollsumr for align='right' by default
    # k = window_width, fill = NA
    mutate(
      "{aggregated_climate_col_name}" := zoo::rollsumr(.data$monthly_value,
                                                       k = window_width,
                                                       fill = NA)
    ) %>%
    # Remove rows where rolling sum couldn't be calculated (start of series)
    filter(!is.na(.data[[aggregated_climate_col_name]]))

  if (nrow(climate_long) == 0) {
    warning("No climate data remaining after calculating rolling sum (period too short?).", call.=FALSE)
    return(NULL)
  }

  # --- Find Common Years (after rolling sum calculation) ---
  common_years <- intersect(growth_filt$year, climate_long$year)
  if (length(common_years) < 3) {
    warning("Fewer than 3 overlapping years between growth and climate data after filtering and rolling sum calculation.", call. = FALSE)
    return(NULL)
  }
  message("Found ", length(common_years), " overlapping years for analysis (", min(common_years), "-", max(common_years), ").")

  # --- Filter data to final common years ---
  growth_final <- growth_filt %>%
    filter(.data$year %in% common_years) %>%
    select(year, all_of(growth_var))

  climate_final <- climate_long %>%
    filter(.data$year %in% common_years)

  # --- Merge Data ---
  merged_data <- inner_join(climate_final, growth_final, by = "year")

  # --- Calculate Correlations by Month ---
  message("Calculating correlations for each ending month...")
  correlation_results <- merged_data %>%
    group_by(end_month = .data$month) %>% # Group by the ending month of the window
    summarise(
      # Calculate correlation test within summarise, handle errors
      test_result = list(tryCatch(
        stats::cor.test(x = .data[[aggregated_climate_col_name]],
                        y = .data[[growth_var]],
                        method = cor_method,
                        use = "complete.obs"), # Handle potential NAs
        error = function(e) {
          warning("Could not calculate correlation for month ", unique(.data$end_month), ": ", e$message, call. = FALSE)
          return(NULL) # Return NULL if cor.test fails
        }
      )),
      .groups = 'drop' # Drop grouping after summarise
    ) %>%
    # Extract results from the list column
    mutate(
      correlation = sapply(.data$test_result, function(res) if(!is.null(res)) res$estimate else NA_real_),
      p_value = sapply(.data$test_result, function(res) if(!is.null(res)) res$p.value else NA_real_),
      # Calculate N based on actual pairs used by cor.test (approximated here)
      n_overlap = sapply(.data$end_month, function(m) {
        sub_data <- merged_data %>% filter(.data$month == m)
        sum(stats::complete.cases(sub_data[[aggregated_climate_col_name]], sub_data[[growth_var]]))
      })
    ) %>%
    # Add month abbreviation and final cleanup
    mutate(month_abbr = month.abb[.data$end_month],
           climate_var = aggregated_climate_col_name) %>%
    select(end_month, month_abbr, correlation, p_value, n_overlap, climate_var) %>%
    arrange(end_month) # Ensure output is ordered by month

  # Check if any results were generated
  if(nrow(correlation_results) == 0 || all(is.na(correlation_results$correlation))) {
    warning("Correlation calculation failed for all months.", call. = FALSE)
    return(NULL)
  }

  message("Correlation calculation complete.")
  return(correlation_results)
}
