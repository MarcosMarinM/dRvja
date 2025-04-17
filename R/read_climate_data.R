#' Read and Process Climate Data File (Detecting Variable Type)
#'
#' Reads climate data, attempts auto-detection of format (daily/monthly) and
#' variable type (Precipitation/Temperature) based on columns, headers, and filename.
#' Standardizes column names and types while preserving the original temporal frequency.
#' Returns a data frame with attributes "frequency" ('daily'/'monthly') and
#' "variable_type" ('P'/'T'/'unknown'/'ambiguous').
#'
#' @details
#' Standardized Output Columns & Frequency:
#' \itemize{
#'   \item Monthly Wide (13 col): `year`(int), `1`..`12`(num). Freq="monthly".
#'   \item Daily YMDV (4 col): `year`(int), `month`(int), `day`(int), `value`(num). Freq="daily".
#'   \item Daily Date-Value (2 col): `date`(Date), `value`(num). Freq="daily".
#'   \item Monthly YM-Value (2 col): `year`(int), `month`(int), `value`(num). Freq="monthly".
#' }
#' Use `attr(result, "frequency")` and `attr(result, "variable_type")`.
#'
#' @param file_path Path to the climate data file.
#' @param var_type Specify variable type: 'auto' (default), 'P', 'T'.
#' @param file_format Force format: 'auto'(default), 'monthly_wide', 'daily_ymd_v', 'daily_date_v'.
#' @param date_format For 'daily_date_v', specify date format string.
#' @param daily_col_names Col names for 2-col daily if no header. Default: `c("date", "value")`.
#' @param value_col_name Name for value column in daily/monthly_long. Default: "value".
#' @param na_strings Strings for NA. Default: `c("NA", "N/A", "", "-999", "-9999", "-999.9")`. (Removed " ")
#' @param skip Lines to skip. Default: 0.
#' @param verbose Print messages? Default: TRUE.
#'
#' @return A data frame with standardized columns and attributes "frequency" and "variable_type", or NULL.
#' @export
#' @importFrom data.table fread := .SD as.data.table copy setnames
#' @importFrom lubridate parse_date_time year month day ymd floor_date is.Date
#' @importFrom tools file_path_sans_ext
#' @import dplyr
#' @import tidyr
#' @importFrom utils head tail type.convert
#'
#' @examples
#' # (Examples remain the same as before)
#' # --- Setup ---
#' make_temp_file <- function(data, fname = tempfile(fileext = ".txt"), ...) {
#'   write.table(data, fname, sep = "\t", row.names = FALSE, quote=FALSE, ...); return(fname)
#' }
#' # --- Test Files ---
#' mw_data <- data.frame(year=1990); for(m in 1:12) mw_data[[as.character(m)]]<-1; f1 <- make_temp_file(mw_data)
#' r1 <- read_climate_data(f1); print(str(r1)); print(attributes(r1)); unlink(f1)
#' dy_data <- data.frame(Y=1990,M=1,D=1:2,V=c(1,2)); f2 <- make_temp_file(dy_data)
#' r2 <- read_climate_data(f2); print(str(r2)); print(attributes(r2)); unlink(f2)
#' dd_data <- data.frame(Date="19900101", Val=3); f3 <- make_temp_file(dd_data, col.names=FALSE)
#' r3 <- read_climate_data(f3); print(str(r3)); print(attributes(r3)); unlink(f3)
#' my_data <- data.frame(YM="199001", Tavg=4); f4 <- make_temp_file(my_data, "monthly_Tavg.txt")
#' r4 <- read_climate_data(f4); print(str(r4)); print(attributes(r4)); unlink(f4)
#'
read_climate_data <- function(file_path,
                              var_type = "auto",
                              file_format = "auto",
                              date_format = NULL,
                              daily_col_names = c("date", "value"),
                              value_col_name = "value",
                              na_strings = c("NA", "N/A", "", "-999", "-9999", "-999.9"), # Removed " "
                              skip = 0,
                              verbose = TRUE) {

  # --- Dependency Checks & File Exists ---
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' required.", call. = FALSE)
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package 'lubridate' required.", call. = FALSE)
  if (!requireNamespace("tools", quietly = TRUE)) stop("Package 'tools' required.", call. = FALSE)
  if (!file.exists(file_path)) { warning("File not found: ", file_path, call. = FALSE); return(NULL) }

  if(verbose) message("\n--- Processing climate file: ", basename(file_path), " ---")

  # --- Reading with fread ---
  header_detected <- TRUE
  raw_dt <- tryCatch({
    data.table::fread(file_path, na.strings = na_strings, skip = skip, header = TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE, data.table = TRUE,
                      sep = "auto", fill = TRUE, blank.lines.skip = TRUE)
  }, error = function(e_head) {
    if(verbose) message("Info: Could not read with header, trying without header.")
    header_detected <<- FALSE
    tryCatch({
      data.table::fread(file_path, na.strings = na_strings, skip = skip, header = FALSE,
                        check.names = FALSE, stringsAsFactors = FALSE, data.table = TRUE,
                        sep = "auto", fill = TRUE, blank.lines.skip = TRUE)
    }, error = function(e_nohead){
      warning("Error reading file: ", basename(file_path), ". Error: ", e_nohead$message, call. = FALSE)
      return(NULL) })
  })

  if (is.null(raw_dt) || nrow(raw_dt) == 0 || ncol(raw_dt) == 0) {
    warning("File is empty or unreadable: ", basename(file_path), call. = FALSE); return(NULL)
  }

  # --- Basic Cleaning ---
  raw_dt <- raw_dt[, colSums(is.na(raw_dt)) < nrow(raw_dt), with = FALSE]
  raw_dt <- raw_dt[rowSums(is.na(raw_dt)) < ncol(raw_dt), ]
  if (nrow(raw_dt) == 0 || ncol(raw_dt) == 0) {
    warning("No valid data after initial NA cleaning: ", basename(file_path), call. = FALSE); return(NULL)
  }

  num_cols <- ncol(raw_dt)
  col_names_orig <- colnames(raw_dt)
  final_df <- NULL; frequency <- "unknown"; variable_type_detected <- "unknown"; input_format_detected <- "unknown"

  # ================================
  # --- Format Detection/Routing ---
  # ================================
  if(file_format != "auto") {
    if(verbose) message("User specified format: ", file_format)
    input_format_detected <- file_format
  } else {
    if (num_cols == 13) input_format_detected <- "monthly_wide"
    else if (num_cols == 4) input_format_detected <- "daily_ymd_v"
    else if (num_cols == 2) input_format_detected <- "daily_or_monthly_2col"
    else input_format_detected <- "ambiguous_cols"
    if(verbose) message("Auto-detected potential format based on ", num_cols, " columns: ", input_format_detected)
  }

  # ==============================
  # --- Processing by Format ---
  # ==============================
  dt_copy <- data.table::copy(raw_dt)

  # == Format 1: Monthly Wide (13 columns) ==
  if (input_format_detected == "monthly_wide" && num_cols == 13) {
    # (Keep previous logic for this format)
    data.table::setnames(dt_copy, c("year", as.character(1:12)), skip_absent=TRUE)
    dt_copy[, year := suppressWarnings(as.integer(year))]
    for(m_col in as.character(1:12)) { if (m_col %in% names(dt_copy)) dt_copy[, (m_col) := suppressWarnings(as.numeric(utils::type.convert(as.character(.SD[[m_col]]), as.is=TRUE)))] }
    dt_copy <- dt_copy[!is.na(year) & rowSums(is.na(dt_copy[, intersect(as.character(1:12), names(dt_copy)), with=FALSE])) < 12]
    if (nrow(dt_copy) > 0 && is.integer(dt_copy$year) && all(sapply(dt_copy[, intersect(as.character(1:12), names(dt_copy)), with=FALSE], is.numeric))) {
      if(verbose) message("Processing as: Monthly Wide (Year, 1-12).")
      final_df <- as.data.frame(dt_copy)
      frequency <- "monthly"
    } else if (input_format_detected == "monthly_wide") { warning("13 columns found, but failed interpretation as numeric Year + 12 months.", call. = FALSE) }
  }

  # == Format 2: Daily Year-Month-Day-Value (4 columns) ==
  else if (input_format_detected == "daily_ymd_v" && num_cols == 4) {
    # (Keep previous logic for this format)
    temp_dt <- dt_copy; temp_dt[, (colnames(temp_dt)) := lapply(.SD, function(x) suppressWarnings(utils::type.convert(as.character(x), as.is=TRUE)))]
    col1_num = is.numeric(temp_dt[[1]]); col2_num = is.numeric(temp_dt[[2]]); col3_num = is.numeric(temp_dt[[3]]); col4_num = is.numeric(temp_dt[[4]])
    col1_yr = col1_num && !all(is.na(temp_dt[[1]])) && all(temp_dt[[1]] > 1000 & temp_dt[[1]] < 3000, na.rm=T);
    col2_mo = col2_num && !all(is.na(temp_dt[[2]])) && all(temp_dt[[2]] >= 1 & temp_dt[[2]] <= 12, na.rm=T);
    col3_dy = col3_num && !all(is.na(temp_dt[[3]])) && all(temp_dt[[3]] >= 1 & temp_dt[[3]] <= 31, na.rm=T)
    if (col1_yr && col2_mo && col3_dy && col4_num) {
      if(verbose) message("Processing as: Daily YMDV (Year, Month, Day, Value).")
      data.table::setnames(temp_dt, c("year", "month", "day", value_col_name))
      temp_dt[, year := as.integer(year)]; temp_dt[, month := as.integer(month)]; temp_dt[, day := as.integer(day)]; temp_dt[, (value_col_name) := as.numeric(get(value_col_name))]
      final_df <- as.data.frame(temp_dt[!is.na(year)&!is.na(month)&!is.na(day)&!is.na(get(value_col_name)), c("year", "month", "day", value_col_name), with=FALSE])
      frequency <- "daily"
    } else if (input_format_detected == "daily_ymd_v") { warning("4 columns found, but failed interpretation as numeric Y, M, D, Value.", call. = FALSE) }
  }

  # == Format 3: Daily or Monthly (2 columns) - v5.2 REVISION ==
  else if (input_format_detected %in% c("daily_date_v", "monthly_ym_v", "daily_or_monthly_2col") && num_cols == 2) {
    processed_daily <- FALSE
    processed_monthly_ym <- FALSE
    # Always make a fresh copy for attempts
    temp_dt <- data.table::copy(dt_copy)
    data.table::setnames(temp_dt, c("col1", "col2")) # Use generic internal names

    # --- 3a. Attempt Daily Date-Value Interpretation ---
    if(input_format_detected %in% c("daily_date_v", "daily_or_monthly_2col")) {
      if(verbose) message("Attempting interpretation as Daily Date-Value (2 cols)...")
      try({
        original_dates <- as.character(temp_dt$col1)
        parsed_dates <- NULL
        # Expanded format list, prioritizing unambiguous ones
        common_daily_formats <- c("Ymd", "Y-m-d", "Y/m/d", "dmY", "d-m-Y", "d/m/Y", "mdY", "m-d-Y", "m/d/Y")

        if (!is.null(date_format)) parsed_dates <- suppressWarnings(lubridate::parse_date_time(original_dates, orders = date_format, quiet = TRUE))
        else parsed_dates <- suppressWarnings(lubridate::parse_date_time(original_dates, orders = common_daily_formats, quiet = TRUE))

        if (!is.null(parsed_dates) && !all(is.na(parsed_dates))) {
          temp_dt[, date := parsed_dates]
          # Ensure value column (col2) is numeric *before* filtering
          temp_dt[, value_num := suppressWarnings(as.numeric(utils::type.convert(as.character(col2), as.is=TRUE)))]
          temp_dt_filt <- temp_dt[!is.na(date) & !is.na(value_num)]

          # Check if *any* rows remain after filtering
          if(nrow(temp_dt_filt) > 0) {
            if(verbose) message("Processed as: Daily Date-Value.")
            data.table::setnames(temp_dt_filt, c("col1", "col2", "date", value_col_name)) # Rename final value column
            final_df <- as.data.frame(temp_dt_filt[, c("date", value_col_name), with=FALSE])
            frequency <- "daily"
            processed_daily <- TRUE
          }
        }
      }, silent = TRUE)
      if(!processed_daily && verbose && input_format_detected != "monthly_ym_v") message("Failed interpretation as Daily Date-Value.")
    } # End daily attempt

    # --- 3b. Attempt Monthly YM-Value Interpretation (if daily failed or user forced monthly) ---
    # Note: Forcing monthly_ym_v is less direct now, auto-detection relies on daily failing
    if (!processed_daily && input_format_detected %in% c("daily_or_monthly_2col", "monthly_ym_v")) {
      if(verbose) message("Attempting interpretation as Monthly YearMonth-Value (2 cols)...")
      temp_dt_monthly <- data.table::copy(raw_dt) # Use fresh copy
      data.table::setnames(temp_dt_monthly, c("col1", "col2"))

      original_ym_vals <- as.character(temp_dt_monthly$col1)
      parsed_ym_dates <- NULL
      common_monthly_formats <- c("Ym","Y-m","Y/m") # Keep this specific
      parsed_ym_dates <- suppressWarnings(lubridate::parse_date_time(original_ym_vals, orders = common_monthly_formats, quiet = TRUE))

      if (!is.null(parsed_ym_dates) && !all(is.na(parsed_ym_dates))) {
        temp_dt_monthly[, date := lubridate::floor_date(parsed_ym_dates, "month")]
        temp_dt_monthly[, value_num := suppressWarnings(as.numeric(utils::type.convert(as.character(col2), as.is=TRUE)))]
        temp_dt_monthly <- temp_dt_monthly[!is.na(date) & !is.na(value_num)]

        if(nrow(temp_dt_monthly) > 0) {
          if(verbose) message("Processed as: Monthly YearMonth-Value.")
          temp_dt_monthly[, year := lubridate::year(date)]
          temp_dt_monthly[, month := lubridate::month(date)]
          temp_dt_agg <- temp_dt_monthly[, .(value = data.table::first(value_num)), keyby = .(year, month)] # Use explicit data.table::first
          data.table::setnames(temp_dt_agg, "value", value_col_name)
          # Pivot to wide format here
          final_df_monthly <- temp_dt_agg %>%
            tidyr::pivot_wider(names_from = month, values_from = all_of(value_col_name), names_prefix = "") %>%
            as.data.frame()
          final_df <- final_df_monthly
          frequency <- "monthly"
          processed_monthly_ym <- TRUE
        }
      }
      if(!processed_monthly_ym && verbose) message("Failed interpretation as Monthly YearMonth-Value.")
    } # End monthly attempt
  } # End 2-column handling


  # --- Handle Failed / Ambiguous Cases ---
  if (is.null(final_df)) {
    failed_format_msg <- ifelse(grepl("fail$", input_format_detected), sub("_fail$","",input_format_detected), input_format_detected)
    warning("Failed to process file: ", basename(file_path),
            ". Could not validate structure for format '", failed_format_msg,
            "'. Check file or try specifying 'file_format'.", call. = FALSE)
    return(NULL)
  }

  # --- Standardize Final Monthly Output (if needed) ---
  if (frequency == "monthly" && !all(c("year", as.character(1:12)) %in% colnames(final_df))) {
    if(verbose) message("Standardizing monthly data to wide format (Year, 1-12)...")
    all_months_chr <- as.character(1:12)
    missing_months <- setdiff(all_months_chr, colnames(final_df))
    if (length(missing_months) > 0) { final_df[, missing_months] <- NA_real_ }
    final_df$year <- as.integer(final_df$year)
    final_df <- final_df %>% select(year, all_of(all_months_chr)) %>% arrange(.data$year)
  } else if (frequency == "monthly") {
    final_df$year <- as.integer(final_df$year)
    final_df <- final_df %>% select(year, all_of(as.character(1:12)))
  }


  # --- Final check for numeric month columns in monthly output ---
  if(frequency == "monthly") {
    month_cols_numeric <- sapply(final_df[, -1], is.numeric)
    if (!all(month_cols_numeric)) {
      warning("Some monthly columns ('1'..'12') contain non-numeric data after final processing: ", basename(file_path), call. = FALSE)
      final_df <- final_df %>% mutate(across(all_of(all_months_chr), ~suppressWarnings(as.numeric(as.character(.)))))
    }
  }

  # =================================
  # --- Variable Type Detection ---
  # =================================
  var_type_clean <- toupper(trimws(var_type))
  if (var_type_clean %in% c("P", "PRECIPITATION", "PRECIP")) {
    variable_type_detected <- "P"
    if(verbose) message("Variable type specified by user: Precipitation (P)")
  } else if (var_type_clean %in% c("T", "TEMPERATURE", "TEMP")) {
    variable_type_detected <- "T"
    if(verbose) message("Variable type specified by user: Temperature (T)")
  } else if (var_type_clean == "AUTO") {
    search_string <- ""; if(header_detected) search_string <- paste(tolower(col_names_orig), collapse=" ")
    search_string <- paste(search_string, tolower(tools::file_path_sans_ext(basename(file_path))))
    p_keywords <- c("prec", "rain", "prcp", "rr"); t_keywords <- c("temp", "tavg", "tmax", "tmin", "tas")
    found_p <- any(sapply(p_keywords, function(k) grepl(k, search_string, ignore.case = TRUE)))
    found_t <- any(sapply(t_keywords, function(k) grepl(k, search_string, ignore.case = TRUE)))
    if (found_p && !found_t) variable_type_detected <- "P"
    else if (found_t && !found_p) variable_type_detected <- "T"
    else if (found_p && found_t) variable_type_detected <- "ambiguous"
    else variable_type_detected <- "unknown"

    if(verbose){
      if(variable_type_detected == "P") message("Auto-detected variable type: Precipitation (P) based on keywords.")
      else if(variable_type_detected == "T") message("Auto-detected variable type: Temperature (T) based on keywords.")
      else if(variable_type_detected == "ambiguous") message("Auto-detection ambiguous (found P and T keywords). Type set to 'ambiguous'.")
      else message("Auto-detection failed (no keywords found). Type set to 'unknown'.")
    }
    if(variable_type_detected %in% c("ambiguous", "unknown")) {
      warning("Could not reliably auto-detect variable type (P/T). Specify 'var_type'. File: ", basename(file_path), call.=FALSE)
    }
  } else {
    warning("Invalid 'var_type' specified: '", var_type,"'. Using 'unknown'. Specify 'auto', 'P', or 'T'.", call. = FALSE)
    variable_type_detected <- "unknown"
  }

  # --- Add Attributes and Return ---
  if(!is.data.frame(final_df)) final_df <- as.data.frame(final_df)
  attr(final_df, "frequency") <- frequency
  attr(final_df, "variable_type") <- variable_type_detected
  if(verbose) message("File processing successful. Final Frequency: ", frequency, ", Variable Type: ", variable_type_detected)
  return(final_df)
}
