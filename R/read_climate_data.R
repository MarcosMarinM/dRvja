#' Read and Process Climate Data File (Monthly or Daily)
#'
#' Reads a climate data file (e.g., precipitation) and standardizes it
#' into a monthly format with columns 'year' and '1' through '12'.
#' Handles both 13-column (year + 12 months) and 2-column (daily: date, value)
#' formats, commonly found in climate datasets. Assumes tab-separated values.
#'
#' @param file_path Character string. The full path to the climate data file (.txt).
#' @param daily_date_format Character string. The format of the date column if the
#'   input is daily data. Defaults to "%Y%m%d".
#' @param daily_col_names Character vector of length 2. Expected column names
#'   for daily data if no header is detected or standard names ('date', 'value')
#'   are not found. Defaults to `c("date", "value")`. The second element will be
#'   treated as the climate variable.
#' @param na_strings Character vector. Strings to be interpreted as NA values
#'   when reading the file. Defaults to `c("NA", "", " ", "-999", "-9999")`.
#'
#' @return A data frame with columns 'year' and '1' through '12' containing
#'   the monthly climate values, or NULL if the file cannot be processed.
#'   Issues warnings for detected problems.
#' @export
#' @importFrom utils read.table head tail
#' @importFrom readr read_tsv
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' # --- Example with a temporary monthly file ---
#' monthly_data <- data.frame(
#'   year = 1990:1992,
#'   Jan = runif(3, 0, 50), Feb = runif(3, 0, 50), Mar = runif(3, 0, 50),
#'   Apr = runif(3, 0, 50), May = runif(3, 0, 50), Jun = runif(3, 0, 50),
#'   Jul = runif(3, 0, 50), Aug = runif(3, 0, 50), Sep = runif(3, 0, 50),
#'   Oct = runif(3, 0, 50), Nov = runif(3, 0, 50), Dec = runif(3, 0, 50)
#' )
#' # Fix column names to be 1-12 for the function
#' colnames(monthly_data) <- c("year", 1:12)
#' temp_monthly_file <- tempfile(fileext = ".txt")
#' write.table(monthly_data, temp_monthly_file, sep = "\t", row.names = FALSE, quote = FALSE)
#' processed_monthly <- read_climate_data(temp_monthly_file)
#' print("Processed Monthly:")
#' print(head(processed_monthly))
#' unlink(temp_monthly_file)
#'
#' # --- Example with a temporary daily file ---
#' start_date <- as.Date("1990-01-01")
#' end_date <- as.Date("1990-03-31")
#' daily_dates <- seq(start_date, end_date, by = "day")
#' daily_data <- data.frame(
#'   date = format(daily_dates, "%Y%m%d"),
#'   precipitation = runif(length(daily_dates), 0, 5)
#' )
#' temp_daily_file <- tempfile(fileext = ".txt")
#' # Write *with* header for this example
#' write.table(daily_data, temp_daily_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
#' processed_daily <- read_climate_data(temp_daily_file)
#' print("Processed Daily:")
#' print(head(processed_daily))
#' unlink(temp_daily_file)
#'
read_climate_data <- function(file_path,
                              daily_date_format = "%Y%m%d",
                              daily_col_names = c("date", "value"),
                              na_strings = c("NA", "", " ", "-999", "-9999")) {

  if (!file.exists(file_path)) {
    warning("File not found: ", file_path, call. = FALSE)
    return(NULL)
  }

  message("\n--- Processing climate file: ", basename(file_path), " ---")

  # Try reading with read.table, suppress warnings initially
  precip_raw <- suppressWarnings(tryCatch({
    utils::read.table(file_path, header = FALSE, sep = "\t",
                      na.strings = na_strings, stringsAsFactors = FALSE,
                      comment.char = "", fill = TRUE, blank.lines.skip = TRUE)
  }, error = function(e) {
    warning("Error reading file: ", basename(file_path), ". Error: ", e$message, call. = FALSE)
    return(NULL)
  }))

  if (is.null(precip_raw) || nrow(precip_raw) == 0) {
    warning("File is empty or could not be read: ", basename(file_path), call. = FALSE)
    return(NULL)
  }

  # Remove rows that are entirely NA
  precip_raw <- precip_raw[rowSums(is.na(precip_raw)) != ncol(precip_raw), , drop = FALSE]
  if (nrow(precip_raw) == 0) {
    warning("File contains no valid data after removing empty rows: ", basename(file_path), call. = FALSE)
    return(NULL)
  }

  num_cols <- ncol(precip_raw)
  processed_df <- NULL

  # --- Monthly Format (13 columns) ---
  if (num_cols == 13) {
    message("Detected monthly format (13 columns).")
    precip_df <- precip_raw
    colnames(precip_df) <- c("year", as.character(1:12))

    # Convert to numeric, handling potential errors
    precip_df <- precip_df %>%
      mutate(across(everything(), ~suppressWarnings(as.numeric(as.character(.))))) %>%
      filter(!is.na(.data$year)) # Remove rows where year conversion failed

    if (nrow(precip_df) == 0) {
      warning("No valid years found in monthly format: ", basename(file_path), call. = FALSE)
      return(NULL)
    }
    processed_df <- precip_df

    # --- Daily Format (2 columns) ---
  } else if (num_cols == 2) {
    message("Detected daily format (2 columns). Converting to monthly...")

    # Try to intelligently detect header based on first row content
    first_row_char <- as.character(precip_raw[1, ])
    # Simple check: if first row looks non-numeric (or matches common headers) assume header
    has_header <- !all(suppressWarnings(!is.na(as.numeric(first_row_char)))) ||
      any(grepl("date|fecha|year|a\u00f1o|anjo|anyo|ano|prec|temp|value|precipitation|temperature|tmean|mean|tmax|tmin|t_mean|pp|rr", first_row_char, ignore.case = TRUE))

    if (has_header) {
      message("Attempting to read with header.")
      # Re-read with header=TRUE
      precip_daily <- suppressWarnings(tryCatch({
        utils::read.table(file_path, header = TRUE, sep = "\t",
                          na.strings = na_strings, stringsAsFactors = FALSE,
                          comment.char = "", check.names = FALSE) # check.names=F to keep original names
      }, error = function(e) NULL)) # Return NULL on error

      if(is.null(precip_daily) || ncol(precip_daily) != 2) {
        warning("Failed to re-read with header or unexpected number of columns. Assuming no header.", call. = FALSE)
        has_header <- FALSE # Fallback to no header
        precip_daily <- precip_raw # Use original read
        colnames(precip_daily) <- daily_col_names
      } else {
        colnames(precip_daily) <- tolower(colnames(precip_daily))
        # Check if standard names are present, otherwise use provided default
        if (!all(daily_col_names %in% colnames(precip_daily))) {
          warning("Column names '", paste(daily_col_names, collapse="' and '"),
                  "' not found in header. Using columns 1 and 2 as date and value.", call.=FALSE)
          colnames(precip_daily) <- daily_col_names
        } else {
          # Ensure correct order if names are standard but possibly swapped
          precip_daily <- precip_daily[, daily_col_names]
        }
      }
    } else {
      message("No header detected or assumed. Using default names: '", daily_col_names[1], "', '", daily_col_names[2], "'.")
      precip_daily <- precip_raw
      colnames(precip_daily) <- daily_col_names
    }

    # Rename consistently for processing
    names(precip_daily) <- c("date_col", "value_col")

    # Process dates
    precip_daily$date_col_str <- as.character(precip_daily$date_col)
    # Basic check for expected format (can be enhanced)
    if (!all(grepl("^\\d+$", precip_daily$date_col_str))) {
      warning("Non-numeric characters detected in date column. Attempting to filter.", call. = FALSE)
      precip_daily <- precip_daily[grepl("^\\d+$", precip_daily$date_col_str), ]
    }
    precip_daily$date <- as.Date(precip_daily$date_col_str, format = daily_date_format)

    # Process values
    precip_daily$value <- suppressWarnings(as.numeric(as.character(precip_daily$value_col)))

    # Filter valid rows
    precip_daily <- precip_daily %>% filter(!is.na(.data$date) & !is.na(.data$value))

    if (nrow(precip_daily) == 0) {
      warning("No valid daily data rows remain after conversion attempts in: ", basename(file_path), call. = FALSE)
      return(NULL)
    }

    # Aggregate to monthly
    precip_monthly_agg <- precip_daily %>%
      mutate(year = as.integer(format(.data$date, "%Y")),
             month = as.integer(format(.data$date, "%m"))) %>%
      group_by(.data$year, .data$month) %>%
      summarise(monthly_value = sum(.data$value, na.rm = TRUE), .groups = 'drop')

    if (nrow(precip_monthly_agg) == 0) {
      warning("Monthly aggregation yielded no results for: ", basename(file_path), call. = FALSE)
      return(NULL)
    }

    # Pivot to wide format (year, 1, 2, ..., 12)
    precip_df_wide <- precip_monthly_agg %>%
      tidyr::pivot_wider(names_from = .data$month,
                         values_from = .data$monthly_value,
                         names_prefix = "") # Keep names as numbers

    # Ensure all month columns 1-12 exist, fill missing with NA
    all_months_chr <- as.character(1:12)
    missing_months <- setdiff(all_months_chr, colnames(precip_df_wide))
    for (m in missing_months) {
      precip_df_wide[[m]] <- NA_real_
    }

    # Order columns: year, 1, 2, ..., 12
    processed_df <- precip_df_wide %>%
      select(year, all_of(all_months_chr)) %>%
      arrange(.data$year)

  } else {
    warning("Unrecognized format (expected 2 or 13 columns, found ", num_cols, ") for file: ", basename(file_path), call. = FALSE)
    return(NULL)
  }

  message("Climate file processing completed for: ", basename(file_path))
  return(as.data.frame(processed_df)) # Ensure final output is data.frame
}
