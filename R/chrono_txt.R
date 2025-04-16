#' Write a chronology object to a tab-separated text file
#'
#' Saves a chronology object (typically the output of `dplR::chron`)
#' to a text file with columns for 'year' and selected chronology series
#' (e.g., 'std', 'res', 'samp.depth').
#'
#' @param chronology_object A data frame or object coercible to data frame,
#'   typically the output from `dplR::chron`. It must have row names
#'   that can be converted to numeric years.
#' @param fname Character string. The path and filename for the output .txt file.
#' @param chronology_types Character vector or NULL. Which columns from the
#'   chronology object to save. If NULL (default), attempts to save common
#'   types like "std", "res", "samp.depth" if they exist. Specify column
#'   names explicitly if needed (e.g., `c("std", "samp.depth")`).
#'
#' @return Invisibly returns `NULL`. Throws an error if writing fails.
#' @export
#' @importFrom utils write.table
#' @importFrom stats na.omit # Used for basic NA removal before writing
#'
#' @examples
#' if (requireNamespace("dplR", quietly = TRUE)) {
#'   data(co021, package = "dplR")
#'
#'   # Step 1: Calculate RWI (e.g., using dplR directly or your det_32)
#'   # For this example, we use dplR's detrend:
#'   rwi_co021_df <- dplR::detrend(co021, method = "Spline", nyrs = 32)
#'
#'   # Step 2: Calculate the chronology object
#'   chron_co021 <- dplR::chron(rwi_co021_df, prewhiten = TRUE)
#'
#'   # Step 3: Use chrono_txt to write the chronology to a temporary file
#'   temp_file <- tempfile(fileext = ".txt")
#'   chrono_txt(chron_co021, fname = temp_file)
#'
#'   print(paste("Chronology data written to temporary file:", temp_file))
#'   # You can inspect the file content: readLines(temp_file, n = 5)
#'
#'   # Clean up the temporary file
#'   unlink(temp_file)
#' }
#'
chrono_txt <- function(chronology_object, fname, chronology_types = NULL) {

  # --- Input Validation ---
  if (missing(fname) || !is.character(fname) || length(fname) != 1 || fname == "") {
    stop("'fname' must be provided as a single, non-empty character string path.", call. = FALSE)
  }
  if (!is.data.frame(chronology_object)) {
    chronology_df <- tryCatch(as.data.frame(chronology_object),
                              error = function(e) stop("'chronology_object' could not be coerced to a data.frame.", call.=FALSE))
  } else {
    chronology_df <- chronology_object
  }
  if (is.null(rownames(chronology_df))) {
    stop("'chronology_object' must have row names (representing years).", call. = FALSE)
  }
  years_num <- suppressWarnings(as.numeric(rownames(chronology_df)))
  if(any(is.na(years_num))) {
    stop("Row names of 'chronology_object' could not be converted to numeric years.", call. = FALSE)
  }

  # --- Select Columns ---
  available_cols <- colnames(chronology_df)
  if (is.null(chronology_types)) {
    default_types <- c("std", "res", "samp.depth") # Common columns from dplR::chron
    cols_to_save <- intersect(default_types, available_cols)
    if (length(cols_to_save) == 0) {
      warning("Default columns ('std', 'res', 'samp.depth') not found in input. Saving all available columns.")
      cols_to_save <- available_cols
    }
  } else {
    cols_to_save <- intersect(chronology_types, available_cols)
    missing_cols <- setdiff(chronology_types, available_cols)
    if(length(missing_cols) > 0) {
      warning("Specified 'chronology_types' not found: ", paste(missing_cols, collapse=", "))
    }
    if (length(cols_to_save) == 0) {
      stop("None of the specified 'chronology_types' were found as columns in the input.", call. = FALSE)
    }
  }
  message("Columns selected for saving: ", paste(cols_to_save, collapse=", "))

  # --- Prepare Data Frame for Writing ---
  output_df <- data.frame(year = years_num)
  output_df <- cbind(output_df, chronology_df[, cols_to_save, drop = FALSE])

  # Remove rows where year might be NA or all other selected columns are NA
  output_df <- output_df[stats::complete.cases(output_df$year), ] # Ensure year is not NA
  # Optional: remove rows where *all* value columns are NA? Depends on desired output
  # output_df <- output_df[rowSums(is.na(output_df[, -1, drop = FALSE])) < length(cols_to_save), ]


  # --- Write to File ---
  tryCatch(
    utils::write.table(output_df,
                       file = fname,
                       sep = "\t",
                       row.names = FALSE,
                       col.names = TRUE,
                       quote = FALSE,
                       na = "NA"), # Represent missing values as NA
    error = function(e) {
      stop("Failed to write the file '", fname, "'. Check path and permissions.\nOriginal error: ", e$message, call. = FALSE)
    }
  )

  invisible(NULL)
}
