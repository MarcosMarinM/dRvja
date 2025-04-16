#' Write Ring Width Index (RWI) data to a file
#'
#' Saves a data frame containing RWI series (typically the output of a
#' detrending function like `det_32`) to a file in a standard
#' dendrochronology format, commonly the Tucson (.rwl) format used for RWI.
#'
#' @param rwi_object A data.frame with years as row names and series IDs
#'   as column names, containing the RWI values. This is typically the
#'   output from a function like `det_32`.
#' @param fname Character string. The path and filename for the output file
#'   (e.g., "my_study_area_rwi.rwl").
#' @param format Character string. The format for saving the file. Passed
#'   to `dplR::write.rwl`. Common options include "tucson" (default),
#'   "compact". See `?dplR::write.rwl` for details.
#' @param ... Optional. Other arguments passed on to `dplR::write.rwl`, such
#'   as `append = FALSE` or `prec = 0.001`.
#'
#' @return Invisibly returns `NULL` upon successful completion. Throws an
#'   error if writing fails (e.g., invalid path, permissions issues).
#' @export
#' @importFrom dplR as.rwl write.rwl
#'
#' @examples
#' \dontrun{
#' # This example assumes you have run det_32 previously
#' # and have the result in 'my_rwi_data'
#'
#' # Example: Generate some dummy RWI data (replace with actual data)
#' years <- 1950:2000
#' series1 <- 1 + rnorm(length(years), 0, 0.2)
#' series2 <- 1 + rnorm(length(years), 0, 0.3)
#' my_rwi_data <- data.frame(series1, series2, row.names = years)
#'
#' # Define the output file path
#' output_file <- "example_output.rwl"
#'
#' # Write the RWI data using the function
#' write_rwi(my_rwi_data, fname = output_file)
#'
#' # Check if the file exists
#' file.exists(output_file)
#'
#' # Clean up the created file
#' unlink(output_file)
#'
#' # Example writing in compact format
#' # write_rwi(my_rwi_data, fname = "example_compact.txt", format = "compact")
#' # unlink("example_compact.txt")
#' }
#'
#' # Example using dplR data and temporary file (safer for examples)
#' if (requireNamespace("dplR", quietly = TRUE)) {
#'   data(co021, package = "dplR")
#'   # Assuming det_32 is loaded from your package (dRvja)
#'   # rwi_co021 <- det_32(co021) # Need to run load_all() first
#'   # Using dplR's detrend directly here for a self-contained example
#'   rwi_co021_df <- dplR::detrend(co021, method = "Spline", nyrs = 32)
#'
#'
#'   temp_file <- tempfile(fileext = ".rwl")
#'   write_rwi(rwi_co021_df, fname = temp_file)
#'   print(paste("RWI data written to temporary file:", temp_file))
#'   # You can check the content, e.g., head(dplR::read.rwl(temp_file))
#'   unlink(temp_file) # Clean up
#' }
#'
write_rwi <- function(rwi_object, fname, format = "tucson", ...) {

  # --- Input Validation (Basic) ---
  if (!is.data.frame(rwi_object)) {
    stop("'rwi_object' must be a data.frame.", call. = FALSE)
  }
  if (missing(fname) || !is.character(fname) || length(fname) != 1 || fname == "") {
    stop("'fname' must be provided as a single, non-empty character string path.", call. = FALSE)
  }
  # Could add more checks: does the directory exist? do rownames look like years?

  # --- Core Logic ---
  # Convert the data.frame to an rwl object using dplR::as.rwl
  # This function expects years as rownames and series as columns.
  rwl_to_write <- tryCatch(
    dplR::as.rwl(rwi_object),
    error = function(e) {
      stop("Failed to convert 'rwi_object' to an 'rwl' object using dplR::as.rwl.\nCheck if input is a data.frame with years as row names.\nOriginal error: ", e$message, call. = FALSE)
    }
  )

  # Write the rwl object to the specified file using dplR::write.rwl
  # Pass the format and any other arguments (...)
  tryCatch(
    dplR::write.rwl(rwl.df = rwl_to_write, fname = fname, format = format, ...),
    error = function(e) {
      stop("Failed to write the file '", fname, "'.\nPlease check the path, permissions, and format.\nOriginal error: ", e$message, call. = FALSE)
    }
  )

  # Return NULL invisibly on success, standard for writing functions
  invisible(NULL)
}
