#' Apply adaptive detrending to an rwl object
#'
#' This function applies different detrending methods (Spline or ModNegExp)
#' based on the length of each individual series within an rwl object, using
#' preset thresholds and spline stiffness common in certain dendrochronological
#' workflows (e.g., targeting a 32-year spline for longer series).
#'
#' @param rwl_object An object of class 'rwl' (e.g., as read by dplR::read.rwl).
#' @param threshold_long Length threshold to use spline_nyrs_long (default: 100).
#' @param threshold_medium Length threshold to use spline_nyrs_medium (default: 50).
#' @param threshold_short Length threshold to use spline_nyrs_short (default: 30).
#' @param spline_nyrs_long Years for the spline in long series (default: 32).
#' @param spline_nyrs_medium Years for the spline in medium series (default: 20).
#' @param spline_nyrs_short Years for the spline in short series (default: 10).
#' @param spline_nyrs_shortest Years for the spline in very short series if ModNegExp fails (default: 5).
#' @param try_modnegexp_short Attempt ModNegExp for short series before spline_nyrs_short (default: TRUE).
#' @param try_modnegexp_shortest Attempt ModNegExp for very short series before spline_nyrs_shortest (default: TRUE).
#'
#' @return A data.frame where each column is a detrended series (Ring Width Index - RWI),
#'   and row names correspond to the years.
#' @export
#' @importFrom dplR detrend.series read.rwl # Added read.rwl here for example usage
#' @importFrom stats sd # Example if you use sd, etc. Add other base package imports if needed
#' @importFrom utils installed.packages # Example
#'
#' @examples
#' if (requireNamespace("dplR", quietly = TRUE)) {
#'   # Load example data from dplR
#'   data(co021, package = "dplR")
#'   # Apply adaptive detrending with default parameters
#'   rwi_co021 <- det_32(co021)
#'   # View the first few rows/columns
#'   print(head(rwi_co021[, 1:5]))
#'
#'   # Example with modified parameters
#'   # rwi_co021_alt <- det_32(co021, threshold_long = 150, spline_nyrs_long = 50)
#' }
#'
det_32 <- function(rwl_object,
                   threshold_long = 100,
                   threshold_medium = 50,
                   threshold_short = 30,
                   spline_nyrs_long = 32,
                   spline_nyrs_medium = 20,
                   spline_nyrs_short = 10,
                   spline_nyrs_shortest = 5,
                   try_modnegexp_short = TRUE,
                   try_modnegexp_shortest = TRUE) {

  # --- Internal helper function ---
  # (Defined inside so it's not visible outside det_32)
  detrend_individual <- function(series_vec) {
    # Use sum(!is.na()) which is more robust than length() for vectors with NAs
    series_length <- sum(!is.na(series_vec))

    # If the series is completely empty, return it as is (or NAs of the same length)
    if (series_length == 0) {
      return(series_vec)
    }

    if (series_length >= threshold_long) {
      detrended_series <- dplR::detrend.series(y = series_vec, method = "Spline", nyrs = spline_nyrs_long)
    } else if (series_length >= threshold_medium) {
      detrended_series <- dplR::detrend.series(y = series_vec, method = "Spline", nyrs = spline_nyrs_medium)
    } else if (series_length >= threshold_short) {
      if (try_modnegexp_short) {
        detrended_series <- tryCatch(
          dplR::detrend.series(y = series_vec, method = "ModNegExp"),
          error = function(e) {
            warning(paste("ModNegExp failed for series of length", series_length, "(short threshold), using Spline(", spline_nyrs_short, "). Message:", e$message), call. = FALSE)
            dplR::detrend.series(y = series_vec, method = "Spline", nyrs = spline_nyrs_short)
          }
        )
      } else {
        detrended_series <- dplR::detrend.series(y = series_vec, method = "Spline", nyrs = spline_nyrs_short)
      }
    } else { # Very short series (< threshold_short)
      if (try_modnegexp_shortest) {
        detrended_series <- tryCatch(
          dplR::detrend.series(y = series_vec, method = "ModNegExp"),
          error = function(e) {
            warning(paste("ModNegExp failed for series of length", series_length, "(very short threshold), using Spline(", spline_nyrs_shortest, "). Message:", e$message), call. = FALSE)
            dplR::detrend.series(y = series_vec, method = "Spline", nyrs = spline_nyrs_shortest)
          }
        )
      } else {
        # Ensure nyrs is not greater than series length - 1 (or similar, check detrend.series docs)
        # Basic adjustment, might need refinement based on spline requirements
        actual_nyrs_shortest <- min(spline_nyrs_shortest, max(1, series_length - 2))
        # Check if spline is feasible (needs at least ~3 non-NA points?)
        if(actual_nyrs_shortest < 3) {
          warning(paste("Series too short (length", series_length, ") to apply spline with nyrs=", spline_nyrs_shortest, ". Returning NAs."), call. = FALSE)
          # Decide what to return in this edge case: original series, NAs, etc.
          # Returning NAs might be safer to avoid propagating potentially meaningless values
          detrended_series <- series_vec
          detrended_series[] <- NA
        } else {
          detrended_series <- dplR::detrend.series(y = series_vec, method = "Spline", nyrs = actual_nyrs_shortest)
        }
      }
    }
    return(detrended_series)
  }
  # --- End helper function ---

  # Extract original years from row names
  original_years <- rownames(rwl_object)

  # Apply the helper function to each column (series) of the rwl object
  # Use lapply on the underlying data.frame structure
  # Ensure input is treated as data.frame for lapply to work on columns
  if (!is.data.frame(rwl_object)) {
    rwl_object_df <- as.data.frame(rwl_object)
  } else {
    rwl_object_df <- rwl_object
  }
  detrended_list <- lapply(rwl_object_df, detrend_individual)

  # Convert the resulting list back to a data frame
  detrended_data <- as.data.frame(detrended_list)

  # Reassign row names (years)
  if (!is.null(original_years) && length(original_years) == nrow(detrended_data)) {
    rownames(detrended_data) <- original_years
  } else if (nrow(detrended_data) > 0) {
    warning("Original years information might be lost or incompatible with detrended data dimensions.")
    # Attempt to infer years if possible, or leave them as default sequence
  } else {
    # Handle case where detrended_data might be empty
    warning("Detrending resulted in an empty data frame.")
  }


  # Return the data frame with detrended series (RWI)
  return(detrended_data)
}
