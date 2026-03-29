#' Set Small Vector Elements to Zero
#'
#' A simple utility function to perform hard thresholding on a numeric vector.
#' Any element whose absolute value is less than or equal to `cut_off` is set to 0.
#'
#' @param u A numeric vector.
#' @param cut_off A single numeric value for the threshold.
#' @return The thresholded numeric vector.
#' @export
vec_cut_off <- function(u, cut_off) {
  u[abs(u) <= cut_off] <- 0
  u
}

#' Rename an R Object
#'
#' Safely renames an object within a specified environment and removes the old object.
#'
#' @param old_name The name of the object to rename, as a string.
#' @param new_name The new name for the object, as a string.
#' @param env The environment in which to look for the object (defaults to the global environment).
#' @return TRUE for success, FALSE for failure (with a warning).
#' @export
rename_object <- function(old_name, new_name, env = globalenv()) {
  if (exists(old_name, envir = env)) {
    assign(new_name, get(old_name, envir = env), envir = env)
    rm(list = old_name, envir = env)
    return(TRUE)
  } else {
    warning(paste0("Object '", old_name, "' not found!"))
    return(FALSE)
  }
}
