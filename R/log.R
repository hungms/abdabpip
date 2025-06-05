#' Log message
#'
#' @return The message
#' @export
log_message <- function(...) {
  
  # Get the parent frame
  parent_call <- sys.call(-1)
  
  # Get the name of the calling function
  caller_name <- as.character(parent_call[[1]])
  
  # Get the arguments from the parent call
  parent_args <- as.list(parent_call)[-1]
  
  # Get other messages
  other_messages <- list(...)
  
  # Print the log messages
  message("LOG : ", caller_name)
  message("#============================================")
  
  # Print each argument's name and value from the calling function
  for (arg_name in names(parent_args)) {
    arg_value <- eval(parent_args[[arg_name]], envir = parent.frame())
    message(paste(arg_name, "=", arg_value, sep = " "))
  }

  # Print other messages
  for (m in other_messages) {
    message(as.character(m))
  }

  message("#============================================")
}


#' Log progress bar
#'
#' @param i The current iteration
#' @param n The total number of iterations
#' @param pb The progress bar
#' @return The progress bar
#' @export
log_progress_bar <- function(i, n, pb = NULL) {
  if (is.null(pb)) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    setTxtProgressBar(pb, i)
    return(pb)
  } else {
    setTxtProgressBar(pb, i)
    if (i == n) close(pb)
    return(invisible(NULL))
  }
}

#' Set up doSNOW cluster with progress bar
#'
#' @param n The total number of iterations
#' @param ncores The number of cores to use
#' @return A list containing the cluster and other necessary objects
#' @export
setup_snow_progress <- function(n, ncores) {
  # Check if doSNOW is available
  if (!requireNamespace("doSNOW", quietly = TRUE)) {
    warning("doSNOW package not available. Using doParallel instead.")
    doParallel::registerDoParallel(ncores)
    return(list(
      cluster = NULL,
      pb = NULL,
      using_snow = FALSE
    ))
  }
  
  # Create a cluster
  cl <- snow::makeCluster(ncores)
  
  # Create a progress bar
  pb <- txtProgressBar(max = n, style = 3)
  
  # Check doSNOW version to determine if progress parameter is supported
  dosnow_version <- packageVersion("doSNOW")
  supports_progress <- dosnow_version >= "1.0.15"  # Adjust this version as needed
  
  if (supports_progress) {
    # Define the progress function
    progress <- function() {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    }
    
    # Register with progress parameter
    doSNOW::registerDoSNOW(cl, progress = progress)
  } else {
    # Fall back to regular doSNOW without progress
    doSNOW::registerDoSNOW(cl)
    message("Your version of doSNOW doesn't support progress bars. Consider updating the package.")
  }
  
  # Return the cluster and progress bar
  return(list(
    cluster = cl,
    pb = pb,
    using_snow = TRUE,
    using_snow_progress = supports_progress
  ))
}