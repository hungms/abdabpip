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