#' Log message
#'
#' @return The message
#' @export
log_message <- function(fun, ...) {
  
  # Capture the calling function's name and arguments
  func_name <- fun

  # Get the arguments
  call <- match.call()
  args_list <- as.list(call)[-1]
  
  # Get other messages
  other_messages <- list(...)
  
  # Print the log messages
  message("LOG : ", func_name)
  message("#============================================")
  
  # Automatically print each argument's name and value
  for (arg_name in names(args_list)) {
    message(paste(arg_name, "=", args_list[[arg_name]], sep = " "))}

  # Print other messages
  for (m in other_messages) {
    message(as.character(m))}

  message("#============================================")
}