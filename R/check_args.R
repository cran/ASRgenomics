#' Check data class
#'
#' @param data_ Data to be checked.
#' @param class_ The expected class of data_.
#'
#' @keywords internal

check.data_ <- function(data_ = NULL,
                        class_ = NULL){

  object_ <- get(data_, envir = parent.frame())

  # Test if data class is compliant.
  if( !any(class(object_) %in% class_) ) # If data has different classes: break
    stop(paste0("The ", data_, " argument should be of class(es) ",
                paste0(class_, collapse = " or ") , "."))

}

#' Check data mode
#'
#' @param data_ Data to be checked.
#' @param mode_ The expected mode of data_.
#'
#' @keywords internal

check.data.mode_ <- function(data_ = NULL,
                        mode_ = NULL){

  object_ <- get(data_, envir = parent.frame())

  # Test if data mode is compliant.
  if( !any(mode(object_) %in% mode_) ) # If data has different classes: break
    stop(paste0("The ", data_, " argument should be of mode(es) \'",
                paste0(mode_, collapse = " or ") , "'."))
}

# #' Check class of objects inside list
# #'
# #' @param data_ Data to be checked.
# #' @param class_ The expected class of data_.
# #'
# #' @keywords internal
#
# check.objects.list_ <- function(data_ = NULL,
#                         class_ = NULL){
#
#   object_ <- get(data_, envir = parent.frame())
#
#   # Test if data class is compliant.
#   if( !any(class(object_) %in% class_) ) # If data has different classes: break
#     stop(paste0("Objects inside list \'", data_, "' should be of class(es) ",
#                 paste0(class_, collapse = " or ") , "."))
#
# }

#' Check logical arguments
#'
#' @param arg_ The boolean argument to be checked.
#'
#' @keywords internal

check.logical_ <- function(arg_ = NULL){

  # Check if logical and stop if not.
  if( !is.logical( get(arg_, envir = parent.frame()) ) )
    stop(paste0("The value of \'", arg_, "' argument should be of class \'logical' (TRUE or FALSE)."))
}

#' Check string arguments
#'
#' @param data_ Data to be checked.
#' @param mandatory_ If the argument is mandatory for the analysis.
#' @param arg_ The string with the name of the function argument (e.g., \code{"gen"}).
#' @param class_ The expected class of the variable in data.
#' @param class.action_ The action to be taken if the variable has the wrong class
#' (e.g., \code{"message"}, \code{"warning"}, \code{"stop"}).
#' @param mutate_ If the argument should be mutated into the desired class.
#' @param message_ If \code{class.action_ == "message"}, write \code{message = "message"}
#' to capture upstream message command.
#'
#' @details This functions uses the \code{get} and \code{assign} which are need access to objects that
#' are one environment up on the hierarchy. The \code{envir} is set to \code{parent.frame}. If the function is looking
#' for something two or more environments up, the arguments of \code{parent.frame} have to be changed.
#'
#' @keywords internal

check.args_ <- function(data_ = NULL,
                        mandatory_ = FALSE,
                        arg_ = NULL,
                        class_ = NULL,
                        class.action_ = NULL,
                        mutate_ = NULL,
                        message_ = message){

  # Get data name.
  data.name_ <- deparse(substitute(data_))

  # Get argument name.
  arg.name_ <- deparse(substitute(arg_))

  # Get class.
  class.fun_ <- paste0("is.", class_)

  # Capture relevant info.
  variable.name_ <- get(arg.name_, envir = parent.frame())

  # Evaluate if variable is not null.
  if(!is.null(variable.name_)){

    # Perform actions for each variable passed.
    for (cur.variable.name_ in variable.name_){

      # Check cur.variable.name_ is in data (mandatory stop).
      if(!cur.variable.name_ %in% names(data_))
        stop(paste0("\'", cur.variable.name_, "' does not correspond to a variable name of \'pheno.data'."))

      # Check class of variable in data.
      if(!getFunction(class.fun_)(data_[[cur.variable.name_]])){

        # Do action unless action is message and message is FALSE and there is no mutation.
        if(!(class.action_ == "message" & isFALSE(message_)) & is.null(mutate_)){
          getFunction(class.action_)(
            paste0("Variable \'", cur.variable.name_, "' must be of class \'", class_, "'."))
        }

        # Mutate variable if requested.
        if (!is.null(mutate_)){
          data_[[cur.variable.name_]] <- getFunction(paste0("as.", class_))(data_[[cur.variable.name_]])
          assign(x = data.name_, value = data_, envir = parent.frame())

          if (message_){
            message("Coercing  \'", cur.variable.name_, "' to class \'", class_, "'.")
          }
        }
      }
    }
  }

  # Evaluate if arg is null.
  if( is.null(variable.name_) & mandatory_ )
    stop(paste0("The argument \'", arg.name_, "' is mandatory."))
}

