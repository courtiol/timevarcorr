#' Display welcome message
#'
#' This function should not be called by the user.
#' It displays a message when the package is being loaded.
#'
#' @param libname argument needed but automatically defined.
#' @param pkgname argument needed but automatically defined.
#'
#' @export
#' @return nothing (invisible NULL).
#'
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(## display message
    'timevarcorr loaded; type ?tcor for help on this package.',
    '\n'
  )
}


#' Determine if the package is being used by pkgdown
#'
#' This function should not be called by the user.
#' It allows to run some examples conditionally to being used by pkgdown.
#' Code copied from [pkgdown::in_pkgdown()].
#'
#' @export
#'
#' @return a logical value (`TRUE` or `FALSE`).
#'
in_pkgdown <- function() {
  identical(Sys.getenv("IN_PKGDOWN"), "true")
}
