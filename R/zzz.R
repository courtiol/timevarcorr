#' Display welcome message
#'
#' This function should not be called by the user.
#' It displays a message when the package is being loaded.
#'
#' @param libname argument needed but automatically defined
#' @param pkgname argument needed but automatically defined
#'
#' @export
#'
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(## display message
    'timevarcorr loaded; type ?tcor for help on this package.',
    '\n'
  )
}
