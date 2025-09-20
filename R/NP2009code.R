#' Display the analysis code from the 2009 Nature protocols paper
#'
#' This function opens an editor displaying the analysis code of the Nature
#' Protocols 2009 paper
#'
#' The [edit()] function uses `getOption("editor")` to select
#' the editor. Use, for instance, `options(editor="emacs")` to set another
#' editor.
#'
#' @author Steffen Durinck, Wolfgang Huber
#' @seealso [edit()]
#' @keywords methods
#'
#' @examplesIf interactive()
#' NP2009code()
#'
#' @export
#' @importFrom utils edit
NP2009code <- function() {
  edit(file = system.file("scripts", "Integration-NP.R", package = "biomaRt"))
}
