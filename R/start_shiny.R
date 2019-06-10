#' @export
start_shiny <- function() {
  appDir <- system.file("shiny_examples", "my_app", package = "BGGM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
