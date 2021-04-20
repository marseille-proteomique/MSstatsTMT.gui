#' @title Run the app to use MSstatsTMT package
#'
#' @description
#'  \code{runGUI} works exactly the same as runExample from \code{\link{shiny}} package.
#'
#' @export
runGUI <- function() {

 appDir <- system.file("shiny-examples", "myapp", package = "MSstatsTMT.gui")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `MSstatsTMT.gui`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
