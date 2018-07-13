#' @export

clustering_GUI <- function(...) {
    shiny::runApp(appDir = file.path(system.file(package = "grappolo"), "shinyGUI"), ...)
}

get_fcs_col_names <- function(f.path) {
    fcs.file <- flowCore::read.FCS(f.path, which.lines = 1)
    params <- flowCore::pData(flowCore::parameters(fcs.file))
    ret <- as.vector(params$desc)
    
    if(any(is.na(ret))) {
        w <- is.na(ret)
        ret[w] <- as.vector(params$name[w])
    }
    
    return(ret)
}