#' Convert a flowFrame to data.frame
#'
#' @param f The \code{flowFrame} to convert
#' @param asinh.cofactor Cofactor for \code{asinh} transformation
#' @param transform.data Whether to apply \code{asinh} transformation to the data
#' @param clip.at.zero Wether to clip negative values (after transformation) at zero
#' @param compensate Wether to compensate the data using the compensation matrix embedded in the \code{flowFrame} (if any)
#'
convert_fcs <- function(f, asinh.cofactor, transform.data = T, clip.at.zero = T, compensate = T) {
    comp <- grep("SPILL", names(flowCore::description(f)), value = T)

    if(compensate && (length(comp) > 0)) {
        message("Found compensation matrix, applying...")
        comp <- flowCore::description(f)[comp][[1]]
        if(is.character(comp)) {
            comp <- strsplit(comp, ",")[[1]]
            num.channels <- as.numeric(comp[1])
            m <- matrix(nrow = num.channels, byrow = T, data = as.numeric(comp[(num.channels + 2):length(comp)]))
            colnames(m) <- comp[2:(1 + num.channels)]
            comp <- m
        }
        f <- flowCore::compensate(f, spillover = comp)
    }
    tab <- flowCore::exprs(f)
    m <- as.matrix(tab)
    if(transform.data)
        m <- asinh(m / asinh.cofactor)

    if(clip.at.zero)
        m[m < 0] <- 0

    tab <- data.frame(m, check.names = F, stringsAsFactors = F)

    colnames(tab) <- flowCore::pData(flowCore::parameters(f))$desc

    if(any(is.na(colnames(tab)))) {
        w <- is.na(colnames(tab))
        colnames(tab)[w] <- flowCore::pData(flowCore::parameters(f))$name[w]
    }

    return(tab)
}

