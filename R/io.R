#' Convert a flowFrame to data.frame
#'
#' This function converts a \code{flowFrame} object to an R data.frame, taking care of data
#' transformation and compensation.
#'
#' The column names of the resulting data.frame will be derived from the channel descriptions (i.e. $PxS FCS keyword)
#' when possible, and the cannel names otherwise ($PxN)
#'
#' @param f The \code{flowFrame} to convert
#' @param asinh.cofactor Cofactor for \code{asinh} transformation
#' @param transform.data Whether to apply \code{asinh} transformation to the data
#' @param clip.at.zero Wether to clip negative values (after transformation) at zero
#' @param compensate Wether to compensate the data using the compensation matrix embedded in the \code{flowFrame} (if any)
#'
#' @return Returns a \code{data.frame} corresponding to the data in \code{flowCore::exprs(f)} after compensation
#'   and transformation
#'
#' @export
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


#' Get the columns that are common to a set of input tabular files
#'
#' @param files.list A vector of input file names. If these are text files, each file should be a tab-separated table,
#'   with the first row representing column headers
#' @param file.type The type of the files
#' @return Returns a vector of column names that are present in all the files in \code{files.list}
#' @export
get_common_columns <- function(files.list, file.type = c("txt", "fcs")) {
    l <- list()
    file.type <- match.arg(file.type)

    for(f in files.list) {
        temp <- NULL

        if(file.type == "txt")
            temp <- read.table(f, header = T, sep = "\t", check.names = F, quote = "", nrows = 1)
        else if(file.type == "fcs") {
            fcs <- flowCore::read.FCS(f, which.lines = 1)
            temp <- convert_fcs(fcs, asinh.cofactor = 1, transform.data = F, compensate = F)
        }

        l <- c(l, list(names(temp)))

    }
    return(Reduce(intersect, l))
}

#' Remove empty files
#'
#' This is an utility function to remove empty FCS files, which can cause problems with downstream analysis
#'
#' @param files.list The list of FCS files to check, files with number of events \code{<= events.threshold}
#'   will be removed
#' @param events.threshold A number. FCS files with number of events \code{<=} than this threshold will be removed
#'
#' @export
remove_empty_files <- function(files.list, events.threshold = 0) {
    for(f in files.list) {
        fcs <- flowCore::read.FCS(f)
        if(dim(fcs)["events"] <= events.threshold) {
            message(sprintf("Removing %s", f))
            file.remove(f)

        }
    }
    return(invisible(NULL))
}



