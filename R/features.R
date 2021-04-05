# The matrix has to contain a single feature for multiple samples
transpose_feature_matrix <- function(m) {
    s <- strsplit(colnames(m), "@")
    v <- unique(sapply(s, "[", 1))
    print(v)
    stopifnot(length(v) == 1)

    m <- t(m)
    colnames(m) <- paste(paste0("cluster", 1:ncol(m)), v, sep = "_")
    row.names(m) <- sapply(s, "[", 2)

    return(m)
}


#' Melt clustering results
#'
#' This function melts pooled clustering results using \code{reshape2::melt}
#'
#' @param tab The input data (i.e. the content of a pooled \code{clustered.txt} file)
#' @param features Optional. The features to extract. If not provided all features
#'   will be extracted. Feature names in the original data are in the format
#'   \code{feature@sample}
#' @param transform.popsize Whether to transform the popsize feature into proportion,
#'   by normalizing for the total number of cells in a sample.
#'
#' @return Returns a \code{data.frame} of molten data
#' @export
melt_cluster_results <- function(tab, features = NULL, transform.popsize = TRUE) {
    col.names <- grep("@", names(tab), value = T)

    if(is.null(features))
        features <- unique(sapply(strsplit(col.names, "@"), "[[", 1))


    m <- as.matrix(tab[, col.names])


    ret <- lapply(features, function(s) {
        temp <- m[, grep(sprintf("%s@", s), colnames(m))]
        temp[is.na(temp)] <- 0

        if(s == "popsize" && transform.popsize) {
            temp <- t(temp)
            temp <- temp / rowSums(temp)
            temp <- t(temp)
        }

        temp <- transpose_feature_matrix(temp)
        temp[!is.finite(temp)] <- 0

        return(temp)

    })
    ret <- do.call(cbind, ret)

    ret <- reshape2::melt(ret)
    names(ret)[1:2] <- c("sample.id", "variable")
    return(ret)

}
