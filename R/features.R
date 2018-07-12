#' Normalize numeric data in multiple steps
#'
#' This function normalizes numeric data in a table in multiple steps, based on groups defined by different combinations of categorical variables. The table
#' is assumed to be in molten format (e.g. see \code{reshape::melt}), with a single value column. This function will then proceed to normalize the data
#' based on the value identified by a categorical variable, then normalize the normalized data again by another value etc. (see below for Details)
#'
#' For this function to work the inputs need to satisfy a number of conditions
#'
#'   \itemize{
#'     \item{The input table needs to be in molten format (i.e. see \code{reshape::melt}). There need to be \code{variable} and
#'       \code{value} columns identifying the variables that will be normalized and their values}
#'     \item{The remaining columns of the input table should all be categorical variables (either strings or factors), identifying different
#'       subsets of rows}
#'     \item{At each step of the normalization the table is grouped using the \code{variable}, \code{subject.var} and all the columns
#'       in \code{names(norm.template)}. After this grouping, for every group, there can be only one row for the value of the current grouping
#'       variable that has been selected as a basis for normalization. In other words the function will not allow you to normalize a vector of values
#'       by another vector of values, it will only allow normalization of a vector by an individual number. This is done to prevent the result to depend
#'       on the ordering of the table.}
#'   }
#' An example should help clarify the working of this function. Assume you have a dataset where different variables have been measured for multiple subjects,
#' under different stimulation conditions, and at different timepoints. For each variable you want the data at each timepoint to be normalized by the value in
#' the "unstim" condition. Then you want this data to be further normalized by the value at the "baseline" timepoint. Assume \code{tab} is in molten
#' format and has the following columns
#'   \itemize{
#'     \item{\code{variable}}: identifies the variable
#'     \item{\code{value}}: the corresponding value of the variable
#'     \item{\code{timepoint}}: categorical variable that identifies the timepoint
#'     \item{\code{condition}}: categorical variable that identified the condition
#'     \item{\code{subject}}: categorical variable that identifies data belonging to the same subject (all the normalization is done within subject)
#'   }
#' To achieve the result described above, you would invoke this function as \code{multistep_normalize(tab, list(condition = "unstim", timepoint = "baseline"), "subject")}.
#' Note that the function would fail if you only specify a single variable (either \code{condition} or \code{timepoint}), because a single variable is not enough
#' to identify a single value for normalization, since you have multiple conditions for each timepoint and viceversa.
#'
#' @param tab The input \code{data.frame} See the details for assumption about its structure
#' @param norm.template A named list idenfying which categorical variables should be used to group data for normalization. The values in the list
#'   represent the value of the corresponding variable that identify the rows that are used as reference for normalization at each step. The data will be normalized
#'   in the same order specified by this list (i.e. data will be normalized according to the first variable, then again according to the second etc.)
#' @param subject.var The name of the column that identifies different subjects in \code{tab}. All normalization operations are done within the subgroups
#'   identified by this variable (i.e. data will never be normalized across subsets idenfied by different values of subject.var)
#'
#'
#'
#' @export
multistep_normalize <- function(tab, norm.template, subject.var) {
    var.names <- names(norm.template)
    var.values <- unlist(norm.template, use.names = F)

    ret <- tab

    for(i in 1:length(var.names)) {
        variable.names <- c("variable", subject.var, var.names[var.names != var.names[i]])
        mutate.s <- sprintf("value / value[%s == '%s']", var.names[i], var.values[i])
        filter.s <- sprintf("%s != '%s'", var.names[i], var.values[i])

        # Check for uniqueness of normalization values
        dplyr::group_by_(ret, .dots = variable.names) %>%
            dplyr::do({
                x <- .[[var.names[i]]]
                if(length(x[x == var.values[i]]) != 1)
                    stop("This combination of variables does not identify a single reference value for normalization")
                .
            })


        ret <- dplyr::group_by_(ret, .dots = variable.names) %>%
            dplyr::mutate_(.dots = setNames(mutate.s, "value")) %>%
            dplyr::filter_(.dots = filter.s)

    }

    return(ret)
}


#' Calculate cluster features for model building
#'
#' This function takes a clustering result and a table of sample metadata, and calculates cluster features to be used for model building
#'
#' This function is designed to work with the results of \code{cluster_data}. These results contain sample-level values for size of the clusters
#' and channel intensities. A column such as \code{foo@sample1} identifies the \code{sample1}-specific value of variable \code{foo} for each cluster. The
#' \code{metadata.tab} must contain a \code{file} column, which matches the names of the samples in \code{tab} (i.e. the part after the \code{@}, "sample1" in the
#' above example). The rest of the columns in \code{metadata.tab} represent file-level metadata, which is used to identify the data corresponding to
#' a given combination of predictors (see below)
#' An example will help clarify the working of this function. Suppose you have collected data from multiple patients at multiple timepoints and under multiple
#' stimulation conditions.
#' In this case the \code{metadata.tab} would look like this
#' \itemize{
#'   \item{\code{file}}{The names of the data files that contain data for each sample. These must match the names in the clustering results (see above)}
#'   \item{\code{timepoint}}{The timepoint information}
#'   \item{\code{condition}}{The stimulation condition}
#'   \item{\code{subject}}{The subjet each file was derived from}
#' }
#' Let's assume a few different scenarios.
#' \enumerate{
#'   \item You have subject level information (e.g. "responder" vs "non-responder") and you want to predict whether any combination of the \code{timepoint} and
#'         \code{condition} information predicts this outcome. In this case you would call the function with \code{predictors = c("condition", "timepoint")} and
#'         \code{endpoint.grouping = "sample"}. The features in the resulting output would look like \code{cluster_1_feature1_condition_timepoint}
#'   \item You have subject and timepoint level information, and you want to see if any of the stimulation conditions predicts it. In this case you would call
#'         the function with \code{predictors = c("condition")} and \code{endpoint.grouping = c("sample", "timepoint")}. The features in the resulting output
#'         would look like \code{cluster_1_feature1_condition}
#' }
#'
#' @param tab A \code{data.frame} representing clustering results, as produced by \code{cluster_data} (see Details)
#' @param metadata.tab A \code{data.frame} containing sample metadata (see Details)
#' @param features.names The name of the features in \code{tab} that are to be included in the output. These names correspond to the portion before the \code{@}
#'   in \code{names(tab)}
#' @param out.format The format of the return value, see below for detail
#' @param predictors Only used if \code{out.format == "tidy"}. Columns in \code{metadata.tab} that identify predictors.
#' @param endpoint.grouping Only used if \code{out.format == "tidy"}. Columns in \code{metadata.tab} that identify the grouping of the response variable
#'   (see Details). The combination of \code{predictors} and \code{endpoint.grouping} must uniquely identify every row in \code{metadata.tab}.
#'   The function will throw an error if this is not the case.
#'
#' @return Returns a data frame whose format depends on the value of the \code{format} parameter
#'   \itemize{
#'     \item{table}: each row corresponds to a combination of the levels of the variables specified in \code{endpoint.grouping}, and the columns are
#'     cluster features, which are combinations of the levels of the \code{predictors} for each feature specified in \code{features.names}
#'     \item{tidy}: there is a single numeric column, and all the other columns represent variables whose combinations uniquely identify each observation (i.e. each row)
#'   }
#' @export

get_cluster_features <- function(tab, metadata.tab, features.names, out.format = "table", predictors = NULL, endpoint.grouping = NULL) {
    out.format <- match.arg(out.format, c("table", "tidy"))
    m <- reshape_cluster_features(tab, features.names)

    df <- reshape::melt(m, varnames = c("file", "variable"))

    df <- merge(df, metadata.tab, by = "file")
    ret <- NULL

    if(out.format == "tidy")
        ret <- df
    else {
        if(is.null(predictors) || is.null(endpoint.grouping))
            stop("You must specifiy predictors and endpoint.grouping when out.format == 'table'")

        formula.exp <- as.formula(sprintf("%s ~ %s", paste(endpoint.grouping, collapse = "+"),
                                          paste(c("variable", predictors), collapse = "+")))

        #Stop if the combination of response grouping and predictors does not uniquely identify each file
        if(nrow(df) != nrow(unique(df[, c("variable", predictors, endpoint.grouping)])))
            stop("The combination of predictors and endpoint.grouping does not uniquely identify every row in metadata.tab")
        ret <- reshape::cast(df, formula.exp)
    }

    return(ret)
}


run_test <- function() {
    tab <- read.table("C:/Users/fgherardini/temp/standalone_citrus/data/Patient20_diseased_unstim.fcs.clustered.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)
    metadata.tab <- read.table("C:/Users/fgherardini/temp/standalone_citrus/data/test_metadata.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)

    df <- grappolo:::get_cluster_features(tab, metadata.tab, c("FunctionalMarker1", "FunctionalMarker2", "popsize"),
                                                  predictors = c("condition"), endpoint.grouping = c("day", "sample"))

}

# The matrix has to contain a single feature for multiple samples
transpose_feature_matrix <- function(m) {
    s <- strsplit(colnames(m), "@")
    v <- unique(sapply(s, "[", 1))
    stopifnot(length(v) == 1)

    m <- t(m)
    colnames(m) <- paste("cluster", 1:ncol(m), v, sep = "_")
    row.names(m) <- sapply(s, "[", 2)

    return(m)
}


reshape_cluster_features <- function(input.tab, features) {
    col.names <- sapply(features, paste, "@", sep = "")
    col.names <- paste(col.names, collapse = "|")
    col.names <- grep(col.names, names(input.tab), value = T)

    m <- as.matrix(input.tab[, col.names])


    ret <- lapply(features, function(s) {
        temp <- m[, grep(s, colnames(m))]
        temp[is.na(temp)] <- 0

        if(s == "popsize") {
            temp <- t(temp)
            temp <- temp / rowSums(temp)
            temp <- t(temp)
        }

        temp <- transpose_feature_matrix(temp)
        temp[!is.finite(temp)] <- 0

        return(temp)

    })
    ret <- do.call(cbind, ret)
    return(ret)

}
