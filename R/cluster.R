#' Calculate sample-level cluster statistics
#'
#'
#' @param tab A \code{data.frame}. Must contain two columns called \code{cellType} and \code{sample}, respectively indicating
#'   the cluster membership and the sample of origin
#'
#' @return Returns a \code{data.frame} with the same columns as \code{tab}, plus sample-level statistics for every column in \code{tab}.
#'   If \code{tab} contains column \code{foo}, the returned value will contain columns \code{foo@sample1}, \code{foo@sample2} etc., containing the values
#'   of foo, as calculated only on the rows of \code{tab} that belong to \code{sample1}, \code{sample2} etc.
#'
get_stats_by_sample <- function(tab) {
    tab.medians <- plyr::ddply(tab, ~cellType, plyr::colwise(median, is.numeric))
    tab.medians.by.sample <- plyr::ddply(tab, ~cellType * sample, plyr::colwise(median, is.numeric))
    pop.size <- plyr::ddply(tab, ~cellType, nrow)
    names(pop.size) <- gsub("V1", "popsize", names(pop.size))
    pop.size.by.sample <- plyr::ddply(tab, ~cellType * sample, nrow)
    names(pop.size.by.sample) <- gsub("V1", "popsize", names(pop.size.by.sample))
    tab.medians <- merge(tab.medians, pop.size, by = "cellType")
    tab.medians.by.sample <- merge(tab.medians.by.sample, pop.size.by.sample, by = c("cellType", "sample"), all.x = T)


    #Rotate the by.sample table
    temp <- reshape::melt(tab.medians.by.sample, id = c("cellType", "sample"))
    temp$variable <- paste(temp$variable, temp$sample, sep = "@")
    temp$sample <- NULL
    temp <- reshape::cast(temp, cellType~variable)


    ret <- merge(tab.medians, temp, by = "cellType", all.x = T)

    return(ret)

}

#' Clusters a table of data
#'
#' This function clusters an input \code{data.frame} using the \code{cluster::clara} function
#'
#' @param tab The input \code{data.frame}
#' @param col.names A vector specifying which columns of \code{tab} should be used for clustering
#' @param k The desired number of clusters
#' @param ... Additional arguments to be passed to \code{cluster::clara}
#'
#' @return Returns \code{tab}, with an added column called \code{groups}, indicating cluster membership
#'
cluster_data <- function(tab, col.names, k, ...) {
    m <- as.matrix(tab[, col.names])

    message("Clustering started")
    flush.console()
    groups <- cluster::clara(m, k, ...)$clustering
    message("Clustering done")
    flush.console()

    tab <- data.frame(tab, groups, check.names = F, stringsAsFactors = F)
    return(tab)
}


#' Process a group of files for clustering
#'
#' @param files A vector of strings. The first string in the vector corresponds to the name to be used for the clustering output,
#'   the remaining strings are the paths of the files that will be pooled together for clustering
#' @inheritParams cluster_fcs_files_groups
#' @inheritParams cluster_fcs_files
#'
#'
process_files_groups <- function(files, col.names, num.clusters, num.samples, asinh.cofactor, downsample.to, output.dir) {
    tab <- NULL
    orig.data <- NULL

    out.name <- files[1]
    files <- files[2:length(files)]

    message("Loading data")
    flush.console()

    for(f in files) {
        message(sprintf("Loading %s", f))
        flush.console()

        fcs.file <- flowCore::read.FCS(f)
        temp.orig.data <- flowCore::exprs(fcs.file)
        temp.tab <- convert_fcs(fcs.file, asinh.cofactor)

        if(downsample.to > 0) {
            x <- NULL
            if(nrow(temp.tab) <= downsample.to) {
                message("Number of events smaller than downsampling target, taking all events")
                x <- 1:nrow(temp.tab)
            }
            else {
                message(sprintf("Predownsampling to %d events", downsample.to))
                x <- sample(1:nrow(temp.tab), size = downsample.to)
            }
            temp.tab <- temp.tab[x,]
            temp.orig.data <- temp.orig.data[x,]
        }

        temp.tab <- as.data.frame(temp.tab, check.names = F, stringsAsFactors = F)

        temp.tab <- data.frame(temp.tab, sample = basename(f), check.names = F, stringsAsFactors = F)
        temp.orig.data <- data.frame(temp.orig.data, sample = basename(f), check.names = F, stringsAsFactors = F)
        tab <- rbind(tab, temp.tab)
        orig.data <- rbind(orig.data, temp.orig.data)
    }

    m <- grappolo:::cluster_data(tab, col.names, k = num.clusters, sampsize = min(nrow(tab), 1000), samples = num.samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- data.frame(orig.data, cellType = m[, "cellType"], check.names = FALSE, stringsAsFactors = FALSE)

    temp <- get_stats_by_sample(m)

    m <- data.frame(m, check.names = F, stringsAsFactors = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = F, check.names = F)

    write_clustering_output(temp, out.name, m, output.dir)
    return(invisible(NULL))
}

#' Process an individual file for clustering
#'
#' @param f The file path
#' @inheritParams cluster_fcs_files
#'
process_file <- function(f, col.names, num.clusters, num.samples, asinh.cofactor, output.dir) {
    fcs.file <- flowCore::read.FCS(f)
    orig.data <- flowCore::exprs(fcs.file)
    tab <- convert_fcs(fcs.file, asinh.cofactor)

    m <- grappolo:::cluster_data(tab, col.names, k = num.clusters, sampsize = min(nrow(tab), 1000), samples = num.samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- data.frame(orig.data, cellType = m[, "cellType"], check.names = FALSE, stringsAsFactors = FALSE)

    tab.medians <- plyr::ddply(m, ~cellType, plyr::colwise(median))

    pop.size <- plyr::ddply(m, ~cellType, nrow)

    temp <- data.frame(tab.medians, sample = f, popsize = pop.size[tab.medians$cellType, "V1"], check.names = F, stringsAsFactors = F)

    m <- data.frame(m, check.names = F, stringsAsFactors = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)

    write_clustering_output(temp, basename(f), m, output.dir)
    return(invisible(NULL))
}


#' Write clustering output
#'
#' @param tab.medians A \code{data.frame} containing median values for each cluster
#' @param base.name The base name for naming output files
#' @param clustered.data A \code{data.frame} containing the original data, with an exta column called \code{cellType} indicating
#'   clustering membership
#'
#'
write_clustering_output <- function(tab.medians, base.name, clustered.data, output.dir) {
    saveRDS(clustered.data, file.path(output.dir, paste(base.name, ".clustered.all_events.rds", sep = "")))
    write.table(tab.medians, file.path(output.dir, paste(base.name, ".clustered.txt", sep = "")), row.names = F, sep = "\t", quote = F)

    return(invisible(NULL))
}



#' Cluster all FCS files contained in a directory
#'
#' @param wd The directory containing the FCS files
#' @inheritDotParams cluster_fcs_files -files.list
#'
#' @return Returns the list of files that have been clustered
#'
#' @export
cluster_fcs_files_in_dir <- function(wd, ...) {
    files.list <- list.files(path = wd, pattern = "*.fcs$", full.names = TRUE, ignore.case = TRUE)
    cluster_fcs_files(files.list, ...)
    return(files.list)
}

#' Cluster FCS files
#'
#' Cluster individual FCS files, using multiple CPU cores if possible
#'
#' This function can produce two types of output:
#'  \itemize{
#'    \item{\code{"file"}}: {Two files will be written, with names derived by appending \code{".clustered.txt"} and
#'      \code{".all_events.rds"} to the file names in \code{files.list}. The \code{".clustered.txt"} file is a tab-separated
#'      table of median values for each cluster. The \code{".all_events.rds"} file is an RDS file (readable with \code{base::readRDS})
#'      containing a \code{data.frame} with all the rows in the original input file and an additional column called
#'      \code{cellType}, indicating cluster membership.
#'    }
#'    \item{\code{"directory"}}: {In addition to the \code{".clustered.txt"} file described above, this mode will create a folder
#'      called \code{"clusters_data"}. This folder will contain a sub-folder for each input file, containing separate RDS files
#'      with the data in each cluster
#'    }
#'  }
#'
#' @param files.list The files to cluster
#' @param num.cores Number of CPU cores to use
#' @param col.names A vector of column names indicating which columns should be used for clustering
#' @param num.clusters The desired number of clusters
#' @param asinh.cofactor Cofactor for asinh transformation. If this is \code{NULL} no transformation is performed (see \code{convert_fcs})
#' @param num.samples Number of samples to be used for the CLARA algorithm (see \code{cluster::clara})
#' @param output.dir The name of the output directory, it will be created if it does not exist
#' @return Returns either \code{NULL} or a \code{try-error} object if some error occurred during the computation
#' @export
cluster_fcs_files <- function(files.list, num.cores, col.names, num.clusters, asinh.cofactor, num.samples = 50, output.dir = ".") {
    if(!dir.exists(output.dir))
        dir.create(output.dir, recursive = TRUE, showWarnings = TRUE)

    parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
                       process_file, col.names = col.names, num.clusters = num.clusters,
                       num.samples = num.samples, asinh.cofactor = asinh.cofactor, output.dir = output.dir)
}


#' Pool FCS files and cluster them
#'
#' @inheritParams cluster_fcs_files_in_dir
#' @param files.list A named list of vectors detailing how the files should be pooled before clustering. Files in the same vector will
#'   be pooled together. The name of the output is going to correspond to the name of the corresponding list element.
#' @param downsample.to The number of events that should be randomly sampled from each file before pooling. If this is 0, no sampling is performed
#' @param output.dir The name of the output directory, it will be created if it does not exist
#' @inheritParams cluster_fcs_files
#'
#' @return Returns either \code{NULL} or a \code{try-error} object if some error occurred during the computation
#'
#' @export
cluster_fcs_files_groups <- function(files.list, num.cores, col.names, num.clusters, asinh.cofactor,
                                        num.samples = 50, downsample.to = 0, output.dir = ".") {

    files.list <- lapply(names(files.list), function(x) {
        c(x, files.list[[x]])
    })

    if(!dir.exists(output.dir))
        dir.create(output.dir, recursive = TRUE, showWarnings = TRUE)

    parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
                       process_files_groups, col.names = col.names, num.clusters = num.clusters, num.samples = num.samples,
                       asinh.cofactor = asinh.cofactor, downsample.to = downsample.to, output.dir = output.dir)

}



