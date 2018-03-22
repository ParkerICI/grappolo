

get_stats_by_sample <- function(tab) {
    tab.medians <- plyr::ddply(tab, ~cellType, colwise(median, is.numeric))
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


cluster_data <- function(tab, col.names, k, algorithm = "", ...) {
    m <- as.matrix(tab[, col.names])

    groups <- clara(m, k, ...)$clustering
    print("Clustering done")
    tab <- cbind(tab, groups)
    return(tab)
}

process_files_groups <- function(files, wd, col.names, num_clusters, num_samples, asinh.cofactor, downsample.to, output_type, output.dir) {
    setwd(wd)
    tab <- NULL
    orig.data <- NULL

    for(f in files) {
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

        temp.tab <- data.frame(temp.tab, sample = f, check.names = F, stringsAsFactors = F)
        temp.orig.data <- data.frame(temp.orig.data, sample = f, check.names = F, stringsAsFactors = F)
        tab <- rbind(tab, temp.tab)
        orig.data <- rbind(orig.data, temp.orig.data)
    }

    m <- scfeatures:::cluster_data(tab, col.names, k = num_clusters, sampsize = min(nrow(tab), 1000), samples = num_samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- cbind(orig.data, cellType = m[, "cellType"])

    temp <- get_stats_by_sample(m)

    m <- data.frame(m, check.names = F, stringsAsFactors = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = F, check.names = F)

    write_clustering_output(f, temp, m, output_type, output.dir)
    #my_save(orig.data, paste(f, ".clustered.all_events.orig_data.RData", sep = ""))
}



process_file <- function(f, wd, col.names, num_clusters, num_samples, asinh.cofactor, output_type, output.dir) {
    setwd(wd)

    fcs.file <- flowCore::read.FCS(f)
    orig.data <- flowCore::exprs(fcs.file)
    tab <- convert_fcs(fcs.file, asinh.cofactor)

    m <- scfeatures:::cluster_data(tab, col.names, k = num_clusters, algorithm = "clara", sampsize = min(nrow(tab), 1000), samples = num_samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- cbind(orig.data, cellType = m[, "cellType"])

    tab.medians <- plyr::ddply(m, ~cellType, plyr::colwise(median))

    pop.size <- plyr::ddply(m, ~cellType, nrow)

    temp <- data.frame(tab.medians, sample = f, popsize = pop.size[tab.medians$cellType, "V1"], check.names = F, stringsAsFactors = F)

    m <- data.frame(m, check.names = F, stringsAsFactors = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)

    write_clustering_output(f, temp, m, output_type, output.dir)
    #my_save(orig.data, paste(f, ".clustered.all_events.orig_data.RData", sep = ""))
}

write_clustering_output <- function(base.name, tab.medians, clustered.data, output.type, output.dir) {
    if(output.type == "legacy") {
        write.table(tab.medians, paste(base.name, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
        my_save(clustered.data, paste(base.name, ".clustered.all_events.RData", sep = ""))
    }
    else if(output.type == "directory") {
        clustered.data.dir <- "clustered.data"
        txt.file.name <- paste(base.name, ".clustered.txt", sep = "")
        full.path <- file.path(output.dir, clustered.data.dir, txt.file.name)
        dir.create(full.path, recursive = T)

        write.table(tab.medians, file.path(output.dir, txt.file.name, sep = ""),
                    row.names = F, sep = "\t", quote = F)
        plyr::ddply(clustered.data, ~cellType, function(x) {
            saveRDS(x, file = file.path(full.path, sprintf("cluster_%d.RData", x$cellType[1])))
        })
    }
}

#' @export
cluster_fcs_files_in_dir <- function(wd, num.cores, col.names, num_clusters, num_samples, asinh.cofactor, output_type = "legacy") {
    files.list <- list.files(path = wd, pattern = "*.fcs$")
    output.dir <- NULL
    if(output_type == "directory")
        output.dir <- sprintf("%s.clustering_run", gsub(".fcs$", "", files.list[[1]]))

    parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
                       process_file, wd = wd, col.names = col.names, num_clusters = num_clusters,
                       num_samples = num_samples, asinh.cofactor = asinh.cofactor, output_type = output_type, output.dir = output.dir)
    return(files.list)
}

#' @export
cluster_fcs_files_groups <- function(wd, files.list, num.cores, col.names, num_clusters, num_samples,
                                     asinh.cofactor, downsample.to, output_type = "legacy", output_dir = NULL) {

    if(output_type == "directory") {
        if(is.null(output_dir))
            output_dir <- sprintf("%s.clustering_run", gsub(".fcs$", "", names(files.list)[1]))
        else
            output_dir <- sprintf("%s.clustering_run", output_dir)
    }
    parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
                       process_files_groups, wd = wd, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples,
                       asinh.cofactor = asinh.cofactor, downsample.to = downsample.to, output_type = output_type, output.dir = output_dir)

    return(files.list)
}



