Please use github issues to report bugs and for feature requests

## Installation

1. install the `flowCore` package
```R
# If using a version of R >= 3.6
install.packages("BiocManager")
BiocManager::install("flowCore")
# else
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```

2. make sure `devtools` is installed on your system.

```R
install.packages("devtools")
```

3. install `grappolo` with the following command

```R
devtools::install_github("ParkerICI/grappolo")
```

## Usage
This is an R package for clustering single-cell flow cytometry data and generate features to be used in mode building. The output of this clustering can be used to generate different types of visualizations using the [vite](https://github.com/ParkerICI/vite) package

The following snippets provide an example usage, documentation for all functions can be accessed directly in R.

### Clustering

Given a set of FCS files, two modes of clustering are possible:
- Each file is clustered separately
- Data from multiple files is pooled together before clustering

The choice between these two possibilities has very important implications for feature generation and model building (see below).

Assuming an input directory called `foo` that contains four files:
```
- A.fcs
- B.fcs
- C.fcs
- D.fcs
```
This is how a clustering run is setup, in case you wanted to cluster each file individually

```R
# These are the names of the columns in the FCS files that you want to use for clustering. 
# The column descriptions from the FCS files are used as name when available (corresponding
# to the $PxS FCS keyword). When descriptions are missing the channel names are used
# instead ($PxN keyword)

col.names <- c("Marker1", "Marker2", "Marker3")

# Please refer to the documentation of this function for an explanation of the parameters
# and for a description of the output. The output is saved on disk, and the function
# simply return the list of files that have been clustered
cluster_fcs_files_in_dir("foo", num.cores = 1, col.names = col.names, num.clusters = 200,
    asinh.cofactor = 5)

# You can also specify a list of files directly using the cluster_fcs_files function,
# which takes essentially the same arguments
files.list <- c("foo/A.fcs", "foo/B,fcs")
cluster_fcs_files(files.list, num.cores = 1, col.names = col.names, num.clusters = 200,
    asinh.cofactor = 5)
```

If instead you wanted to pool some files together, you would setup the run as follows

```R
# Assuming for instance that you wanted to pool A.fcs and B.fcs in group 1, and C.fcs
# and D.fcs in group2 (once again please refer to the documentation for details)
files.groups <- list(
    group1 = c("foo/A.fcs", "foo/B.fcs")
    group2 = c("foo/C.fcs", "foo/D.fcs")
)

cluster_fcs_files_groups(files.groups, num.cores = 1, col.names = col.names, 
    num.clusters = 200, asinh.cofactor = 5)
```

#### Using the GUI

A GUI is available to launch a clustering run. The GUI allows you to specify all the input options in a graphical environment, instead of having to write R code.

To launch the GUI type the following in your R console

```R
grappolo::clustering_GUI()
```

When the GUI starts you will be prompted to select a working directory. This directory must contain all the files that you want to include in the analysis. Select any file in that directory, and the directory that contains the file will be selected as working directory.


#### Output

Both clustering functions ouptut two types of data:
- A summary table of per-cluster statistics
- One or more RDS (R binary format) files containing cluster memberships for every cell event

The summary table contains one row for each cluster, and one column for each channel in the original FCS files, with the table entries representing the median intensity of the channel in the corresponding cluster.
If multiple files have been pooled together this table also contains columns in the form `Marker1@A.fcs`, which contain the median expression of `Marker1`, calculated only on the cells in that cluster that came from sample `A.fcs`

The RDS files contain R data frames, where each row represents a different cell, and the columns the intensity of different markers. A special column called `cellType` indicates cluster membership

---

Copyright 2018. Parker Institute for Cancer Immunotherapy










