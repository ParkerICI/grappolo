## Usage
This is an R package for clustering single-cell flow cytometry data and generate features to be used in mode building.

The following snippets provide an example usage, documentations for all functions can be accessed directly in R

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
# and for a description of the output type. The output is saved on disk, and the function
# simply return the list of files that have been clustered
cluster_fcs_files_in_dir("foo", num.cores = 1, col.names = col.names, num.clusters = 200,
    asinh.cofactor = 5, output.type = "directory")

# You can also specify a list of files directly using the cluster_fcs_files function,
# which takes essentially the same arguments
files.list <- c("foo/A.fcs", "foo/B,fcs")
cluster_fcs_files("foo", num.cores = 1, col.names = col.names, num.clusters = 200,
    asinh.cofactor = 5, output.type = "directory")
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
    num.clusters = 200, asinh.cofactor = 5, output.type = "directory")
```





