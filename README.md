**--> WARNING: This package is in Beta and under active development. Things may break or change without notice. Use it at your own risk and makre sure you have a backup copy of your data <--**

Please use github issues to report bugs and for feature requests



## Installation

1. install the `flowCore` package
```R
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```

2. make sure `devtools` is installed on your system.

```R
install.packages("devtools")
```

3. install `scfeatures` with the following command

```R
devtools::install_github("ParkerICI/scfeatures")
```

## Usage
This is an R package for clustering single-cell flow cytometry data and generate features to be used in mode building.

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
# and for a description of the output type. The output is saved on disk, and the function
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
    num.clusters = 200, asinh.cofactor = 5, output.type = "directory")
```

## Output

Both clustering functions ouptut two types of data:
- A summary table of per-cluster statistics
- One or more RDS (R binary format) files containing cluster memberships for every cell event

The details of the RDS output depend on the `output.type` option, please refer to the R documentation for more details. The summary table contains one row for each cluster, and one column for each channel in the original FCS files, with the table entries representing the median intensity of the channel in the corresponding cluster.

If multiple files have been pooled together this table also contains columns in the form `Marker1@A.fcs`, which contain the median expression of `Marker1`, calculated only on the cells in that cluster that came from sample `A.fcs`

### Features generation

This package also contains functions to rearrange the clustering output to calculate cluster features that can be used to build a predictive model (similar to the approach used in the Citrus package). These functions operate on the clusters summary table described above, and require data to have been pooled together before clustering (i.e. the clustering should have been run with the `cluster_fcs_files_groups` function). In other words, if you want to build a model that includes data from the four files in the example above, you need to cluster them as a single group.

The general approach for features generation for model building, is that you want to generate a table where each row represents a cluster feature (e.g. the abundance of a cluster, or the expression of a marker in a cluster), and each column represent an observation (e.g. a different sample), for which you have a categorical or continuous endpoint of interest that you want to predict using the cluster features.

This package allows you to gather data that has been collected in multiple FCS files and, after clustering, integrate all the different pieces together to generate that matrix.

The two main functions are (please refer to the R documentation for all the details):
- `get_cluster_features`: this function takes a model specification and rearranges the clustering output to produce features that are suitable for model building
- `multistep_normalize`: this function can be used to do complex normalization operations on the features

 Suppose that you have the following dataset

|file   |timepoint  |condition  |subject    |label  |tumor_size |
|-------|-----------|-----------|-----------|-------|-----------|
|A.fcs  |baseline   |stim1      |subject1   |R      |0.1        |
|B.fcs  |baseline   |stim2      |subject1   |R      |0.1        |
|C.fcs  |baseline   |unstim     |subject1   |R      |0.1        |
|D.fcs  |week8      |stim1      |subject1   |R      |1.5        |
|E.fcs  |week8      |stim2      |subject1   |R      |1.5        |
|F.fcs  |week8      |unstim     |subject1   |R      |1.5        |
|G.fcs  |baseline   |stim1      |subject2   |NR     |0.2        |
|H.fcs  |baseline   |stim2      |subject2   |NR     |0.2        |
|I.fcs  |baseline   |unstim     |subject2   |NR     |0.2        |
|L.fcs  |week8      |stim1      |subject2   |NR     |3.2        |
|M.fcs  |week8      |stim2      |subject2   |NR     |3.2        |
|N.fcs  |week8      |unstim     |subjcet2   |NR     |3.2        |

There are multiple way that you could envision leveraging this data for model construction. The two key parameters of `get_cluster_features` that allow you to specify different models are:
- `predictors`: this specifies which variables are going to be used as predictors
- `endpoint.grouping`: this specifies which variables are used to group together files that are associated with the same endpoint

A few examples should clarify how to use this function

#### Example model 1

Let's assume that you want to see if any of the data you have collected predicts the *label* variable, which is a subject-level variable. in this case you would have the following:
- `predictors`: timepoint, condition
- `endpoint.grouping`: subject

In this case the feature matrix will have the following format:
- rows: cluster x feature x timepoint x condition. So assuming a total of 200 clusters and 3 features per cluster (abundance, expression of Marker1, and expression of Marker2), there will be 200 * 3 * 2 * 3 rows in the matrix, corresponding to all the possible combinations of clusters, features, timepoints and conditions
- columns: subject. There will then be 2 columns in the matrix corresponding to the two subjects

#### Example model 1

Let's assume that you want to see if any of the data predicts the *tumor_size* variable, which is a *subject* x *timepoint*-level feature (i.e. there is a *tumor_size* measurement for each combination of *subject* and *timepoint*). In this case you would call `get_cluster_features` with the following parameters:
- `predictors`: condition
- `endpoint.grouping`: subject, timepoint

In this case the feature matrix would have the following format:
- rows: cluster x feature x condition, so 200 * 3 * 3 rows
- columns: subject, timepoint, so 2 * 3 columns


---

Copyright 2018. Parker Institute for Cancer Immunotherapy










