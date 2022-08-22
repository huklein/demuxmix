# demuxmix

demuxmix is a package for demultiplexing single-cell sequencing experiments
of pooled cells labeled with barcode oligonucleotides. The package
implements methods to fit regression mixture models for a probabilistic
classification of cells, including multiplet detection. Demultiplexing
error rates can be estimated, and methods for quality control are provided.


## Installation

The package is available at [Bioconductor](https://bioconductor.org)
and can be installed via BiocManager::install:

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("demuxmix")
```
The package only needs to be installed once. Load the package into an R
session with

``` r
library(demuxmix)
```


## Vignette

Please consult the package vignette for example usages of the methods
in demuxmix. The vignette can be accessed by entering the following
command and clicking on the resulting HTML link:

``` r
browseVignettes("demuxmix")
```
