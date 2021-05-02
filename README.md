# muon accessibility R package

Muon for R is an experimental package to provide functionality to save R objects with multimodal data to `.h5mu` files that can be further integrated into workflows in multiple programming languages, including [`muon` Python library](https://github.com/gtca/muon), as well as to read `.h5mu` files into R objects.

`muon.r` implements a generic `WriteH5MU` function that currently works for Seurat objects (v3 and above) and MultiAssayExperiment objects. `ReadH5MU` function allows to make a Seurat or a MultiAssayExperiment object with the data from the `.h5mu` file.

## Installation

```R
remotes::install_github("gtca/muon.r")
```

## Quick start

### Seurat objects

Start with an existing dataset, e.g. a [Seurat](https://github.com/satijalab/seurat) object with CITE-seq data:

```R
library(SeuratData)
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
```
`muon.r` allows to save the object into a `.h5mu` file:

```R
library(muon)
WriteH5MU(bm, "bmcite.h5mu")
```

### MultiAssayExperiment objects

Start with a dataset, e.g. the `myMultiAssay` object from [the MAE vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html):

```R
library(muon)
WriteH5MU(myMultiAssay, "myMultiAssay.h5mu")
```

### Reading files

Data from `.h5mu` files can be read into `Seurat` or `MultiAssayExperiment` objects:

```R
bm <- ReadH5MU("bmcite.h5mu", "seurat")
```

```R
ReadH5MU("myMultiAssay.h5mu", "mae")
```


## Relevant projects

Other R packages for multimodal I/O include:

- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)
