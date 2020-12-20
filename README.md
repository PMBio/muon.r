# muon accessibility R package

Muon for R is an experimental package to provide converters that allow to save R objects with multimodal data to `.h5mu` files that can be further integrated into workflows in multiple programming languages, including [`muon` Python library](https://github.com/gtca/muon).

`rmuon` implements a generic `WriteH5MU` function that currently works for Seurat objects (v3 and above) and MultiAssayExperiment objects.

## Installation

```R
remotes::install_github("gtca/rmuon")
```

## Quick start

Start with an existing dataset, e.g. a Seurat object with CITE-seq data:

```R
library(SeuratData)
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
```
`rmuon` allows to save the object into a `.h5mu` file:

```R
library(muon)
WriteH5MU(muon, "bmcite.h5mu")
```

## Relevant projects

Other R packages for multimodal I/O include:

- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)
