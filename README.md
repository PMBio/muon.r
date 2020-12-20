# muon accessibility R package

Muon for R package (`rmuon`) is an experimental package to provide converters that allow to save R objects with multimodal data to `.h5mu` files that can be further integrated into workflows in multiple programming languages, including [`muon` Python library](https://github.com/gtca/muon).

`rmuon` implements a generic `WriteH5MU` function that currently works for Seurat objects (v3 and above) and MultiAssayExperiment objects.

## Relevant projects

Other R packages for multimodal I/O include:

- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)
