⚠️ `muon.r` has been deprecated: please use 

- [MuData](https://github.com/PMBio/MuDataMAE) for bioconductor objects (`MultiAssayExperiment`) or 
- [MuDataSeurat](https://github.com/PMBio/MuDataSeurat) for `Seurat` objects.

Having different libraries for different ecosystems allows to reduce the number of extra dependencies (e.g. `rhdf5` vs `hdf5r`) and to improve compatibility with other packages, such as to store matrices in a `MultiAssayExperiment` object as `DelayedArray`s on disk.
