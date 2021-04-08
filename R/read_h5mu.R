#' @details MultiAssayExperiment-helpers
#'
#' @description Create a MultiAssayExperiment or a Seurat object from the .h5mu file
#'
#' @import hdf5r
#'
#' @exportMethod ReadH5MU
ReadH5MU <- function(file, as) {
  if (tolower(as) %in% c("mae", "multiassayexperiment")) {

    library(MultiAssayExperiment)

    # Connect to the the file
    h5 <- H5File$new(file, mode="r")

    # Check all the assays are written
    assays <- h5[['mod']]$names

    read_with_index <- function(dataset) {
      table <- dataset$read()
      dataset_attr <- h5attributes(dataset)
      
      indexcol <- "_index"
      if ("_index" %in% names(dataset_attr)) {
        indexcol <- dataset_attr$`_index`
      }

      if (indexcol %in% colnames(table)) {
        rownames(table) <- table[,indexcol,drop=TRUE]
        table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
      }
      table
    }

    # Create global colData
    metadata <- read_with_index(h5[['obs']])
    
    # Create an experiments list
    modalities <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]
      X <- view[['X']]$read()

      var <- read_with_index(view[['var']])

      obs <- read_with_index(view[['obs']])
      if (is("obs", "data.frame"))
        rownames(obs) <- paste(mod, rownames(obs), sep="-")

      se <- SummarizedExperiment(assays=SimpleList(counts=X), rowData=var, colData=obs)
      se
    })
    names(modalities) <- assays

    # Create sampleMap
    mapping <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]
      view_attr <- h5attributes(view[["obs"]])
      indexcol <- "_index"
      if ("_index" %in% names(view_attr)) {
        indexcol <- view_attr$`_index`
      }
      obs_names <- view[['obs']]$read()[,indexcol,drop=TRUE]
      sm <- data.frame(primary = obs_names,
        colname = rownames(colData(modalities[[mod]])),
        stringsAsFactors = FALSE)
      sm
    })

    obsmap <- do.call(rbind, mapping)
    obsmap["assay"] <- rep(assays, times=vapply(mapping, nrow, 1))
    obsmap <- obsmap[,c("assay", "primary", "colname")]

    # Close the connection
    h5$close_all()

    # Create a MAE object
    mae <- MultiAssayExperiment(modalities, metadata, obsmap)

    mae
  } else if (as %in% c("seurat")) {
    stop("Reading to a Seurat object is not implemented yet")
  } else {
    stop("Provide 'mae' or 'seurat' as the second argument.")
  }
}

