
read_with_index <- function(dataset) {
  if ("H5Group" %in% class(dataset)) {
    # Table is saved as a group rather than a dataset
    dataset_attr <- tryCatch({
      h5attributes(dataset)
    }, error = function(e) {
      list("_index" = "_index")
    })
    indexcol <- "_index"
    if ("_index" %in% names(dataset_attr)) {
      indexcol <- dataset_attr$`_index`
    }

    columns <- names(dataset)
    columns <- columns[columns != "__categories"]

    col_list <- lapply(columns, function(name) dataset[[name]]$read())
    table <- data.frame(Reduce(cbind, col_list))
    colnames(table) <- columns

    if (indexcol %in% colnames(table)) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }
  } else {
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
  }
  table
}

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

      if ("obsm" %in% names(view)) {
        obsm <- lapply(names(view[["obsm"]]), function(space) {
          view[["obsm"]][[space]]$read()
        })
        names(obsm) <- names(view[["obsm"]])
        se <- SingleCellExperiment(assays=SimpleList(counts=X), rowData=var, colData=obs, reducedDims=obsm)
      } else {
        se <- SummarizedExperiment(assays=SimpleList(counts=X), rowData=var, colData=obs)
      }

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
    
    library(Seurat)

    # Connect to the the file
    h5 <- H5File$new(file, mode="r")

    # Check all the assays are written
    assays <- h5[['mod']]$names

    # Create global metadata
    metadata <- read_with_index(h5[["obs"]])

    # Create embeddings
    if ("obsm" %in% names(h5)) {
      obsm_names <- names(h5[["obsm"]])
      obsm_names <- obsm_names[!(obsm_names %in% assays)]
      embeddings <- lapply(obsm_names, function(space) {
        emb <- t(h5[["obsm"]][[space]]$read())
        rownames(emb) <- rownames(metadata)
        emb
      })
      names(embeddings) <- obsm_names
    } else {
      embeddings <- list()
    }
    
    # mod/.../X
    modalities <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]
      X <- view[['X']]$read()

      var <- read_with_index(view[['var']])

      obs <- read_with_index(view[['obs']])
      if (is("obs", "data.frame"))
        rownames(obs) <- paste(mod, rownames(obs), sep="-")

      colnames(X) <- rownames(obs)
      rownames(X) <- rownames(var)

      assay <- CreateAssayObject(counts = X)

      assay
    })
    names(modalities) <- assays

    # mod/.../obsm
    mod_obsm <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]

      obs <- read_with_index(view[['obs']])
      if (is("obs", "data.frame"))
        rownames(obs) <- paste(mod, rownames(obs), sep="-")

      if ("obsm" %in% names(view)) {
        obsm <- lapply(names(view[["obsm"]]), function(space) {
          emb <- t(view[["obsm"]][[space]]$read())
          rownames(emb) <- rownames(obs)
          emb
        })
        # these will be added when concatenating lists
        # names(obsm) <- paste(mod, names(view[["obsm"]]), sep="_")
        names(obsm) <- names(view[["obsm"]])
        obsm
      }

      obsm
    })
    names(mod_obsm) <- assays
    mod_obsm <- do.call(c, mod_obsm)
    
    embeddings <- c(embeddings, mod_obsm)

    # Only common observations can be read
    obs_names <- Reduce(intersect, lapply(modalities, colnames))

    # Create a Seurat object
    srt <- CreateSeuratObject(modalities[[1]][,obs_names], assay=names(modalities)[1])
    for (modality in names(modalities)[2:length(modalities)]) {
      srt[[modality]] <- subset(modalities[[modality]], cells = obs_names)
    }

    # Add embeddings
    for (emb in names(embeddings)) {
      srt@reductions[[emb]] <- embeddings[[emb]][obs_names,,drop=FALSE]
    }

    # Close the connection
    h5$close_all()

    srt
  } else {
    stop("Provide 'mae' or 'seurat' as the second argument.")
  }
}

