
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

    col_list <- lapply(columns, function(name) {
      values <- dataset[[name]]$read()
      values_attr <- tryCatch({
        h5attributes(dataset[[name]])
      }, error = function(e) {
        list()
      })
      if (length(values_attr) > 0) {
        if ("categories" %in% names(values_attr)) {
          # Make factors out of categorical data
          ref <- values_attr$categories
          values_labels <- ref$dereference(obj = NULL)[[1]]
          values <- factor(as.integer(values), labels = values_labels$read())
        }
      }
      values
    })
    table <- data.frame(Reduce(cbind, col_list))
    colnames(table) <- columns

    if (indexcol %in% colnames(table)) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }

    # DEPRECATED:
    # For consistency with other tools, this is done via references (see above)
    # Make factors out of categorical data
    # if ("__categories" %in% names(dataset)) {
    #   cats <- dataset[["__categories"]]
    #   for (cat in names(cats)) {
    #     table[[cat]] <- factor(as.integer(table[[cat]]) + 1, labels = cats[[cat]]$read())
    #   }
    # }

    # Fix column order
    if ("column-order" %in% names(dataset_attr)) {
      ordered_columns <- dataset_attr[["column-order"]]
      # Do not consider index as a column
      ordered_columns <- ordered_columns[ordered_columns != indexcol]
      table <- table[,ordered_columns[ordered_columns %in% columns],drop=FALSE]
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

read_matrix <- function(dataset) {
  if ("data" %in% names(dataset) && "indices" %in% names(dataset) && "indptr" %in% names(dataset)) {
      i <- dataset[["indices"]]$read()
      p <- dataset[["indptr"]]$read()
      x <- dataset[["data"]]$read()
      if ("shape" %in% h5attr_names(dataset)) {
        X_dims <- h5attr(dataset, "shape")
      } else {
        X_dims <- c(max(i), max(p))
      }
      X <- Matrix(0, X_dims[1], X_dims[2])
      X@i <- i
      X@p <- p
      X@x <- x
      t(X)
    } else {
      dataset$read()
    }
}


OBSM2VARM <- list("X_pca" = "PCs", "X_mofa" = "LFs")

#' @details MultiAssayExperiment-helpers
#'
#' @description Create a MultiAssayExperiment or a Seurat object from the .h5mu file
#'
#' @import hdf5r, Matrix
#'
#' @exportMethod ReadH5MU
ReadH5MU <- function(file, as) {
  if (tolower(as) %in% c("mae", "multiassayexperiment")) {

    library(MultiAssayExperiment)

    # Connect to the the file
    h5 <- open_and_check_mudata(file)

    # Check all the assays are written
    assays <- h5[['mod']]$names

    # Create global colData
    metadata <- read_with_index(h5[['obs']])

    # Create an experiments list
    modalities <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]
      X <- read_matrix(view[['X']])

      var <- read_with_index(view[['var']])

      obs <- read_with_index(view[['obs']])
      if (is("obs", "data.frame"))
        rownames(obs) <- paste(mod, rownames(obs), sep="-")

      if ("obsm" %in% names(view)) {
        obsm <- lapply(names(view[["obsm"]]), function(space) {
          t(view[["obsm"]][[space]]$read())
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
    h5 <- open_and_check_mudata(file)

    # Check all the assays are written
    assays <- h5[['mod']]$names

    # Create global metadata
    metadata <- read_with_index(h5[["obs"]])
    ft_metadata <- read_with_index(h5[["var"]])

    # Get embeddings
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

    # Get loadings
    if ("varm" %in% names(h5)) {
      varm_names <- names(h5[["varm"]])
      varm_names <- varm_names[!(varm_names %in% assays)]
      loadings <- lapply(varm_names, function(space) {
        ld <- t(h5[["varm"]][[space]]$read())
        rownames(ld) <- rownames(ft_metadata)
        ld
      })
      names(loadings) <- varm_names
    } else {
      loadings <- list()
    }
    
    # mod/.../X
    modalities <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]

      X <- read_matrix(view[['X']])

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
      } else {
	       obsm <- list()
      }

      obsm
    })
    names(mod_obsm) <- assays

    # mod/.../varm
    mod_varm <- lapply(assays, function(mod) {
      view <- h5[['mod']][[mod]]

      var <- read_with_index(view[['var']])
      if (is("var", "data.frame"))
        rownames(var) <- paste(mod, rownames(var), sep="-")

      if ("varm" %in% names(view)) {
        varm <- lapply(names(view[["varm"]]), function(space) {
          emb <- t(view[["varm"]][[space]]$read())
          rownames(emb) <- rownames(var)
          emb
        })
        names(varm) <- names(view[["varm"]])
      } else {
        varm <- list()
      }

      varm
    })
    names(mod_varm) <- assays

    # Only common observations can be read
    obs_names <- Reduce(intersect, lapply(modalities, colnames))

    # Create a Seurat object
    srt <- CreateSeuratObject(modalities[[1]][,obs_names], assay=names(modalities)[1])
    for (modality in names(modalities)[2:length(modalities)]) {
      srt[[modality]] <- subset(modalities[[modality]], cells = obs_names)
    }

    # Add joint embeddings
    for (emb in names(embeddings)) {
      emb_name <- toupper(gsub('X_', '', emb))

      maybe_loadings <- matrix()
      if (emb %in% names(OBSM2VARM)) {
        varm_key = OBSM2VARM[[emb]]
        maybe_loadings <- loadings[[varm_key]]
      } 
      srt[[emb_name]] <- CreateDimReducObject(
        embeddings = embeddings[[emb]][obs_names,,drop=FALSE], 
        loadings = maybe_loadings,
        key = paste0(emb_name, "_"),
        assay = DefaultAssay(srt),  # this is not true but an existing assay must be provided
      )
    }

    # Add modality-specific embeddings
    for (mod in names(mod_obsm)) {
      mod_embeddings <- mod_obsm[[mod]]
      for (emb in names(mod_embeddings)) {
        emb_name <- paste(mod, toupper(gsub('X_', '', emb)), sep = "")

        maybe_loadings <- matrix()
        if (emb %in% names(OBSM2VARM)) {
          varm_key = OBSM2VARM[[emb]]
          maybe_loadings <- mod_varm[[mod]][[varm_key]]
        }

        srt[[emb_name]] <- CreateDimReducObject(
          embeddings = mod_embeddings[[emb]][obs_names,,drop=FALSE],
          loadings = maybe_loadings,
          key = paste0(emb_name, "_"), 
          assay = mod,
        )
      }
    }

    # Close the connection
    h5$close_all()

    srt
  } else {
    stop("Provide 'mae' or 'seurat' as the second argument.")
  }
}

