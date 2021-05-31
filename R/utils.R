.mudataversion <- "0.1.0"
.anndataversion <- "0.1.0"
.name <- paste0(getPackageName(), ".r")
.version <- as.character(packageVersion(getPackageName()))

#' @import hdf5r
open_h5 <- function(filename) {
    h5p_create <- H5P_FILE_CREATE$new()
    h5p_create$set_userblock(512)
    H5File$new(filename, mode="w", file_create_pl=h5p_create)
}

#' @import hdf5r
finalize_mudata <- function(h5) {
    h5$create_attr("encoding-type", "MuData", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .mudataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))

    filename <- h5$get_filename()
    h5$close_all()
    h5 <- file(filename, "r+b")
    writeChar(paste0("MuData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
    close(h5)
}

#' @import hdf5r
finalize_anndata_internal <- function(h5) {
    h5$create_attr("encoding-type", "AnnData", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .anndataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))
}

#' @import hdf5r
open_and_check_mudata <- function(filename) {
    if (readChar(filename, 6) != "MuData") {
        if (is_hdf5(filename)) {
            warning("The HDF5 file was not created by muon, we can't guarantee that everything will work correctly", call.=FALSE)
        } else (
            stop("The file is not an HDF5 file", call.=FALSE)
        )
    }

    H5File$new(filename, mode="r")
}
