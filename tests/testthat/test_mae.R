context("Creating .h5mu files from MultiAssayExperiment objects")
library(MultiAssayExperiment)
library(SingleCellExperiment)
library(muon)
library(hdf5r)
library(fs)  # for file_temp()

fileh5mu <- paste0(file_temp(), ".h5mu")

test_that("a model can be created from a simple MAE object", {
    # This is adapted from
    # https://www.bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html
    
    # colData
    patient.data <- data.frame(sex=c("M", "F", "M", "F"),
        age=38:41,
        row.names=c("Jack", "Jill", "Bob", "Barbara"))
    
    # assays
    arraydat <- matrix(seq(101, 108), ncol=4,
                       dimnames=list(c("ENST00000294241", "ENST00000355076"),
                                     c("array1", "array2", "array3", "array4")))
    coldat <- data.frame(slope53=rnorm(4),
                         row.names=c("array1", "array2", "array3", "array4"))
    exprdimred <- list("pca"=matrix(10:17, ncol=2))
    exprdat <- SingleCellExperiment(arraydat, colData=coldat, reducedDims=exprdimred)
    exprmap <- data.frame(primary=rownames(patient.data)[c(1, 2, 4, 3)],
                          assay=c("array1", "array2", "array3", "array4"),
                          stringsAsFactors = FALSE)

    methyldat <- matrix(1:8, ncol=4,
                        dimnames=list(c("ENST00000355077", "ENST00000383706"),
                                      c("methyl1", "methyl2", "methyl3",
                                        "methyl4")))
    methylmap <- data.frame(primary = c("Jack", "Jill", "Barbara", "Bob"),
                            assay = c("methyl1", "methyl2", "methyl3", "methyl4"),
                            stringsAsFactors = FALSE)

    microdat <- matrix(201:212, ncol=3,
                       dimnames=list(c("hsa-miR-21", "hsa-miR-191",
                                       "hsa-miR-148a", "hsa-miR148b"),
                                     c("micro1", "micro2", "micro3")))
    micromap <- data.frame(primary = c("Jack", "Barbara", "Bob"),
                           assay = c("micro1", "micro2", "micro3"),
                           stringsAsFactors = FALSE)

    nrows <- 5; ncols <- 4
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    rowRanges <- GRanges(rep(c("chr1", "chr2"), c(2, nrows - 2)),
                         IRanges(floor(runif(nrows, 1e5, 1e6)), width=100),
                         strand=sample(c("+", "-"), nrows, TRUE),
                         feature_id=sprintf("ID\\%03d", 1:nrows))
    names(rowRanges) <- letters[1:5]
    colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 2),
                         row.names= c("mysnparray1", "mysnparray2",
                                      "mysnparray3", "mysnparray4"))
    rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                rowRanges=rowRanges, colData=colData)
    rangemap <- data.frame(primary = c("Jack", "Jill", "Bob", "Barbara"),
                           assay = c("mysnparray1", "mysnparray2", "mysnparray3",
                                      "mysnparray4"), stringsAsFactors = FALSE)

    # sampleMap
    listmap <- list(exprmap, methylmap, micromap, rangemap)
    names(listmap) <- c("Affy", "Methyl 450k", "Mirna", "CNV gistic")
    dfmap <- listToMap(listmap)
    colnames(dfmap) <- c("assay", "primary", "colname")

    # objlist
    objlist <- list("Affy" = exprdat, "Methyl 450k" = methyldat,
                    "Mirna" = microdat, "CNV gistic" = rse)

    # MAE
    myMultiAssay <- MultiAssayExperiment(objlist, patient.data, dfmap)

    # Writing
    outfile <- fileh5mu
    result <- WriteH5MU(myMultiAssay, outfile)

    # Assert the data is saved
    expect_true(result)

    # Read back
    h5 <- H5File$new(outfile, mode="r")

    # Check all the assays are written
    assays <- h5[['mod']]$names
    assays_orig <- sort(names(assays(myMultiAssay)))
    expect_equal(assays, assays_orig)

    h5$close_all()
})


test_that("a MAE object can be created from an .h5mu file", {
    mae <- ReadH5MU(fileh5mu, "mae")
    expect_equal(names(mae)[1], "Affy")
    expect_equal(length(reducedDims(mae[["Affy"]])), 1)
})
