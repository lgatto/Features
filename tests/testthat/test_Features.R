data(feat1)

test_that("empty QFeatures", {
    feat0 <- QFeatures()
    expect_true(validObject(feat0))
    expect_true(isEmpty(feat0))
    expect_null(show(feat0))
})

test_that("Manual QFeatures", {
    ## code from inst/scripts/test_data.R
    ## used to generate feat1
    psms <- matrix(1:20, ncol = 2)
    colnames(psms) <- paste0("S", 1:2)
    rowdata <- DataFrame(Sequence = c("SYGFNAAR", "SYGFNAAR", "SYGFNAAR", "ELGNDAYK",
                                      "ELGNDAYK", "ELGNDAYK", "IAEESNFPFIK",
                                      "IAEESNFPFIK", "IAEESNFPFIK", "IAEESNFPFIK"),
                         Protein = c("ProtA", "ProtA", "ProtA", "ProtA", "ProtA",
                                     "ProtA", "ProtB", "ProtB", "ProtB", "ProtB"),
                         Var = 1:10)
    rownames(rowdata) <- rownames(psms) <- paste0("PSM", 1:10)
    coldata <- DataFrame(Group = 1:2)
    rownames(coldata) <- colnames(psms)
    psms <- SummarizedExperiment(psms, rowData = rowdata)
    feat2 <- QFeatures(list(psms = psms), colData = coldata)
    expect_true(validObject(feat2))
    ## subsetting
    expect_null(show(feat2))
    expect_equal(psms, feat2[[1]])
    expect_equal(psms, feat2[["psms"]])
    expect_equal(feat2, feat2[1:10, 1:2, 1])
    ## compare to serialised data
    ## data(feat1)
    ## expect_equivalent(feat1, feat2)
})


test_that("addAssay", {
    data(feat1)
    assay2 <- feat1[[1]]
    feat1 <- addAssay(feat1, assay2, name = "psms2")
    expect_identical(names(feat1), c("psms", "psms2"))
    expect_identical(feat1[[1]], feat1[[2]])
})

test_that("[,QFeatures", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
    expect_true(expect_warning(validObject(feat1[, , "psms"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , "peptides"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 1]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 2]),
                               regexp = "'experiments' dropped; see 'metadata'"))
})


test_that("RowData", {
    data(feat2)
    rd <- rowData(feat2)
    expect_identical(names(rd), paste0("assay", 1:3))
    expect_identical(rd[[1]], rowData(feat2[[1]]))
    expect_identical(rd[[2]], rowData(feat2[[2]]))
    expect_identical(rd[[3]], rowData(feat2[[3]]))
})


test_that("RowData<-", {
    data(feat2)
    feat3 <- feat4 <- feat5 <- feat2
    value <- rowData(feat2)
    value[["assay1"]]$Prot <- letters[seq_len(nrow(value[["assay1"]]))]
    value[["assay1"]] <- value[["assay1"]][, "Prot",  drop = FALSE]
    value[["assay1"]]$foo <- rep("bar", nrow(value[["assay1"]]))
    rowData(feat3) <- value[-3]
    ## assay not in value are untouched
    expect_identical(feat3[[2]], feat2[[2]])
    ## replacing by untouched rowData leads to the same rowData
    expect_identical(feat3[[3]], feat2[[3]])
    ## rowvars in value and in rowData are replaced
    expect_identical(rowData(feat3[[1]])$Prot, letters[seq_len(nrow(feat3[[1]]))])
    ## rowvars in value but not in rowData are added
    expect_true("foo" %in% colnames(rowData(feat3[[1]])))
    expect_identical(rowData(feat3[[1]])$foo, rep("bar", nrow(feat3[[1]])))
    ## rowvars in rowData but not in value are untouched
    expect_identical(rowData(feat3[[1]])$x, rowData(feat2[[1]])$x)
    ## The value is not a DataFrameList
    value2 <- lapply(value, as.data.frame)
    rowData(feat4) <- value2
    expect_identical(feat3, feat4)
    ## invalide value leads to warning
    names(value) <- NULL
    expect_warning(rowData(feat5) <-  value, 
                   regexp = "Could not find a common assay")
    expect_identical(feat5, feat2)
})


test_that("rowDataNames", {
    rdn <- rowDataNames(feat1)
    expect_identical(length(feat1), length(rdn))
    expect_identical(names(feat1), names(rdn))
    for (i in seq_along(length(feat1)))
        expect_identical(rdn[[i]], names(rowData(feat1[[i]])))
})


test_that("selectRowData", {
    x <- c("Sequence", "Protein")
    ft <- selectRowData(feat1, x)
    expect_identical(length(ft), length(feat1))
    expect_identical(names(ft), names(feat1))
    expect_identical(rowDataNames(ft)[[1]], x)
    expect_error(selectRowData(feat1))
    expect_message(ft <- selectRowData(feat1,
                                       c("Sequence", "Protein",
                                         "Var", "location", "pval",
                                         "var_not_found")))
    expect_identical(feat1, ft)
})


test_that("rbindRowData", {
    data(feat2)
    ## Rbind rowData from 1 assay
    rd <- rbindRowData(feat2, 1)
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1])
    expect_identical(colnames(rd), c("assay", "rowname", 
                                     colnames(rowData(feat2)[[1]])))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1]))
    ## Get all common variable from all assays
    rd <- rbindRowData(feat2, seq_along(feat2))
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2))
    expect_identical(colnames(rd), c("assay", "rowname", "Prot", "x"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, ]))
    ## Warning no common variables 
    a <- feat2[[1]]
    rowData(a) <- DataFrame(foo = "bar")
    feat3 <- addAssay(feat2, a, name = "assay4")
    expect_warning(rd <- rbindRowData(feat3, seq_along(feat3)),
                   regexp = "No common columns")
    expect_identical(DataFrame(), rd)
})


test_that("renaming", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides", fun = colSums)
    feat1 <- aggregateFeatures(feat1, "peptides", fcol = "Protein",
                               name = "proteins", fun = colSums)
    expect_true(validObject(feat1))
    feat2 <- feat1
    names(feat2) <- LETTERS[1:3]
    expect_true(validObject(feat2))
    # Expect errors
    feat2@assayLinks[[1]]@name <- "foo"
    expect_error(validObject(feat2),
                 regexp = "@names not valid")
    feat2 <- feat1
    feat2@assayLinks[[1]]@from <- "bar"
    expect_error(validObject(feat2),
                 regexp = "@from not valid")
    expect_error(names(feat2) <- 1:3,
                 regexp = "must be a character")
    expect_error(names(feat2) <- letters[c(1,2,2)],
                 regexp = "is duplicated$")
})

## This does not seem to work when runing R CMD check...
# test_that("plotting", {
#     data(feat2)
#     feat2 <- joinAssays(feat2, i = 1:3)
#     feat2 <- aggregateFeatures(feat2, 4, "Prot", name = "proteins")
#     ## Plot QFeautres with interactive = FALSE
#     ## expect_doppelganger creates a snapshot and compares to a test 
#     ## snapshot in /_snaps/
#     set.seed(1234)
#     vdiffr::expect_doppelganger("qFeatures-plot", plot(feat2, interactive = FALSE))
# })

test_that("assays must have unique rownames", {
    hlpsms <- hlpsms[1:10, ]
    ft1 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms", fname = "Sequence")
    rownames(ft1[[1]][1:2]) <- rep("1", 2)
    expect_error(validObject(ft1))
})


test_that("longFormat", {
    colData(feat2)$X <- 1:12
    ## Select a single colvars and rowvars
    lt <- longFormat(feat2, rowvars = "Prot", colvars = "X")
    ## Check dimensions
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Check content
    expect_identical(lt$Prot,
                     unname(unlist(lapply(rowData(feat2), 
                                          ## Repeat 4x because 4 samples
                                          function(x) rep(x$Prot, 4)))))
    ## Select a single colvars and no rowvars (make sure that 
    ## the implementation does not break the MAE implementation)
    lt <- longFormat(feat2, colvars = "X")
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+1L)
    ## Select multiple rowvars
    lt <- longFormat(feat2, rowvars = c("Prot", "x"))
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Test error: rowvars is missing in rowData
    expect_error(longFormat(feat2, rowvars = "y"),
                 regexp = "not found")
})


