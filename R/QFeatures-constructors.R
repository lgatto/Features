##' @title QFeatures from tabular data
##'
##' @description
##'
##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `QFeatures` object.
##'
##' @param table File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`. It is
##'     advised not to encode characters as factors.
##'
##' @param ecol A `numeric` indicating the indices of the columns to
##'     be used as expression values. Can also be a `character`
##'     indicating the names of the columns. Caution must be taken if
##'     the column names are composed of special characters like `(`
##'     or `-` that will be converted to a `.` by the `read.csv`
##'     function. If `ecol` does not match, the error message will
##'     display the column names as seen by the `read.csv` function.
##'
##' @param fnames An optional `character(1)` or `numeric(1)`
##'     indicating the column to be used as feature names.
##'
##' @param ... Further arguments that can be passed on to `read.csv`
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @param name An `character(1)` to name assay in the `QFeatures`
##'     object. If not set, `features` is used.
##'
##' @return An instance of class [QFeatures] or [SummarizedExperiment].
##'
##' @author Laurent Gatto
##'
##' @describeIn readQFeatures See description.
##'
##' @importFrom utils read.csv
##' @importFrom methods new validObject
##' @import SummarizedExperiment
##'
##' @seealso The [QFeatures] class for an example on how to use
##'     `readQFeatures` and how to further manipulate the resulting data.
##'
##' @md
##' @aliases readSummarizedExperiment
##' @export
##'
##' @examples
##'
##' ## Load a data.frame with PSM-level data
##' data(hlpsms)
##'
##' ## Create the QFeatures object
##' fts2 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms")
##' fts2
readQFeatures <- function(table, ecol, fnames, ..., name = NULL)  {
    se <- readSummarizedExperiment(table, ecol, fnames, ...)
    if (anyDuplicated(rownames(se))) {
        message("Making assay rownames unique.")
        rownames(se) <- make.unique(rownames(se))
    }
    cd <- DataFrame(row.names = colnames(se))
    if (is.null(name))
        name <- "features"
    el <- structure(list(se), .Names = name[1])
    al <- AssayLinks(AssayLink(name = name[1]))
    ans <- MultiAssayExperiment(el, colData = cd)
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = al)
}
##' @describeIn readQFeatures Convert tabular data from a spreadsheet or a
##' `data.frame` into a `SummarizedExperiment` object.
##' @export
readSummarizedExperiment <- function(table, ecol, fnames, ...) {
    if (is.data.frame(table)) xx <- table
    else {
        args <- list(...)
        args$file <- table
        if ("rownames" %in% names(args)) {
            if (missing(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
        args$stringsAsFactors <- FALSE
        xx <- do.call(read.csv, args)
    }
    if (is.character(ecol) || is.factor(ecol)) {
        ecol0 <- ecol
        ecol <- match(ecol0, colnames(xx))
        if (any(is.na(ecol)))
            stop("Column identifiers ",
                 paste(ecol0[is.na(ecol)], collapse = ", "),
                 " not recognised among\n",
                 paste(colnames(xx), paste = ", "))
    } else if (is.logical(ecol)) {
        if (length(ecol) != length(xx)) stop("Length of 'ecol' and 'table' do not match.")
        ecol <- which(ecol)
    }
    assay <- as.matrix(xx[, ecol, drop = FALSE])
    fdata <- DataFrame(xx[, -ecol, drop = FALSE])

    if (!missing(fnames)) {
        fnames <- fnames[1]
        if (is.numeric(fnames))
            fnames <- colnames(xx)[fnames]
        if (is.na(match(fnames, colnames(xx))))
            stop(fnames, " not found among\n",
                 paste(colnames(xx), paste = ", "))
        rownames(fdata) <- rownames(assay) <- fdata[, fnames]
    } else {
        rownames(fdata) <- rownames(assay) <- seq_len(nrow(assay))
    }
    SummarizedExperiment(assay, rowData = fdata)
}



##' @export
##' 
##' @rdname QFeatures-class
##' 
##' @param assayLinks An optional [AssayLinks] object.
QFeatures <- function(..., assayLinks = NULL) {
    ans <- MultiAssayExperiment(...)
    if (isEmpty(ans)) assayLinks <- AssayLinks()
    else {
        if (is.null(assayLinks))
            assayLinks <- AssayLinks(names = names(ans))
    }
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = assayLinks)
}

