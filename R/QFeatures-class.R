##' @title Quantitative MS QFeatures
##'
##' @description
##'
##' Conceptually, a `QFeatures` object holds a set of *assays*, each
##' composed of a `matrix` (or `array`) containing quantitative data
##' and row annotations (meta-data).  The number and the names of the
##' columns (samples) must always be the same across the assays, but
##' the number and the names of the rows (features) can vary. The
##' assays are typically defined as `SummarizedExperiment` objects. In
##' addition, a `QFeatures` object also uses a single `DataFrame` to
##' annotate the samples (columns) represented in all the matrices.
##'
##' The `QFeatures` class extends the
##' [MultiAssayExperiment::MultiAssayExperiment] and inherits all
##' the functionality of the
##' [MultiAssayExperiment::MultiAssayExperiment] class.
##'
##' A typical use case for such `QFeatures` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values. The largest assay (the one with the
##' highest number of features, PSMs in the example above) is
##' considered the main assay.
##'
##' The recommended way to create `QFeatures` objects is the use the
##' `readQFeatures()` function, that creates an instance from tabular
##' data. The `QFeatures` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##' @section Constructors:
##'
##' - `QFeatures(..., assayLinks)` allows the manual construction of
##'   objects. It is the user's responsability to make sure these
##'   comply. The arguments in `...` are those documented in
##'   [MultiAssayExperiment::MultiAssayExperiment()]. For details
##'   about `assayLinks`, see [AssayLinks]. An example is shown below.
##'
##' - The [readQFeatures()] function constructs a `QFeatures` object
##'   from text-based spreadsheet or a `data.frame` used to generate
##'   an assay. See the function manual page for details and an
##'   example.
##'
##' @section Accessors:
##'
##' - The `QFeatures` class extends the
##'   [MultiAssayExperiment::MultiAssayExperiment] class and inherits
##'   all its accessors and replacement methods.
##'
##' - The `rowData` method returns a `DataFrameList` containing the
##'   `rowData` for each assay of the `QFeatures` object. On the other
##'   hand, `rowData` can be modified using `rowData(x) <- value`,
##'   where `value` is a list of tables that can be coerced to `DFrame`
##'   tables. The names of `value` point to the assays for
##'   which the `rowData` must be replaced. The column names of each
##'   table are used to replace the data in the existing `rowData`. If
##'   the column name does not exist, a new column is added to the
##'   `rowData`.
##'
##' - The `rbindRowData` functions returns a `DFrame` table that
##'   contains the row binded `rowData` tables from the selected
##'   assays. Only rowData variables that are common to all assays are
##'   kept.
##'
##' - The `rowDataNames` accessor returns a list with the `rowData`
##'   variable names.
##'
##' - The `longFormat` accessor takes a `QFeatures` object and returns
##'   it in a long format `DataFrame`. Each quantitative value is
##'   reported on a separate line. `colData` and `rowData` data can
##'   also be added. This function is an extension of the `longFormat`
##'   function in the [MultiAssayExperiment::MultiAssayExperiment].
##'
##' @section Adding assays:
##'
##' - The [aggregateFeatures()] function creates a new assay by
##'   aggregating features of an existing assay.
##'
##' - `addAssay(x, y, name, assayLinks)`: Adds a new assay (or
##'   list of assays) `y` to the `QFeatures` instance `x`. `name`
##'   is a `character(1)` naming the single assay (default is
##'   `"newAssay"`), and is ignored if `y` is a list of
##'   assays. `assayLinks` is an optional [AssayLinks].
##'
##' @section Subsetting:
##'
##' - QFeatures object can be subset using the `x[i, j, k, drop =
##'   TRUE]` paradigm. See the argument descriptions for details.
##'
##' - The [subsetByFeature()] function can be used to subset a
##'   `QFeatures` object using one or multiple feature names that will
##'   be matched across different assays, taking the aggregation
##'   relation between assays.
##'
##' - The `selectRowData(x, rowvars)` function can be used to
##'   select a limited number of `rowData` columns of interest named
##'   in `rowvars` in the `x` instance of class `QFeatures`. All other
##'   variables than `rowvars` will be dropped. In case an element in
##'   `rowvars` isn't found in any `rowData` variable, a message is
##'   printed.
##'
##' @param object An instance of class [QFeatures].
##' 
##' @param x An instance of class [QFeatures].
##' 
##' @param y A single assay or a *named* list of assays. 
##' 
##' @param i `character()`, `integer()`, `logical()` or `GRanges()`
##'     object for subsetting by rows.
##'
##' @param j `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by `colData` rows.
##'
##' @param k `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by assays
##'
##' @param drop logical (default `TRUE`) whether to drop empty assay
##'     elements in the `ExperimentList`.
##'
##' @param ... See `MultiAssayExperiment` for details. For `plot`, 
##'     further arguments passed to `igraph::plot.igraph`.
##'
##'
##' @return See individual method description for the return value.
##'
##' @seealso
##'
##' - The [readQFeatures()] constructor and the [aggregateFeatures()]
##'   function. The *QFeatures* vignette provides an extended example.
##'
##' - The [QFeatures-filtering] manual page demonstrates how to filter
##'   features based on their rowData.
##'
##' - The [missing-data] manual page to manage missing values in
##'   `QFeatures` objects.
##'
##' - The [QFeatures-processing] and [aggregateFeatures()] manual pages
##'   and *Processing* vignette describe common quantitative data
##'   processing methods using in quantitative proteomics.
##'
##' @import MultiAssayExperiment ProtGenerics
##'
##' @name QFeatures
##'
##' @aliases QFeatures QFeatures-class class:QFeatures addAssay dims,QFeatures-method show,QFeatures-method [,QFeatures,ANY,ANY,ANY-method [,QFeatures,character,ANY,ANY-method
##'
##' @aliases rowDataNames selectRowData
##'
##' @rdname QFeatures-class
##'
##' @exportClass QFeatures
##'
##' @author Laurent Gatto
##'
##' @examples
##' ## ------------------------
##' ## An empty QFeatures object
##' ## ------------------------
##'
##' QFeatures()
##'
##' ## -----------------------------------
##' ## Creating a QFeatures object manually
##' ## -----------------------------------
##'
##' ## two assays (matrices) with matching column names
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' sample_names <- paste0("S", 1:4)
##' colnames(m1) <- colnames(m2) <- sample_names
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[1:4]
##'
##' ## two corresponding feature metadata with appropriate row names
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = rownames(m1))
##' df2 <- DataFrame(row.names = rownames(m2))
##'
##' (se1 <- SummarizedExperiment(m1, df1))
##' (se2 <- SummarizedExperiment(m2, df2))
##'
##' ## Sample annotation (colData)
##' cd <- DataFrame(Var1 = rnorm(4),
##'                 Var2 = LETTERS[1:4],
##'                 row.names = sample_names)
##'
##' el <- list(assay1 = se1, assay2 = se2)
##' fts1 <- QFeatures(el, colData = cd)
##' fts1
##' fts1[[1]]
##' fts1[["assay1"]]
##'
##' ## Rename assay
##' names(fts1) <- c("se1", "se2")
##'
##' ## Add an assay
##' fts1 <- addAssay(fts1, se1[1:2, ], name = "se3")
##'
##' ## Get the assays feature metadata
##' rowData(fts1)
##'
##' ## Keep only the Fa variable
##' selectRowData(fts1, rowvars = "Fa")
##'
##' ## -----------------------------------
##' ## See ?readQFeatures to create a
##' ## QFeatures object from a data.frame
##' ## or spreadsheet.
##' ## -----------------------------------
##'
NULL


## ----------------------------------
## QFeatures Class ChangeLog
##
## Version 0.1:
##  - Contains MatchedAssayExperiment
## Version 0.2:
##  - Contains MultiAssayExperiment (see issue 46)
## Version 0.3:
##  - Rename to QFeatures (see issue 89)

setClass("QFeatures",
         contains = "MultiAssayExperiment",
         slots = c(version = "character",
                   assayLinks = "AssayLinks"),
         prototype = prototype(
             version = "0.3"))


##' @rdname QFeatures-class
##' 
##' @exportMethod show
setMethod("show", "QFeatures",
          function(object) {
              if (isEmpty(object)) {
                  cat(sprintf("A empty instance of class %s", class(object)), "\n")
                  return(NULL)
              }
              n <- length(object)
              cat(sprintf("An instance of class %s", class(object)), "containing", n, "assays:")
              el <- experiments(object)
              o_class <- class(el)
              elem_cl <- vapply(el, class, character(1L))
              o_len <- length(el)
              o_names <- names(el)
              featdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[1]
              }, FUN.VALUE = integer(1L))
              sampdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[2]
              }, FUN.VALUE = integer(1L))
              if (n <= 7) {
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len), o_names, elem_cl,
                              featdim, sampdim), "\n")
              } else {
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len)[1:3], o_names[1:3], elem_cl[1:3],
                              featdim[1:3], sampdim[1:3]), "\n")
                  cat(" ...")
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len)[(n-2):n], o_names[(n-2):n], elem_cl[(n-2):n],
                              featdim[(n-2):n], sampdim[(n-2):n]), "\n")
              }
          })


## Function that creates a plotly network graph from an igraph object
## @param graph An i graph object.
## @param coords A n vertices by 2 matrix with the coordinates of the 
##     nodes.
##' @importFrom igraph get.edgelist layout_as_tree plot.igraph add_edges
##' @importFrom grDevices rgb
##' 
.plotlyGraph <- function(graph, coords) {
    stopifnot(inherits(graph, "igraph"))
    ## Initialize plotly
    pl <- plotly::plot_ly()
    ## Add edges
    el <- get.edgelist(graph)
    edge_coords <- sapply(1:nrow(el), function(i) {
        edge_coord <- c(coords[el[i, 1], 1],
                        coords[el[i, 1], 2],
                        coords[el[i, 2], 1],
                        coords[el[i, 2], 2])
    })
    pl <- plotly::add_segments(pl, 
                               x = edge_coords[1, ],
                               y = edge_coords[2, ],
                               xend = edge_coords[3, ],
                               yend = edge_coords[4, ],
                               line = list(color = "grey", width = 0.3,
                                           showarrow = TRUE))
    ## Add nodes
    pl <- plotly::add_markers(pl,
                              x = coords[, 1], y = coords[, 2], 
                              marker = list(color = rgb(0.8, 0.8, 0.8),
                                            size = 40),
                              text = names(V(graph)), 
                              hoverinfo = "text")
    ## Add labels
    pl <- plotly::add_text(pl,
                           x = coords[, 1], y = coords[, 2], 
                           hoverinfo = "text",
                           text = names(V(graph)))
    
    ## Edit plot plot
    axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
    plotly::layout(pl,
                   xaxis = axis,
                   yaxis = axis,
                   showlegend = FALSE)
}

##' @rdname QFeatures-class
##' 
##' @param interactive A `logical(1)`. If `TRUE`, an interactive graph
##'     is generated using `plotly`. Else, a static plot using `igraph`
##'     is generated. We recommend interactive exploration when the
##'     `QFeatures` object contains more than 50 assays.
##' 
##' @importFrom igraph make_graph layout_as_tree plot.igraph add_edges V
##' @export
plot.QFeatures <- function (x, interactive = FALSE, ...) {
    ## Check arguments
    if (!interactive & length(x) > 50) 
        warning("The QFeatures object contains many assays. You may ",
                "want to consider creating an interactive plot (set ",
                "'interactive = TRUE')")
    ## Create the network graph
    graph <- make_graph(edges = character(0), 
                                isolates = names(x))
    ## Add the edges = links between assays
    roots <- c()
    for (child in names(x)) {
        parents <- assayLink(x, child)@from
        for(parent in parents) {
            if (!is.na(parent)) {
                graph <- add_edges(graph, c(parent, child))
            } else {
                roots <- c(roots, child)
            }
        }
    }
    ## Tree layout
    coords <- layout_as_tree(graph, root = roots)
    rownames(coords) <- names(V(graph))
    ## Add an interleaved offset on the y coord for better rendering 
    ## when many assays have to be drawn
    interleaved <- c(seq(1, 5, 2), seq(2, 6, 2))
    offset <- rep_len(interleaved * 0.05 * max(coords[, 2]), 
                       length.out = nrow(coords))
    coords[, 2] <- coords[, 2] + offset
    ## Perform plotting 
    if (!interactive) {
        plot.igraph(graph, layout = coords, ...)
        return(invisible(NULL))
    } else {
        return(.plotlyGraph(graph, coords))
    }
}


## Function that prunes a `Hits` object from an `AssayLink` object, 
## making sure that the `Hits` object is still valid with respect to 
## a parent assay and its corresponding assay (self). The validity is
## ensured by removing missing features. 
## 
## @param hits A `Hits` object
## @param parent A `SummarizedExperiment` object or any object that 
##     inherits from it. This is the assay that `hits` links from.
## @param self A `SummarizedExperiment` object or any object that 
##     inherits from it. This is the assay that `hits` links to.
## 
.pruneHits <- function(hits, parent, self) {
    ## Get the feature names in the parent and self assay
    featnParent <- rownames(parent)
    featnSelf <- rownames(self)
    ## Check which links are still in parent and self
    inParent <- mcols(hits)$names_from %in% featnParent
    inSelf <- mcols(hits)$names_to %in% featnSelf
    ## Remove lost feature links 
    hits[inParent & inSelf, ]
}

## Function that prunes an `AssayLink` of a `QFeatures` object, 
## making sure that the `AssayLink` object is still valid with respect
## to a given `QFeatures` object.
## 
## @param al An `AssayLink` object
## @param object A `QFeatures` object. `al` will be adapted so that it
##     becomes valid when contained in `object`.
## 
.pruneAssayLink <- function(al, object) {
    ## Identify lost assays that need to be pruned
    lost <- !al@from %in% names(object)
    ## Prune the links to assays in `@name`. When an assay is lost, 
    ## the links to that assay are removed. 
    if (all(lost)) { ## If all parent assays are lost, return an 
        ## empty AssayLink
        al <- AssayLink(name = al@name)
        return(al)
    } else if (any(lost)) { ## If some parents are lost (only in case 
        ## of multiple parents), remove the link to the lost parent(s) 
        ## (`@from`) and subset the corresponding `Hits` object(s)
        al@from <- al@from[!lost]
        al@hits <- al@hits[!lost]
    }
    ## Prune the links to features in `@hits`. Even when an assay is 
    ## not lost, some of its features might be lost and the 
    ## corresponding `Hits` must be adapted.
    if (inherits(al@hits, "List")) { ## If the AssayLink contains a 
        ## HitsList object, iterate through each `Hits` element.
        al@hits <- mendoapply(function(hits, parent) {
            .pruneHits(hits, parent, self = object[[al@name]])
        }, hits = al@hits, parent = experiments(object)[al@from])
    } else { ## If the AssayLink contains a single Hits object
        al@hits <- .pruneHits(hits = al@hits, 
                              parent = object[[al@from]], 
                              self = object[[al@name]])
    }
    al
}

## Function that prunes the `AssayLinks` of a `QFeatures` object, 
## making sure that the `AssayLinks` object is still valid after 
## `QFeatures` subsetting
## 
## @param object A `QFeatures` object
## 
.pruneAssayLinks <- function(object) {
    ## Subset the AssayLinks
    object@assayLinks <- object@assayLinks[names(object)]
    ## Removed lost links in each AssayLink object
    object@assayLinks <- endoapply(object@assayLinks, 
                                   .pruneAssayLink, object = object)
    object
}

##' @rdname QFeatures-class
##' 
##' @importFrom methods callNextMethod
##' 
##' @exportMethod [
setMethod("[", c("QFeatures", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              ## Subset the assays
              ans <- callNextMethod(x, i, j, ..., drop)
              
              ## Prune the AssayLinks so that the `QFeatures` object
              ## remains valid
              ans <- .pruneAssayLinks(ans)
              
              ## Check new object
              if (validObject(ans))
                  return(ans)
          })


##' @rdname QFeatures-class
##' @importFrom BiocGenerics dims
##' @exportMethod dims
setMethod("dims", "QFeatures",
          function(x) vapply(experiments(x), dim, integer(2)))


##' @rdname QFeatures-class
setMethod("[", c("QFeatures", "character", "ANY", "ANY"),
          function(x, i, j, k, ..., drop = TRUE) {
              if (missing(j)) j <- TRUE
              if (missing(k)) k <- TRUE
              subsetByFeature(x, i)[, j, k]
          })

##' @rdname QFeatures-class
##'
##' @param use.names A `logical(1)` indicating whether the rownames of
##'     each assay should be propagated to the corresponding `rowData`.
##'
setMethod("rowData", "QFeatures",
          function(x, use.names = TRUE, ...) {
              List(lapply(experiments(x), function(xx)
                  mcols(xx, use.names = use.names, ...)))
          })

##' @rdname QFeatures-class
##'
##' @export
setReplaceMethod("rowData", c("QFeatures", "DataFrameList"),
                 function(x, value) {
                     i <- intersect(names(value), names(x))
                     if (!length(i)) {
                         warning("Could not find a common assay between ",
                                 "'names(value)' and names(object)")
                         return(x)
                     }
                     el <- experiments(x)
                     for (ii in i)
                         rowData(el[[ii]])[, colnames(value[[ii]])] <-
                         value[[ii]]
                     BiocGenerics:::replaceSlots(x,
                                                 ExperimentList = el,
                                                 check = FALSE)
                 })

##' @rdname QFeatures-class
##'
##' @export
setReplaceMethod("rowData", c("QFeatures", "ANY"),
                 function(x, value) {
                     value <- endoapply(value, as, "DataFrame")
                     value <- as(value, "List")
                     rowData(x) <- value
                     x
                 })

##' @rdname QFeatures-class
##'
##' @export
rbindRowData <- function(object, i)  {
    ## Extract the rowData and column names from the desired assay(s)
    rdlist <- rowData(object)[i]
    rdNames <- rowDataNames(object)[i]
    ## Get the common variables between the selected rowData
    commonCols <- Reduce(intersect, rdNames)
    if (!length(commonCols)) {
        warning("No common columns between rowData tables were found.")
        return(DataFrame())
    }
    ## Add assay and rowname to the rowData
    rdlist <- lapply(names(rdlist),
                     function(x) cbind(assay = x,
                                       rowname = rownames(rdlist[[x]]),
                                       rdlist[[x]][, commonCols]))
    ## Row bind all tables in one DataFrame
    rdlist <- do.call(rbind, rdlist)
    rownames(rdlist) <- NULL
    rdlist
}



##' @rdname QFeatures-class
##'
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay.
##'
##' @export
selectRowData <- function(x, rowvars) {
    stopifnot(inherits(x, "QFeatures"))
    rowvars <- as.character(rowvars)
    allvars <- unique(unlist(rowDataNames(x)))
    missingvars <- setdiff(rowvars, allvars)
    if (length(missingvars))
        message(length(missingvars), " missing/mis-typed rowvars.")
    for (i in seq_len(length(x))) {
        rd <- rowData(x[[i]])
        rowData(x[[i]]) <- rd[, colnames(rd) %in% rowvars]
    }
    x
}


##' @rdname QFeatures-class
##'
##' @importFrom Biobase fData
##'
##' @export
rowDataNames <- function(x) {
    stopifnot(inherits(x, "MultiAssayExperiment"))
    CharacterList(lapply(experiments(x),
                         function(xx) {
                             if (inherits(xx, "SummarizedExperiment"))
                                 colnames(rowData(xx))
                             else if (inherits(xx, "eSet"))
                                 colnames(Biobase::fData(xx))
                             else NA_character_
                         }))
}


##' @rdname QFeatures-class
##'
##' @param value The values to use as a replacement. See the
##'     corresponding section in the documentation for more details.
##'
##' @exportMethod names<-
setReplaceMethod("names", c("QFeatures", "character"),
                 function(x, value) {
                     key_vals <- cbind(names(x), value)
                     x <-  callNextMethod(x, value)
                     names(x@assayLinks) <- value
                     for (i in seq_len(length(x))) {
                         al <- x@assayLinks[[i]]
                         al@name  <- unname(key_vals[key_vals[, 1] == al@name, 2])
                         if (!is.na(al@from))
                             al@from <- unname(key_vals[key_vals[, 1] == al@from, 2])
                         x@assayLinks[[i]] <- al
                     }
                     x
                 })


##' @rdname QFeatures-class
##'
##' @param colvars A `character()` that selects column(s) in the
##'     `colData`.
##' @param index The assay indicator for `SummarizedExperiment`
##'     objects. A vector input is supported in the case that the
##'     `SummarizedExperiment` object(s) has more than one assay
##'     (default `1L`)
##'
##' @importFrom MultiAssayExperiment longFormat
##'
##' @export
longFormat <- function(object,
                       colvars = NULL,
                       rowvars = NULL,
                       index = 1L) {
    if (!is.null(rowvars)) {
        rdNames <- rowDataNames(object)
        misNames <- sapply(rdNames,
                           function (x) any(!rowvars %in% x))
        ## Check that all required
        if (any(misNames))
            stop("Some 'rowvars' not found in assay(s): ",
                 paste0(names(misNames)[misNames], collapse = ", "))
        ## Get long format table with quantification values and colvars
        longDataFrame <-
            MultiAssayExperiment::longFormat(object, colvars, index)
        ## Get the required rowData
        rds <- lapply(rowData(object),
                      function(rd) rd[, rowvars, drop = FALSE])
        rds <- do.call(rbind, rds)
        ## Merge the rowData to the long table
        cbind(longDataFrame,
              rds[longDataFrame$rowname, , drop = FALSE])
    } else {
        ## If rowvars is null, return the MAE longFormat output
        MultiAssayExperiment::longFormat(object, colvars, index)
    }
}


##' @param name A `character(1)` naming the single assay (default is
##'     `"newAssay"`). Ignored if `y` is a list of assays.
##' @param assayLinks An optional [AssayLinks].
##'
##' @md
##'
##' @rdname QFeatures-class
##'
##' @export
addAssay <- function(x,
                     y,
                     name = "newAssay",
                     assayLinks = AssayLinks(names = name)) {
    stopifnot(inherits(x, "QFeatures"))
    el0 <- x@ExperimentList@listData
    if (is.list(y)) el1 <- y
    else el1 <- structure(list(y), .Names = name[1])
    el <- ExperimentList(c(el0, el1))
    smap <- MultiAssayExperiment:::.sampleMapFromData(colData(x), el)
    if (inherits(assayLinks, "AssayLink"))
        assayLinks <- AssayLinks(assayLinks)
    new("QFeatures",
        ExperimentList = el,
        colData = colData(x),
        sampleMap = smap,
        metadata = metadata(x),
        assayLinks = append(x@assayLinks,
                            assayLinks))
}
