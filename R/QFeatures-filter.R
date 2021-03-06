##' @title Filter features based on their rowData
##'
##' @description
##'
##' The `filterFeatures` methods enables users to filter features
##' based on a variable in their `rowData`. The features matching the
##' filter will be returned as a new object of class `QFeatures`. The
##' filters can be provided as instances of class `AnnotationFilter`
##' (see below) or as formulas.
##'
##' @section Variable filters:
##'
##' The variable filters are filters as defined in the
##' [AnnotationFilter] package. In addition to the pre-defined filter,
##' users can arbitrarily set a field on which to operate. These
##' arbitrary filters operate either on a character variables (as
##' `CharacterVariableFilter` objects) or numerics (as
##' `NumericVariableFilters` objects), which can be created with the
##' `VariableFilter` constructor.
##'
##' @seealso The [QFeatures] man page for subsetting and the `QFeatures`
##'     vignette provides an extended example.
##'
##' @return An filtered `QFeature` object.
##'
##' @author Laurent Gatto
##'
##' @name QFeatures-filtering
##'
##' @rdname QFeatures-filtering
##'
##' @aliases filterFeatures filterFeatures,QFeatures,formula-method filterFeatures,QFeatures,AnnotationFilter-method CharacterVariableFilter NumericVariableFilter VariableFilter
##'
##' @examples
##'
##' ## ----------------------------------------
##' ## Creating character and numberic
##' ## variable filters
##' ## ----------------------------------------
##'
##' VariableFilter(field = "my_var",
##'                value = "value_to_keep",
##'                condition = "==")
##'
##' VariableFilter(field = "my_num_var",
##'                value = 0.05,
##'                condition = "<=")
##'
##' example(aggregateFeatures)
##'
##' ## ----------------------------------------------------------------
##' ## Filter all features that are associated to the Mitochondrion in
##' ## the location feature variable. This variable is present in all
##' ## assays.
##' ## ----------------------------------------------------------------
##'
##' ## using the forumla interface, exact mathc
##' filterFeatures(feat1, ~  location == "Mitochondrion")
##'
##' ## using the forumula intefrace, martial match
##' filterFeatures(feat1, ~startsWith(location, "Mito"))
##'
##' ## using a user-defined character filter
##' filterFeatures(feat1, VariableFilter("location", "Mitochondrion"))
##'
##' ## using a user-defined character filter with partial match
##' filterFeatures(feat1, VariableFilter("location", "Mito", "startsWith"))
##' filterFeatures(feat1, VariableFilter("location", "itochon", "contains"))
##'
##' ## ----------------------------------------------------------------
##' ## Filter all features that aren't marked as unknown (sub-cellular
##' ## location) in the feature variable
##' ## ----------------------------------------------------------------
##'
##' ## using a user-defined character filter
##' filterFeatures(feat1, VariableFilter("location", "unknown", condition = "!="))
##'
##' ## using the forumula interface
##' filterFeatures(feat1, ~ location != "unknown")
##'
##' ## ----------------------------------------------------------------
##' ## Filter features that have a p-values lower or equal to 0.03
##' ## ----------------------------------------------------------------
##'
##' ## using a user-defined numeric filter
##' filterFeatures(feat1, VariableFilter("pval", 0.03, "<="))
##'
##' ## using the formula interface
##' filterFeatures(feat1, ~ pval <= 0.03)
##'
##' ## you can also remove all p-values that are NA (if any)
##' filterFeatures(feat1, ~ !is.na(pval))
##'
##' ## ----------------------------------------------------------------
##' ## Negative control - filtering for an non-existing markers value
##' ## or a missing feature variable, returning empty results
##' ## ----------------------------------------------------------------
##'
##' filterFeatures(feat1, VariableFilter("location", "not"))
##'
##' filterFeatures(feat1, ~ location == "not")
##'
##' filterFeatures(feat1, VariableFilter("foo", "bar"))
##'
##' filterFeatures(feat1, ~ foo == "bar")
##'
##' ## ----------------------------------------------------------------
##' ## Example with missing values
##' ## ----------------------------------------------------------------
##'
##' data(feat1)
##' rowData(feat1[[1]])[1, "location"] <- NA
##' rowData(feat1[[1]])
##'
##' ## The row with the NA is not removed
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion")[[1]])
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = FALSE)[[1]])
##'
##' ## The row with the NA is removed
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = TRUE)[[1]])
##'
##' ## Note that in situations with missing values, it is possible to
##' ## use the `%in%` operator or filter missing values out
##' ## explicitly.
##'
##' rowData(filterFeatures(feat1, ~ location %in% "Mitochondrion")[[1]])
##' rowData(filterFeatures(feat1, ~ location %in% c(NA, "Mitochondrion"))[[1]])
##'
##' ## Explicit handling
##' filterFeatures(feat1, ~ !is.na(location) & location == "Mitochondrion")
##'
##' ## Using the pipe operator
##' library("magrittr")
##' feat1 %>%
##'    filterFeatures( ~ !is.na(location)) %>%
##'    filterFeatures( ~ location == "Mitochondrion")
NULL



##' @import AnnotationFilter
##' @exportClass CharacterVariableFilter
##' @rdname QFeatures-filtering
setClass("CharacterVariableFilter", contains = "CharacterFilter")

##' @exportClass NumericVariableFilter
##' @rdname QFeatures-filtering
setClass("NumericVariableFilter", contains = "DoubleFilter")


##' @param field `character(1)` refering to the name of the variable
##'     to apply the filter on.
##'
##' @param value `character()` or `integer()` value for the
##'     `CharacterVariableFilter` and `NumericVariableFilter` filters
##'     respectively.
##'
##' @param condition `character(1)` defining the condition to be used in
##'     the filter. For `NumericVariableFilter`, one of `"=="`,
##'     `"!="`, `">"`, `"<"`, `">="` or `"<="`. For
##'     `CharacterVariableFilter`, one of `"=="`, `"!="`,
##'     `"startsWith"`, `"endsWith"` or `"contains"`. Default
##'     condition is `"=="`.
##'
##' @param not `logical(1)` indicating whether the filtering should be negated
##'     or not. `TRUE` indicates is negated (!). `FALSE` indicates not negated.
##'     Default `not` is `FALSE`, so no negation.
##'
##' @export VariableFilter
##' @rdname QFeatures-filtering
VariableFilter <- function(field,
                           value,
                           condition = "==",
                           not = FALSE) {
    if (is.numeric(value))
        new("NumericVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition,
            not = not)
    else if (is.character(value))
        new("CharacterVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition,
            not = not)
    else
        stop("Value type undefined.")
}



##' @param object An instance of class [QFeatures].
##'
##' @param filter Either an instance of class [AnnotationFilter] or a
##'     formula.
##'
##' @param na.rm `logical(1)` indicating whether missing values should
##'     be removed. Default is `FALSE`.
##'
##' @param ... Additional parameters. Currently ignored.
##'
##' @exportMethod filterFeatures
##'
##' @rdname QFeatures-filtering
setMethod("filterFeatures",
          c("QFeatures", "AnnotationFilter"),
          function(object, filter, na.rm = FALSE, ...)
              filterFeaturesWithAnnotationFilter(object, filter, na.rm, ...))

##' @rdname QFeatures-filtering
setMethod("filterFeatures",
          c("QFeatures", "formula"),
          function(object, filter, na.rm = FALSE, ...)
              filterFeaturesWithFormula(object, filter, na.rm, ...))

##' @importFrom BiocGenerics do.call
filterFeaturesWithAnnotationFilter <- function(object, filter, na.rm, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      if (field(filter) %in% names(x))
                          do.call(condition(filter),
                                  list(x[, field(filter)],
                                       value(filter)))
                      else
                          rep(FALSE, nrow(x))
                  })
    sel <- lapply(sel, function(x) {
        x[is.na(x)] <- !na.rm
        x
    })
    if (not(filter)) sel <- lapply(sel, "!")
    object[sel, , ]
}


##' @importFrom lazyeval f_eval
filterFeaturesWithFormula <- function(object, filter, na.rm, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      tryCatch(lazyeval::f_eval(filter, data = as.list(x)),
                               error = function(e) rep(FALSE, nrow(x)))
                  })
    sel <- lapply(sel, function(x) {
        x[is.na(x)] <- !na.rm
        x
    })
    object[sel, , ]
}


## Internal function called by `filterFeaturesWithAnnotationFilter` when
## `condition` is `"contains"`
contains <- function(x, value) {
    ## Replace regex special character by regular character matching
    value <- gsub('([[:punct:]])', '\\[\\1\\]', value)
    ## Return whether elements in x contain value or not
    grepl(pattern = value, x = x)
}
