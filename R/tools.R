# Convenience helper functions that are useful in general and belong somewhere.

#' Customizable summary of data frame variables
#'
#' @param .data a data frame (extension)
#' @param var.group variables to group, either single name or gvars()
#' @param .fns named list of summary functions
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> Selection of variables to summarize on
#'
#' @return a data frame (extension) of summary
#' @export
#'
#' @examples
#' # Single group variable
#' groupThenSummarize(mtcars, cyl, list(m=mean, s=sd), disp:wt)
#' # Multiple group variable
#' groupThenSummarize(mtcars, gvars(cyl,vs), list(m=mean, s=sd), disp:wt)
groupThenSummarize <- function(.data, var.group, .fns, ...){
  # S1. Quote grouping variable
  var.group <- rlang::enexpr(var.group)
  if (is.symbol(var.group)){
    # If symbol (single group variable)
    var.group <- rlang::exprs(!!var.group)
  } else{
    # Otherwise call
    if (rlang::call_name(var.group) != "gvars")
      rlang::abort("Use `gvars()` to specify multiple variables.")
    var.group <- rlang::call_args(var.group)
  }
  # S2. Perform group-then-summarize
  .data |>
    # Grouping with group_by
    #		var.group is a list, splice
    dplyr::group_by(!!!var.group) |>
    # Forwarding `...` selection to `across`
    dplyr::summarize(dplyr::across(
      c(...), .fns,
      .names = "{.fn}_{.col}"
    ), .groups = "drop")
}
