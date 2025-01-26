# Convenience helper functions that are useful in general and belong somewhere.

#' Grouping of data frame without taking up the ellipsis
#'
#' Motivation: the original `dplyr::group_by()` uses ellipsis to allow
#' flexibility in specifying grouping variables. However, sometimes a function
#' might need the ellipsis for other purposes. In those cases, it is desirable
#' to "save" the ellipsis by allowing the user to provide all grouping variables
#' in a single parameter.
#'
#' This function is not intended for end users.
#'
#' @param .data Data frame to perform `dplyr::group_by()`.
#' @param var.group **Quoted** group variable(s). Multiple variables must be
#' inside `gvars()`. To quote group variables use `rlang::enexpr()` in the
#' caller function.
#' @param ... Must be empty.
#'
#' @returns Grouped data frame.
#' @export
#'
#' @examples
#' # Correct use
#' f <- function(df, var.group){
#'     var.group <- rlang::enexpr(var.group)
#'     return(grouper(df, var.group))
#' }
#' # Single group variable
#' f(mtcars, cyl)
#' # Multiple group variable
#' f(mtcars, gvars(cyl,vs))
#'
#' # Wrong use
#' f <- function(df, var.group){ return(grouper(df, var.group)) }
#' # f(mtcars, cyl) will error

grouper <- function(.data, var.group, ...){
  # Check dots
  rlang::check_dots_empty0(...)

  # Check that var.group is already quoted
  if (!rlang::is_expression(var.group))
    rlang::abort("Please quote by `rlang::enexpr()` before calling grouper.")

  # Depending on symbol / call convert into spliceable expression list.
  if (is.symbol(var.group)){
    # If symbol (single group variable)
    var.group <- rlang::exprs(!!var.group)
  } else{
    # Otherwise call
    if (rlang::call_name(var.group) != "gvars")
      rlang::abort("Use `gvars()` to specify multiple variables.")
    var.group <- rlang::call_args(var.group)
  }

  # Perform `group_by()`
  .data <- .data |>
    # Grouping with group_by
    #		var.group is a list, splice
    dplyr::group_by(!!!var.group)

  return(.data)
}

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
#' #   cyl == 4 & vs == 0 group will have NA sd values.
#' #   This is because there is only one row in this group.
#' mtcars[mtcars$cyl == 4 & mtcars$vs == 0,]
groupThenSummarize <- function(.data, var.group, .fns, ...){
  # Quote and group
  var.group <- rlang::enexpr(var.group)
  .data <- grouper(.data, var.group)
  # Summarize by across()
  .data |>
    # Forwarding `...` selection to `across`
    dplyr::summarize(dplyr::across(
      c(...), .fns,
      .names = "{.fn}_{.col}"
    ), .groups = "drop")
}
