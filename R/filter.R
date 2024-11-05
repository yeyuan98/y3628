#' Filtering metaRcrd objects
#'
#' @param mRcrd An S3 vector of class `y3628_metaRcrd`.
#' @param what Filtering criteria. Quoted and evaluated to subset the record.
#'
#' @return Subsetted `y3628_metaRcrd`.
#' @export
#'
#' @examples
#' require(datasets)
#' my_mtcars <- datasets::mtcars[c("cyl", "hp", "mpg")]
#' my_mtcars <- metaRcrd(my_mtcars, meta.fields = "mpg")
#' mR_filter(my_mtcars, mpg > 20)
mR_filter <- function(mRcrd, what){
  what <- rlang::enquo(what)
  mRcrd_proxy <- vctrs::vec_proxy(mRcrd)

  which <- rlang::eval_tidy(what, data = mRcrd_proxy)

  return(mRcrd[which])
}

#' Applying summary function to string splits
#'
#' @param string Input vector, passed to [stringr::str_split()].
#' @param pattern Pattern to look for, passed to [stringr::str_split()].
#' @param summary Summary function that will be applied to each split.
#'
#' @return Vector same length as input vector
#' @export
#'
#' @examples
#' t <- c("0.1,0.3,0.7", "NA", "0.2", "")
#' str_split_summary(t, ",", \(x) mean(as.numeric(x)))
str_split_summary <- function(string, pattern, summary){

  # Split string
  splits <- stringr::str_split(string, pattern)

  # Apply summary function to get length=1 outputs (check correctness)
  splits <- lapply(splits, summary)
  lengths <- lapply(splits, length) |> unlist()
  if (any(lengths != 1))
    rlang::abort("Summary function must return length()=1.")

  # Unlist and return
  return(unlist(splits))
}
