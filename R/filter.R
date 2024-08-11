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
