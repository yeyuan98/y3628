# Table reader that have consistent behavior with different file extensions
#   Motivation:
#     Different people prefer different file formats for tables (csv,tsv,xlsx)
#       even if the file should have the same meaning
#   Design:
#     Programmer specify column names, column types
#     The $reader will take in all common file types
#     The


#' Read 'table' file with flexible file extension
#'
#' @param col_names Column names, character vector
#' @param col_types Column types, same as `readr::read_csv()`. See details.
#' @param skip Number of lines to skip before reading data
#' @inheritParams rlang::args_dots_empty
#'
#' @return Reader function that accepts a file path
#' @export
#'
#' @details
#' Depending on the actual file format provided by the user, column types may
#' not be respected. For example, if the file is excel spreadsheet, then the
#' column types provided is mapped to that supported by the `readxl` package.
#'
#' @examples
#' # TODO
flexTableReader <- function(col_names, col_types, skip, ...){

  force(col_names)
  force(skip)
  col_types <- rlang::maybe_missing(col_types, default = NULL)

  rlang::check_dots_empty()

  ftr <- function(file){
    ext <- tools::file_ext(file)
    reader <- getTableReader(ext)
    params <- getNormalizedParams(ext, file = file,
      col_names = col_names, col_types = col_types, skip = skip)
    do.call(reader, params)
  }
  return(ftr)
}


#' Get filetype-specific table reader
#'
#' @param ext File extension.
#'
#' @return A suitable reader function.
#'
#' @examples
#' # Not exported.
getTableReader <- function(ext){
  switch(
    ext,
    csv = readr::read_csv,
    tsv = readr::read_tsv,
    xls = readxl::read_excel,
    xlsx = readxl::read_excel,
    rlang::abort("File extension not supported.")
  )
}


#' Get specific parameters for a table reader
#'
#' @param ext File extension.
#' @param ... Parameters, see details.
#'
#' @return Named list of normalized parameters.
#'
#' @details
#' Currently the following parameters are supported: col_names, col_types,
#' skip. For excel files, col_types is not supported yet.
#'
#' @examples
#' # Not exported
getNormalizedParams <- function(ext, ...){
  params <- rlang::list2(...)

  # Hard-coded conversions depending on extension type

  #   Excel
  conv_excel_type <- function(types){
    # TODO
    rlang::abort("Excel custom column types not implemented.")
  }
  conv_excel <- function(params){
    converted <- params
    if (!is.null(params[["col_types"]])){
      converted[["col_types"]] <- conv_excel_type(params[["col_types"]])
    }
    if (!is.null(params[["file"]])){
      converted[["path"]] <- params[["file"]]
      converted[["file"]] <- NULL
    }
    return(converted)
  }

  #   Overall convert switch
  convert <- switch(
    ext,
    csv = identity,
    tsv = identity,
    xls = conv_excel,
    xlsx = conv_excel,
    rlang::abort("File extension not supported.")
  )

  return(convert(params))
}
