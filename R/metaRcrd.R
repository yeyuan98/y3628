# S3 vctrs object for metadata upkeeping

#' Internal vctrs methods
#'
#' @import vctrs
#' @import methods
#' @keywords internal
#' @name y3628-vctrs

methods::setOldClass(c("y3628_metaRcrd", "vctrs_rcrd", "vctrs_vctr"))

# constructor
new_metaRcrd <-
  function(fields, meta.fields, class = character()){
    # Check that attribute fields are present
    if (!all(meta.fields %in% names(fields))){
      rlang::abort("Metadata fields must all present.")
    }
    # Check that no field is metaRcrd (no nesting)
    inherit_mR <-
      vapply(X=fields, FUN.VALUE=logical(1),
             FUN= \(o) inherits(o, what = "y3628_metaRcrd"))
    if (any(inherit_mR)){
      rlang::abort("metaRcrd does not allow nesting.")
    }
    #   generate data fields
    data.fields <- names(fields)
    data.fields <- data.fields[!(data.fields %in% meta.fields)]
    new_rcrd(
      fields, class = c(class, "y3628_metaRcrd"),
      meta.fields = meta.fields, data.fields = data.fields
    )
  }

# Apply a function to all numeric fields
.applyNumFields <- function(metaRcrd, FUNC){
  # Get all numeric fields
  all.fields <- fields(metaRcrd)
  # For each numeric field, apply function
  for (fn in all.fields){
    currField <- field(metaRcrd, fn)
    if (is.numeric(currField))
      field(metaRcrd, fn) <- FUNC(currField)
  }
  return(metaRcrd)
}

#' `metaRcrd` record-style vector
#'
#' This creates a vector where each item is a record of data and metadata fields.
#' Data fields are used for record equality/comparison operations.
#' Metadata fields are conceptually "attributes" attached to the data fields.
#'
#' # Note
#'
#' Set operations must use `vctrs` methods (e.g., [vctrs::vec_set_union()]).
#' Base set operations are not generic and hence invalid.
#'
#' @param fields A named list or data.frame where each row is a record. Names
#' of this list are the field names for the record vector.
#' @param meta.fields A character vector giving fields that should be considered
#'  as "metadata" fields.
#' @param ... For future extensions. Must be empty.
#' @return An S3 vector of class `y3628_metaRcrd`.
#' @export
#' @seealso [vctrs::new_rcrd()]
#' @examples
#' require(vctrs)
#' require(tibble)
#'
#' ## Representing metadata of experimental samples
#'
#' today <- Sys.Date()
#' dates <- c(today, today, today, today+1) # when
#' conds <- factor(c("L", "M", "H", "L"), levels = c("L", "M", "H")) # condition
#' who <- c("me", "me", "me", "Ahri") # personnel
#' eg1 <- metaRcrd(list(date=dates, condition=conds, personnel=who), "personnel")
#'
#' dates <- c(today) # when
#' conds <- factor("L", levels = c("L", "M", "H")) # condition
#' who <- "Ahri" # personnel
#' eg2 <- metaRcrd(list(date=dates, condition=conds, personnel=who), "personnel")
#'
#' ### concatenate
#' eg_full <- vec_c(eg1, eg2) # c() works just fine too
#' eg_full
#'
#' ### equality/set operation
#' eg1 == eg2
#' vec_set_difference(eg1, eg2) # use vec_set_* methods
#'
#' ### comparison/sort
#' eg1 <= eg2
#' sort(eg_full)
#'
#' ### use in a tibble
#'
#' proj <- c("fancy project 1", "ok project 2") # what project
#' rec <- list(eg1, eg2) # experiment record for each project
#' project_record <- tibble(project = proj, record = rec)
#' project_record
#'
metaRcrd <-
  function(fields, meta.fields, ...){
    rlang::check_dots_empty()
    # Call new_
    obj <- new_metaRcrd(fields, meta.fields, ...)
    # Cast all numeric fields to double type and return
    return(.applyNumFields(obj, as.double))
  }

# format
#   used by print
#' @export
format.y3628_metaRcrd <-
  function(x, ...){
    # Format all numeric fields
    x <- .applyNumFields(x, \(n) format.default(n, digits = 2))
    # Format object to form (DATA) ~ (METADATA)
    all.data <- vec_data(x)
    data <- all.data[attr(x, "data.fields")]
    meta <- all.data[attr(x, "meta.fields")]

    out_data <- do.call(paste, c(data, sep=","))
    out_meta <- do.call(paste, c(meta, sep=","))
    out <- paste("(", out_data, ")", " ~ ", "(", out_meta, ")", sep = "")

    out
  }

# abbr name
#' @export
vec_ptype_abbr.y3628_metaRcrd <- function(x, ...) {
  "mRcrd"
}

#' @export
obj_print_footer.y3628_metaRcrd <- function(x, ...) {
  cat("Data Fields:", attr(x, "data.fields"), "\n", sep = " ")
  cat("Meta Fields: ", attr(x, "meta.fields"), "\n", sep = " ")
}

# Type coercion

#   Between metaRcrd objects
#       Must possess the same data and metadata fields
#' @export
vec_ptype2.y3628_metaRcrd.y3628_metaRcrd <-
  function(x, y, ...){

    # Check point and attribute fields
    check_data <- setequal(attr(x, "data.fields"), attr(y, "data.fields"))
    check_meta <- setequal(attr(x, "meta.fields"), attr(y, "meta.fields"))
    if (!check_data | !check_meta){
      # Incompatible; error
      why <- paste0(
        "INCOMPATIBLE: ",
        ifelse(check_data, "", "DATA"),
        ifelse(check_meta, "", " METADATA")
      )
      stop_incompatible_type(x, y, details = why, ...)
    }

    vec_ptype(x)
  }
#' @export
vec_cast.y3628_metaRcrd.y3628_metaRcrd <- function(x, to, ...) x

# proxy - equality
#' @export
vec_proxy_equal.y3628_metaRcrd <-
  function(x, ...){
    vec_data(x)[attr(x, "data.fields")]
  }
