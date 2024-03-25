# Adenine Base Editor analysis
# Preprocessing Module

# Helpers

#' Get data frame of primary variants from a VCF file
#'
#' @param vcf VariantAnnotation::VCF object
#' @param variant_ids Identifiers to subset the VCF object
#'
#' @return A data.frame of counts for primary variants of each ID queried.
#'
#' @examples Not applicable. Internal use only.
getPrimaryVariantTable = function(vcf, variant_ids){
  # Subset a VCF to get variants with variant_ids
  #   then wrangle the variant data into a data frame with following columns
  #     variant_id    REF    ALT    countREF    countALT
  #     character character character int int
  # variant_ids must be a character vector
  vcf <- vcf[variant_ids]
  # Grab ALT types
  alts <- stringr::str_sub(variant_ids, start = -3L, end = -1L)
  alts <- stringr::str_split(alts, "/")
  REF <- vapply(X = alts, FUN.VALUE = NA_character_,
                FUN = function(alt) alt[1])
  ALT <- vapply(X = alts, FUN.VALUE = NA_character_,
                FUN = function(alt) alt[2])
  # Grab ALT counts
  counts <- VariantAnnotation::geno(vcf)$AD
  counts <- counts[variant_ids,]
  countREF <- vapply(X = counts, FUN.VALUE = NA_integer_,
                     FUN = function(count) count[1])
  countALT <- vapply(X = counts, FUN.VALUE = NA_integer_,
                     FUN = function(count) count[2])
  # Return table
  data.frame(
    variant_id = variant_ids, REF = REF, ALT = ALT,
    countREF = countREF, countALT = countALT
  )
}

#' Prune a GenomicRanges containing variants to retain only A->G edits
#'
#' @param rangeWithVariants GenomicRanges::GRanges with $variant_ids.
#' Each range can contain multiple variant_ids separated by `;`.
#' Example $variant_ids = "chr2L:1_A/G;chr2L:100_C/A"
#' @param stranded Boolean, whether to consider strand of the ranges
#'
#' @return GenomicRanges::GRanges with only A->G edits
#'
#' @examples Not applicable. Internal use only.
pruneAdenineEditVariants = function(rangeWithVariants, stranded = TRUE){
  # Prune range$variant_ids to include variant_ids that are A->G edits
  # For the effect of `stranded`, refer to analysis$filterAdenineEditVariants
  # TODO: Remove for loops to improve speed
  strands <- as.character(GenomicRanges::strand(rangeWithVariants))
  variants <- rangeWithVariants$variant_ids
  variants <- stringr::str_split(variants, ";")
  nBefore <- length(unique(unlist(variants))) # Num of variants before pruning
  if (stranded){
    for (idx in seq_along(strands)){
      currStrand <- strands[idx]
      currVariants <- variants[[idx]]
      stopifnot(currStrand %in% c("+", "-"))
      if (currStrand == "+"){
        currVariants <- currVariants[stringr::str_detect(currVariants, "A/G")]
      } else{
        currVariants <- currVariants[stringr::str_detect(currVariants, "T/C")]
      }
      variants[[idx]] <- paste(currVariants, collapse = ";")
    }
  } else{
    for (idx in seq_along(strands)){
      currVariants <- variants[[idx]]
      selected <- stringr::str_detect(currVariants, "A/G|T/C")
      currVariants <- currVariants[selected]
      variants[[idx]] <- paste(currVariants, collapse = ";")
    }
  }
  variants <- unlist(variants)
  nAfter <- length(unique(unlist(stringr::str_split(variants, ";"))))
  message(
    sprintf("%d out of %d (%.1f%%) unique variants are A->G edits",
            nAfter, nBefore, nAfter/nBefore*100)
  )
  rangeWithVariants$variant_ids <- variants
  rangeWithVariants
}

# Exported functions

#' Assign overlapping VCF variant IDs to a GRanges
#'
#' @param range GenomicRanges::GRanges object
#' @param vcf VariantAnnotation::VCF object
#' @param sep Separator to use if >1 variants overlap to a range
#'
#' @return Input `range` except adding a $variant_ids column
#' @export
#'
#' @examples Refer to `vignette("Adenine Base Editor analysis")`
assignVariantIds = function(range, vcf, sep = ";"){
  # For EACH range, assign variants overlapping to it
  #   Returns GRanges range with an extra colData column `variants` of VCF ids.

  vcfRanges <- SummarizedExperiment::rowRanges(vcf)
  ovl <- GenomicRanges::findOverlaps(vcfRanges, range)
  ovl <- as.data.frame(ovl)

  message(
    sprintf("%d out of %d (%.1f%%) variants overlap with >=1 ranges",
            length(unique(ovl$queryHits)), length(vcf),
            length(unique(ovl$queryHits)) / length(vcf) * 100)
  )
  message(
    sprintf("%d out of %d (%.1f%%) ranges overlap with >=1 variants",
            length(unique(ovl$subjectHits)), length(range),
            length(unique(ovl$subjectHits)) / length(range) * 100)
  )

  ovl$variant_ids <- names(vcfRanges)[ovl$queryHits]
  ovl |> dplyr::group_by(subjectHits) |>
    dplyr::summarize(variant.ids = paste(variant_ids, collapse = sep),
                     .groups = "drop") -> ovl
  range$variant_ids <-NA_character_
  range$variant_ids[ovl$subjectHits] <- ovl$variant.ids

  range

}

#' Filter variants to include only Adenine base editing
#'
#' @param assignedRange GenomicRanges::GRanges with $variant_ids column.
#' You can assign variants with the function `assignVariantIds()`.
#' @param stranded Boolean whether to consider strandedness. If TRUE,
#' A/G editing for '+' strand, T/C for '-' strand, raising error if any strand
#' takes '*' value.
#' @param pruning Whether to prune values in $variant_ids column. If TRUE,
#' $variant_ids will only retain adenine base editing events. If FALSE,
#' $variant_ids is not changed.
#'
#' @return GenomicRanges::GRanges containing ranges with at least one adenine
#' base editing event.
#' @export
#'
#' @examples Refer to `vignette("Adenine Base Editor analysis")`
filterAdenineEditVariants = function(assignedRange,
                                     stranded = TRUE, pruning = TRUE){
  # Input GRanges must have variant_ids column, annotated by assignVariantIds
  #   Returns a filtered GRanges which contains A->G variants
  # Two modes defined by stranded
  #   stranded == FALSE
  #     range is retained if variant_id holds A/G or T/C alts
  #     strand(assignedRange) is not used during filtering
  #   stranded == TRUE
  #     range is retained if either of the following holds
  #       A/G alts & strand(range) == "+"
  #       T/C alts & strand(range) == "-"
  #       Raise an error if strand(range) == "*"
  # If pruning == TRUE
  #   variant_ids will only retain those that are A->G variants
  #   else, variant_ids will not be modified
  nAll <- length(assignedRange)
  range <- assignedRange[!is.na(assignedRange$variant_ids),]
  nBefore <- length(range)
  if (stranded == FALSE){
    detect <- "A/G|T/C"
    range <- range[stringr::str_detect(range$variant_ids, detect)]
  } else{
    strands <- as.character(GenomicRanges::strand(range))
    vids <- stringr::str_split(range$variant_ids, ";")
    alts <- lapply(vids, function(v) stringr::str_sub(v, -3L, -1L))
    query <- lapply(seq_along(alts),
                    function(idx) paste0(strands[idx], alts[[idx]]))
    query <- vapply(query, FUN.VALUE = NA_character_,
                    FUN = function(s) paste(s, collapse = ";"))
    detect <- "\\+A/G|\\-T/C"
    range <- range[stringr::str_detect(query, detect)]
  }
  nAfter <- length(range)
  message(
    sprintf("Total = %d, w/ variant = %d (%.1f%%), w/ A-edit = %d (%.1f%%)",
            nAll, nBefore, nBefore/nAll*100, nAfter, nAfter/nAll*100)
  )
  if (pruning){
    range <- pruneAdenineEditVariants(range, stranded)
  }
  range
}

#' Return detailed variant information (mutation rate, counts for REF/ALT).
#'
#' @param vcf VariantAnnotation::VCF object
#' @param assignedRange GenomicRanges::GRanges with $variant_ids column.
#' You can assign variants with the function `assignVariantIds()`.
#'
#' @return A tibble::tibble
#' @export
#'
#' @examples
getVariantTable = function(vcf, assignedRange){
  # Take variant ID from the GRanges and return a variant info table with columns:
  #   `variant_id` = unique variants that are found
  #   `mutation.rate` = mutation rate in percentage, countALT / (countREF + countALT)
  #   `alt.count` = countALT
  #   `ref.count` = countREF

  variant_id <- GenomicRanges::mcols(assignedRange)[["variant_ids"]]
  variant_id <- stringr::str_split(variant_id, pattern = ";") |> unlist()
  variant_id <- unique(variant_id)

  getPrimaryVariantTable(vcf, variant_id) |>
    tibble::as_tibble() |>
    dplyr::select(-REF, -ALT) |>
    dplyr::mutate(mutation.perc = countALT / (countREF + countALT) * 100)
}
