# Adenine Base Editor analysis
# Core Module

# Helpers

#' Convert A/G variant ID vector into GenomicRanges::GRanges
#'
#' @param variant_ids Character vector of VCF variant identifier. Example:
#' "chr2L:1_A/G"
#'
#' @return GenomicRanges::GRanges object of the variants. All widths are one
#' as this is record of single base mutations. The following mcols are defined:
#' \itemize{
#'   \item `type`, "A/G" (strand = "+") or "T/C" (strand = "-")
#'   to reflect variant type.
#'   \item `variant_id`, the original input as record
#' }
#' @details Variant type is extracted but NOT checked at all. You must make sure
#' that only "A/G" and "T/C" are present as input.
variantID2GRanges = function(variant_ids){
  # Converts a variant ID character vector into a GRanges object
  #   mcol `type` will be "A/G" or "T/C" to reflect variant type
  #   mcol `variant_id` will be the original input as record

  # Extract genomic information from IDs
  seqnames <- stringr::str_extract(variant_ids, "(chr.*?):", group = 1)
  pos <- stringr::str_extract(variant_ids, ":([0-9]+)_", group = 1)
  pos <- as.numeric(pos)
  type <- stringr::str_extract(variant_ids, "_(.*)$", group = 1)
  strand <- ifelse( type == "A/G", "+", "-" )

  GenomicRanges::GRanges(
    seqnames = seqnames, ranges = IRanges::IRanges(start = pos, end = pos),
    strand = strand, type = type, variant_id = variant_ids
  )
}

#' Format summary of a numeric vector
#'
#' @param nums A numeric vector
#'
#' @return A character(1) showing standard summary statistics (i.e.,
#' 0/25/50/mean/75/100 quantiles).
prettyFormatSummary = function(nums){
  # Returns pretty formatted summary from given nums numeric vector
  sprintf("%.2f", summary(nums)) |>
    paste(collapse = "/")
}

# Exported Functions

#' Get 'confident' variants for downstream analysis
#'
#' @param variantTableList A list of variant tables. Each variant table must have
#' $variant_id column.
#' @param print.stats Boolean, whether to print statistics.
#' @param min.occurence Integer, minimum number of occurences for a variant to
#' be considered as 'confident'.
#'
#' @return A character vector of 'confident' variant IDs.
#' @export
#'
#' @examples
#' vignette("rABE-analysis")
getConfidentVariantIds = function(
    variantTableList, print.stats = TRUE, min.occurence = 2 ){
  # Takes a list of variant tables and return a vector of variant IDs for downstream analysis
  #   Current implementation: ID occurs in >=2 replicates
  lapply( variantTableList, function(varTable) varTable[["variant_id"]] ) |>
    unlist() |>
    table() -> idOccurence

  if (print.stats){
    stats <- table(idOccurence)
    message(
      sprintf("#total unique variants = %d, out of which:\n",
              length(idOccurence)),
      sprintf("    %d (%.2f%%) occurred in %s samples\n",
              as.numeric(stats),
              as.numeric(stats) / length(idOccurence) * 100, names(stats))
    )
  }

  names(idOccurence)[ idOccurence >= min.occurence ]
}

#' Merge a list of variant tables into one table.
#'
#' @param variantTableList list of variant tables.
#' @param variantIDsToKeep which variant IDs to keep in the merged frame
#'
#' @return A single tibble::tibble of merged variant tables with following columns:
#' \itemize{
#'   \item variant_id, character vector
#'   \item variant_data list of tibble::tibble, each with columns $countREF,
#'   $countALT, $mutatin.perc.
#' }
#' @export
#'
#' @examples
#' vignette("rABE-analysis")
mergeVariants = function(variantTableList, variantIDsToKeep){
  # Takes a list of variant tables and return a single table with columns
  #   `variant_id`    character
  #   `variant_data`    tibble with columns
  #       `countREF`    int
  #       `countALT`    int
  #       `mutation.perc`    dbl

  variantTableList <- lapply(
    variantTableList,
    function(varTable) dplyr::filter(varTable,
                                     variant_id %in% variantIDsToKeep))
  varTable <- Reduce(rbind, variantTableList)
  varTable |>
    dplyr::nest_by(variant_id, .key = "variant_data", .keep = FALSE)
}

#' Set difference of two variant tables
#'
#' @param varTableX variant table X
#' @param varTableY variant table Y
#'
#' @return variant table of X \ Y (based on $variant_id)
#' @export
#'
#' @examples
#' vignette("rABE-analysis")
diffVariants = function(varTableX, varTableY){
  # Set difference based on `variant_id` matching X \ Y

  sel <- !(varTableX$variant_id %in% varTableY$variant_id)
  message(
    sprintf("%d out of %d (%.2f%%) variants are found and removed",
            sum(!sel), length(sel), sum(!sel)/length(sel)*100 )
  )

  varTableX[ !(varTableX$variant_id %in% varTableY$variant_id) , ]
}

#' Filter a single merged variant table
#'
#' @param mergedVarTable Single tibble::tibble of merged variants
#' @param perc.limits Lower and upper limits of median mutation percentage.
#' @param plot Boolean of whether to plot a
#' simple histogram of mutation percentages
#'
#' @return Filtered merged variant table
#' @export
#'
#' @examples
#' vignette("rABE-analysis")
filterVariants = function(mergedVarTable, perc.limits, plot = TRUE){
  # Takes a merged variant table (i.e., with `variant_id` and `variant_data` columns)
  #   filter variants based on the following properties:
  #     `perc.limits`: median value of `mutation.perc` is within the specified range
  # if `plot`, provides simple histogram of `mutation.perc`

  mergedVarTable |>
    dplyr::rowwise(variant_id) |>
    dplyr::mutate(
      mutation.perc.median = stats::median(variant_data$mutation.perc)) |>
    dplyr::ungroup() -> mergedVarTable

  mergedVarTable |>
    dplyr::filter(
      dplyr::between(mutation.perc.median, perc.limits[1], perc.limits[2]))  -> results

  if (plot){
    graphics::hist(log10(results$mutation.perc.median),
         xlab = "log10(Median mutation rate)", ylab = "Number of variants",
         main = "Distribution after filtering", xlim = c(-2, 2))

    graphics::hist(log10(mergedVarTable$mutation.perc.median),
         xlab = "log10(Median mutation rate)", ylab = "Number of variants",
         main = "Distribution before filtering", xlim = c(-2, 2))
  }

  message(
    sprintf("%d out of %d (%.2f%%) variants are filtered\n",
            nrow(mergedVarTable) - nrow(results), nrow(mergedVarTable),
            (nrow(mergedVarTable) - nrow(results)) / nrow(mergedVarTable) * 100),
    sprintf("    # below limit = %d\n",
            sum(mergedVarTable$mutation.perc.median < perc.limits[1])),
    sprintf("    # above limit = %d\n",
            sum(mergedVarTable$mutation.perc.median > perc.limits[2])),
    "Median mutation rate summaries:\n",
    sprintf("    below limit: %s\n",
            prettyFormatSummary(
              mergedVarTable$mutation.perc.median[mergedVarTable$mutation.perc.median < perc.limits[1]])),
    sprintf("    above limit: %s\n",
            prettyFormatSummary(
              mergedVarTable$mutation.perc.median[mergedVarTable$mutation.perc.median > perc.limits[2]])),
    sprintf("    in range: %s\n",
            prettyFormatSummary(results$mutation.perc.median))
  )

  results <- results |> dplyr::select(-mutation.perc.median)
  results
}

#' Assign variant IDs to a GRanges object based on overlap
#'
#' @param variant_ids Character vector of variant IDs
#' @param targetGRanges GenomicRanges::GRanges object
#' @param prune Boolean, if TRUE will only rturn rows that overlap to >=1 variants
#'
#' @return GenomicRanges::GRanges with an additional mcol $variant_ids.
#' $variant_ids is a IRanges::CharacterList
#' @export
#'
#' @details This function considers strandedness when finding overlap. "A/G"
#' variants will only be assigned to "+" stranded ranges, and vice versa.
#' However, you must make sure that only "A/G" and "T/C" variants are
#' included in the input.
#' @examples
#' vignette("rABE-analysis")
assignToGRanges = function(variant_ids, targetGRanges, prune = TRUE){
  # For each range specified in targetGRanges,
  #   assign variants that resides in the range
  #   range could be exons, UTRs, etc. Anything that is stranded is ok.
  # If `prune`, will only return rows that overlap to >=1 variants
  # TODO: Cannot execute efficiently in a apply scenario
  #   every execution takes a long overhead time so far!

  # Get overlaps
  variantRanges <- variantID2GRanges(variant_ids)
  ovl <- GenomicRanges::findOverlaps(variantRanges, targetGRanges)

  # Get (target index) -> IDs mapping
  #   Extract results of findOverlaps
  targetIdxs <- S4Vectors::to(ovl)
  variantIdxs <- S4Vectors::from(ovl)
  variantNames <- variant_ids[variantIdxs]
  #   Group by targetIdx, result is tibble with columns
  #     `targetIdx` - int
  #     `variant_idx` - list(character())
  idTable <- tibble::tibble(
    variant_id = variantNames, targetIdx = targetIdxs)
  idTable <- dplyr::nest_by(idTable, targetIdx,
                            .key = "variant_ids", .keep = FALSE)
  idTable$variant_ids <- lapply(idTable$variant_ids, function(df) df$variant_id)
  #   Complete the mapping
  tibble::tibble(
    targetIdx = seq(from = 1, to = length(targetGRanges), by = 1)) |>
    dplyr::left_join(idTable, by = "targetIdx") -> idTable

  # Construct
  GenomicRanges::mcols(targetGRanges)[["variant_ids"]] <-
    IRanges::CharacterList(idTable$variant_ids)

  # Prune if asked
  if (prune){
    counts <-
      vapply(GenomicRanges::mcols(targetGRanges)[["variant_ids"]],
             length, NA_integer_)
    targetGRanges <- targetGRanges[counts > 0]
  }

  targetGRanges
}

#' Convert GRanges with assigned variants into a tibble
#'
#' @param assignedGRanges an "assigned GRanges"
#' (i.e., GenomicRanges::GRanges with names and `variant_ids`).
#'
#' @return A tibble::tibble with columns: \itemize{
#'   \item name, character of `name(assignedGRanges)`
#'   \item variant_ids, list(character) of assignedGRanges$variant_ids
#' }
#' @export
#'
#' @examples
#' vignette("rABE-analysis")
getNameVarTable = function(assignedGRanges){
  # Take an "assigned GRanges" (i.e., GRanges with names and `variant_ids`)
  # Convert into a data frame with columns:
  #   name - character
  #   variant_ids - list(character)

  namesVec <- names(assignedGRanges)
  variant_ids <- as.list(assignedGRanges$variant_ids)

  tibble::tibble(
    name = namesVec, variant_ids = variant_ids
  ) |>
    dplyr::nest_by(name, .key = "variant_ids", .keep = FALSE) -> results
  results$variant_ids <- lapply(
    results$variant_ids,
    function(vids) unlist(vids[["variant_ids"]]))

  results
}
