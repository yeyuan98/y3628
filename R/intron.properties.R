# Functions for analyzing intron properties

# Internal helpers

#' Stranded version of GenomicRanges::shift
#'
#' @param GRanges.input GenomicRanges::GRanges
#' @param shift how many bp to shift each range. For '+' stranded range,
#' a positive shift will be towards 3' of the genomic coordinate. For '-',
#' a positive shift will be towards 5' of the genome.
#'
#' @return Shifted GenomicRanges::GRanges.
.strandedShift <- function(GRanges.input, shift){
  # The GenomicRanges::shift is NOT strand-aware.
  direction <- ifelse(
    GenomicRanges::strand(GRanges.input) == "+",
    1, -1
  )
  GenomicRanges::shift(
    GRanges.input, shift = direction * shift
  )
}

#' Count number of overlapping bases
#'
#' @param query GenomicRanges::Ganges of query.
#' @param subject GenomicRanges::Ganges of subject.
#'
#' @details
#' 'overlapping bases' is counted for each query against ALL subject ranges.
#' If a query overlaps with two subject ranges with 5 and 4 bases, the number
#' reported will be 5+4=9. Implementation is:
#' \enumerate{
#'   \item `GenomicRanges::subtract(query, subject)`
#'   \item Count up the width of the resulting GRangesList
#'   \item (original width of query) - (width of subtracted query)
#' }
#'
#' @return Integer vector of numbers of overlapping bases. Guaranteed to be
#' the same order as the query ranges.
.overlapWidths <-
  function(query, subject){

    subtracted <- GenomicRanges::subtract(query, subject)
    width.subtracted <- GenomicRanges::width(subtracted)
    width.subtracted <-
      vapply(X = width.subtracted, FUN.VALUE = NA_integer_,
             FUN = function(v) sum(v), USE.NAMES = FALSE)
    return(
      GenomicRanges::width(query) - width.subtracted
    )

  }

# Splicing site strength scoring by MaxEntScan
#   Cite: http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html
#' Splicing site strength scoring by MaxEntScan
#'
#' @param path.zip.MES Path to the MaxEntScan perl program (zipped). While not
#' provided by the package, you may obtain a copy from the [original MaxEntScan
#' author](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html).
#' Alternatively, you may get an archived copy from [Github](
#' https://github.com/yeyuan98/archived_external_resources).
#' @param BSgenome BSgenome object from Bioconductor.
#' @param GRange.intron GenomicRanges::GRanges of the introns. Genome must match
#' that of the BSgenome object. This intron ranges must be stranded (i.e.,
#' only contains '+' and '-' strand values.)
#'
#' @return data.frame with the following columns. Rows are guaranteed to match
#' order of the input introns. \itemize{
#'   \item MaxEnt.5ss, score for the 5' splicing site
#'   \item MaxEnt.3ss, score for the 3' splicing site
#' }
#' @export
#'
#' @examples
#' vignette("intron-properties")
MaxEntScan <-
  function(
    path.zip.MES = "burgelab.maxent.zip",
    BSgenome, GRange.intron
  ){

    # Check for perl
    tryCatch(
      expr = system2("perl", args = "--version", stdout = NULL, stderr = NULL),
      warning = function(w) stop("perl not found."),
      error = function(e) stop("perl not found.")
    )

    # Check for MaxEntScan program
    if (!file.exists(path.zip.MES)){
      stop("MaxEntScan program not found.")
    }

    # Introns must be at least 20bp to run the algorithm
    if (
      any(GenomicRanges::width(GRange.intron) < 20)
    ) stop("Introns must be at least 20bp to run MaxEntScan.")

    # Original working directory
    original.wd <- getwd()

    # Unzip the program
    ME.dir = "./maxent"
    utils::unzip(path.zip.MES, exdir = ME.dir)

    # Grab 5' sequences (3bp exon + 6bp intron)
    gr.5shifted <-
      .strandedShift(GRange.intron, shift = -3L)
    gr.5shifted <-
      GenomicRanges::resize(gr.5shifted, width = 9L, fix = "start")
    fiveSS.seq <-
      BSgenome::getSeq(BSgenome, names = gr.5shifted) |> as.character()

    # 5' sequences save to file and run program
    val.5ss <- .MaxEntScanRun("score5.pl", fiveSS.seq, ME.dir)

    # Grab 3' sequences (20bp intron + 3bp exon)
    gr.3shifted <-
      .strandedShift(GRange.intron, shift = 3L)
    gr.3shifted <-
      GenomicRanges::resize(gr.3shifted, width = 23L, fix = "end")
    threeSS.seq <-
      BSgenome::getSeq(BSgenome, names = gr.3shifted) |> as.character()

    # 3' sequences save to file and run program
    val.3ss <- .MaxEntScanRun("score3.pl", threeSS.seq, ME.dir)

    # Clean up the unzipped MaxEntScan
    unlink(ME.dir, recursive = TRUE, force = TRUE)
    # Return data.frame
    data.frame(
      MaxEnt.5ss = val.5ss,
      MaxEnt.3ss = val.3ss
    )
  }

#' System cell of MaxEntScan perl script
#'
#' @param ps.MES What MaxEntScan script to run.
#' @param sequences character vector of sequences.
#' @param ME.dir Path to the unzipped MaxEntScan directory.
#'
#' @return numeric vector of the MaxEntScan return values.
.MaxEntScanRun <- function(ps.MES, sequences, ME.dir){
  # Invokes MaxEntScan and return numeric vector of the values

  original.wd <- getwd()

  fp <- tempfile()
  writeLines(text = sequences, con = fp)
  setwd(ME.dir)
  results <- system2("perl", args = c(ps.MES, fp), stdout = TRUE)
  setwd(original.wd)
  return(stringr::str_extract(results, "[0-9.]+") |> as.numeric())
}



#' Conservation scoring by phastCons
#'
#' @param GRange.intron GenomicRanges::GRanges of the introns.
#' @param bw.phastCons.path Path to a bigwig file of phastCons scores. This
#' may be retrieved from the UCSC (e.g., [dm6](
#' https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phastCons124way/)). Look for
#' the `dm6.phastCons124way.bw` file.
#' @param bed.phastCons.path Path to a bed file of conserved regions annotated
#' by phastCons. This may be retrieved from the UCSC (e.g., [dm6](
#' https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/)). Look for the
#' `phastConsElements124way.txt.gz` file.
#'
#' @return data.frame. Number of rows is the same as number of ranges of
#' `GRange.intron`. Results have the following columns \describe{
#'   \item{mean}{Mean phastCons values over each range}
#'   \item{perc.in.element}{Percentage of bases in
#'   conserved phastCons elements for each range}
#' }
#' @export
#'
#' @details
#' This function can be used for compute scores for any GRanges of interest.
#'
#' @examples
#' vignette("intron-properties")
phastCons <-
  function(
    GRange.intron, bw.phastCons.path, bed.phastCons.path
  ){

    # Mean score
    # Only load the intron scores, essential for memory efficiency
    bw.scores <- rtracklayer::import(
      con = bw.phastCons.path, which = GRange.intron
    )
    #   find overlap. Note that bigwig ranges are of width 1.
    #     therefore, mean will be over individual bases.
    ovl <- GenomicRanges::findOverlaps(GRange.intron, bw.scores)
    scores <-
      as.data.frame(ovl) |>
      dplyr::group_by(queryHits) |>
      dplyr::summarize(
        mean = mean(bw.scores$score[subjectHits])
      )
    #   check that all ranges have scores (otherwise, fill NAs for no scores)
    if (nrow(scores) != length(GRange.intron)){
      all.indices <- data.frame(
        queryHits = 1:length(GRange.intron)
      )
      scores <- dplyr::left_join(
        x = all.indices, y = scores, by = "queryHits"
      )
      warning("No phastCons score for some ranges. Filled with NAs.")
    }

    # Percentage of bases
    #   import BED ranges
    #   first column of the original UCSC data file needs to be removed
    #   that is index column for other databases
    fp <- tempfile()
    readr::read_delim(
      file = bed.phastCons.path, delim = "\t", col_names = FALSE,
      show_col_types = FALSE, progress = FALSE
    ) |> dplyr::select(-1) |> readr::write_delim(
      file = fp, delim = "\t", col_names = FALSE, progress = FALSE
    )
    bed.elements <- rtracklayer::import(
      con = fp, format = "bed"
    )
    #   get number of overlap bases for each range
    scores$perc.in.element <- .overlapWidths(GRange.intron, bed.elements)
    #     normalize by full width
    scores$perc.in.element <-
      scores$perc.in.element / GenomicRanges::width(GRange.intron) * 100

    return(scores |> dplyr::select(-1))
  }

#' Distance from branchpoint to 3' splicing site
#'
#' @param GRanges.intron GenomicRanges::GRanges of the introns.
#' @param branchpoint.motif An universalmotif::universalmotif-class object
#' representing the branch point. Can be loaded using universalmotif::read_*
#' methods. For example, `universalmotif::read_meme()`.
#' @param BSgenome BSgenome object from Bioconductor.
#' @param logodds.threshold logodds threshold used by
#' `universalmotif::scan_sequences()`.
#'
#' @return A numeric vector in order of ranges of GRanges.intron. Ranges that
#' do not have identified branchpoint will take NA values.
#' @export
#'
#' @examples
#' vignette("intron-properties")
BranchPointScan <-
  function(
    GRanges.intron, branchpoint.motif,
    BSgenome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
    logodds.threshold = .5
  ){

    # Get sequence and search for branch point
    seqs <- BSgenome::getSeq(BSgenome, GRanges.intron)
    scan <- universalmotif::scan_sequences(
      motifs = branchpoint.motif, sequences = seqs,
      threshold = logodds.threshold, threshold.type = "logodds"
    )

    # Tidy up
    #   For each sequence, we only take the highest score match
    #   Tie-breaker: select the match most 3' of the sequence

    #   Report statistics
    message(sprintf(
      "Threshold = %d%%, found >=1 branchpoint motifs %d/%d (%.1f%%)",
      as.integer(logodds.threshold*100), length(unique(scan$sequence.i)),
      length(seqs), length(unique(scan$sequence.i))/length(seqs)*100
    ))

    #   Max score for each sequence
    max.scores <-
      data.frame(
        row.index = 1:nrow(scan),
        sequence.i = scan$sequence.i, score = scan$score
      ) |>
      dplyr::group_by(sequence.i) |>
      dplyr::summarize(
        score = max(score)
      )

    #   Filter to retain max score
    scan <- scan |>
      tibble::as_tibble() |>
      dplyr::inner_join(
        max.scores, by = c("sequence.i", "score")
      )

    #   Tie-breaker
    max.start <- scan |>
      dplyr::group_by(sequence.i) |>
      dplyr::summarize(
        start = max(start)
      )
    scan <- scan |>
      dplyr::inner_join(max.start, by = c("sequence.i", "start"))

    # Return final values
    seq.width <- Biostrings::width(seqs)
    scan <- scan |>
      dplyr::mutate(
        sequence.width = seq.width[sequence.i],
        dist2ss3 = sequence.width - stop + 1 # Dist will be from the CTA[A]T
      )
    scan.all <- data.frame(
      sequence.i = seq_along(seqs)
    )
    scan.all <- scan.all |>
      dplyr::left_join(
        scan, by = "sequence.i"
      )

    return(scan.all$dist2ss3)
  }
