# Variable list for dplyr operations

utils::globalVariables(c(
  # ABEtools
  c("targetIdx", "subjectHits", "variant_ids"),
  c("variant_id", "variant_data", "mutation.perc.median"),
  c("name", "REF", "ALT", "countREF", "countALT"),
  # intron.properties
  c("queryHits", "score", "sequence.i", "sequence.width", "start")
))
