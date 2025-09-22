#' Exports getSequence results to FASTA format
#'
#' @param sequences A data.frame that was the output of the [getSequence()]
#' function
#' @param file File to which you want to write the data
#'
#' @author Steffen Durinck
#' @author Hugo Gruson
#'
#' @keywords methods
#'
#' @examplesIf interactive()
#' mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#'
#' seq <- getSequence(
#'   id = "BRCA1",
#'   type = "hgnc_symbol",
#'   seqType = "cdna",
#'   mart = mart
#' )
#' exportFASTA(seq, file = "test.fasta")
#'
#' @export
exportFASTA <- function(sequences, file) {
  if (missing(sequences) || !is.data.frame(sequences)) {
    cli::cli_abort(
      c(
        "No data.frame given to write FASTA.",
        "i" = "The data.frame should be the output of the {.fn getSequence} function."
      )
    )
  }
  if (missing(file)) {
    cli::cli_abort("Please provide filename to write to.")
  }
  if (ncol(sequences) == 2) {
    writeLines(
      sprintf(">%s\n%s", sequences[, 2], sequences[, 1]),
      file,
      sep = "\n\n"
    )
  } else {
    writeLines(
      sprintf(
        ">chromosome_%s_start_%s_end_%s\n%s",
        sequences[, 1],
        sequences[, 2],
        sequences[, 3],
        sequences[, 4]
      ),
      file,
      sep = "\n\n"
    )
  }
}
