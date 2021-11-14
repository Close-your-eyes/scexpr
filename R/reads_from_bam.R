#' Get reads from a bam file
#'
#' to do
#'
#' Read scores: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
#' CellRanger tags: Cell barcode (CR), Cell barcode read quality (CY), Alignment score (AS), UMI (UR),
#' UMI read quality (UR), Query hit index (HI), Number of reported alignments for query (NH), Number of mismatches per pair (nM),
#' Region type (E = exonic, N = intronic, I = intergenic) (RE)
#'
#' @param file_path path to a position-sorted bam file; index file (.bai) has to be in the same directory
#' @param genomic_ranges GRanges object with n genomic ranges
#' @param add_flags a named list of columns to add; each element must have n elements (length of genomic_ranges)
#' @param add_tags tags to extract from bam file, passed to Rsamtools::ScanBamParam(); character(0) for nothing; missing tags do not seem to matter
#' @param read_scores calculate read scores from PhredQuality
#' @param revcomp_minus_strand calculate reverse complement of reads on minus strand
#'
#' @return a data frame with flagged (annotated) reads
#' @export
#'
#' @examples
reads_from_bam <- function(file_path,
                           genomic_ranges,
                           add_tags = c("CR", "CY", "AS", "UR", "UY", "HI", "NH", "nM", "RE"),
                           add_flags = NULL,
                           read_scores = T,
                           revcomp_minus_strand = T) {

  if (!missing(add_flags)) {
    if (any(length(genomic_ranges) != sapply(add_flags, length))) {
      stop(paste0("All add_flags need to have the length of genomic_ranges, which is ", length(genomic_ranges), "."))
    }
  }

  print ("Reading BAM file.")
  params <- Rsamtools::ScanBamParam(which = genomic_ranges, what = Rsamtools::scanBamWhat(), tag = add_tags) # ... #reverseComplement = FALSE --> all seqs refer to the (+)Strand, reads are (-)Strand are provided as rev.comp
  reads <- Rsamtools::scanBam(file_path, param = params)

  if (!is.null(add_flags)) {
    for (k in seq_along(reads)) {
      for (p in seq_along(add_flags)) {
        reads[[k]][[names(add_flags)[p]]] <- add_flags[[p]][k]
        names(add_flags[[p]]) <- names(reads)
      }
    }
  }

  # start of a read always refers to the (+)Strand, so for reads on the (-)Strand start is actually the end, (see IGV browser, read details)
  reads <- dplyr::bind_rows(lapply(names(reads), function(x) {
    temp <- data.frame(readName = reads[[x]][["qname"]],
                       start = reads[[x]][["pos"]],
                       length = reads[[x]][["qwidth"]],
                       strand = reads[[x]][["strand"]],
                       seq = reads[[x]][["seq"]],
                       qual = reads[[x]][["qual"]],
                       stringsAsFactors = FALSE)

    for (i in add_tags) {
      temp[,i] <- reads[[x]][["tag"]][[i]]
    }

    if (!is.null(add_flags)) {
      for (k in seq_along(reads)) {
        for (p in seq_along(add_flags)) {
          temp[,names(add_flags)[p]] <- add_flags[[p]][x]
        }
      }
    }
    return(temp)
  }))

  if (read_scores) {
    print("Calculating read score.")
    reads$minQual <- sapply(methods::as(Biostrings::PhredQuality(reads$qual), "IntegerList"), min)
    reads$meanQual <- sapply(methods::as(Biostrings::PhredQuality(reads$qual), "IntegerList"), mean)
    reads$n_belowQ30 <- sapply(methods::as(Biostrings::PhredQuality(reads$qual), "IntegerList"), function(x) sum(x < 30))
  }

  if (revcomp_minus_strand) {
    print ("Calculating reverse-complement of reads on minus strand.")
    reads[which(reads$strand == "-"),"seq"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reads[which(reads$strand == "-"),"seq"])))
  }

  reads$end <- reads$start + reads$length

  return(reads)
}
