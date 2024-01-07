#' Get reads from a bam file
#'
#' This is basically a wrapper around Rsamtools::ScanBamParam and Rsamtools::scanBam. The output from scanBam is processed to a data frame and additional columns
#' are attached. Providing an exact range has been found to not always work as expected. E.g. there were reads in chr6 outside the exonic regions of HLA-A
#' that could be mapped to HLA-A. This may be an individual problem of the underlying BAM file (mapping). In order to not miss any relevant reads, one may pass
#' a wider genomic range for reads to return (e.g. whole chr6 if HLA loci are of interest, see example).
#'
#' Read scores: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
#' CellRanger tags: Cell barcode (CR), error-corrected Cell barcode (CB), Cell barcode read quality (CY), Alignment score (AS), UMI (UR),
#' UMI read quality (UY), Query hit index (HI), Number of reported alignments for query (NH), Number of mismatches per pair (nM),
#' Region type (E = exonic, N = intronic, I = intergenic) (RE)
#'
#' @param file_path path to a position-sorted bam file; index file (.bai) has to be in the same directory
#' @param genomic_ranges GRanges object with n genomic ranges
#' @param add_flags a named list of columns to add; each element must have n elements (length of genomic_ranges)
#' @param add_tags tags to extract from bam file, passed to Rsamtools::ScanBamParam(); character(0) for nothing; missing tags do not seem to matter
#' @param read_scores calculate read scores from PhredQuality
#' @param revcomp_minus_strand calculate reverse complement of reads on minus strand; passed as reverseComplement to Rsamtools::ScanBamParam;
#' for alignment to reference sequence, set it to FALSE; per default the sequence on + strand is returned (=FALSE), if one wants the
#' sequence from - strand then set to TRUE
#' @param lapply_fun lapply function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return a data frame of reads
#' @export
#'
#' @examples
#' \dontrun{
#' The seqnames of different BAM file may require different formats:
#' GenomicRanges::GRanges(seqnames = "chr1", strand = "+", ranges = IRanges::IRanges(start = start, end = end)) or
#' GenomicRanges::GRanges(seqnames = "1", strand = "+", ranges = IRanges::IRanges(start = start, end = end))
#'
#' # genomic range over part of chromosome 6 (or whole)
#' chr6 <- GenomicRanges::GRanges(seqnames = "6", strand = "+",
#' ranges = IRanges::IRanges(start = 29000000, end = 35000000))
#' # ranges = IRanges::IRanges(start = 1, end = 536870912))
#'
#' # alternatively multiple regions of HLA-A exons (in hg19)
#' # these may have to be obtained from the BAM file, e.g. IGV browser; or the reference genome
#' hlaa <- GenomicRanges::GRanges(
#'   seqnames = "6", strand = "+",
#'   ranges = IRanges::IRanges(
#'     start = c(29910247, 29910534, 29911045, 29911899, 29912277, 29912836, 29913011, 29913228),
#'     end = c(29910403, 29910803, 29911320, 29912174, 29912393, 29912868, 29913058, 29913661)
#'   )
#' )
#'
#' reads <- scexpr::reads_from_bam(
#'   file_path = "my_bam_path",
#'   genomic_ranges = chr6,
#'   lapply_fun = parallel::mclapply, mc.cores = parallel::detectCores()
#' )
#'
#' # passing multiple regions may return reads twice or multiple times
#' # if these reads overlap two or more of the regions (see ?scanBam and ?ScanBamParam)
#' reads <- scexpr::reads_from_bam(
#'   file_path = "my_bam_path",
#'   genomic_ranges = hlaa,
#'   lapply_fun = parallel::mclapply, mc.cores = parallel::detectCores(),
#'   add_flags = list(exon = c(1:8), gene = rep("A", 8))
#' )
#'
#' # filter and process reads
#' reads <- reads[which(reads$minQual >= 27), ]
#' reads <- reads[which(reads$n_belowQ30 <= 3), ]
#' # filter duplicate reads
#' # if additional flags like exons have been passed
#' # these columns will prevent dplyr::distinct from filtering
#' reads <- dplyr::distinct(reads, start, seq, .keep_all = T)
#' # only reads with standard nucleotides
#' reads <- reads[which(!grepl("[^ACTGU]", reads[,"seq",drop=T])),]
#' # readNames were found to be not unique in any case
#' # (same name for reads with different start and different seq)
#' reads$readName <- make.unique(reads$readName)
#' }
reads_from_bam <- function(file_path,
                           genomic_ranges,
                           add_tags = c("CR", "CB", "CY", "AS", "UR", "UY", "HI", "NH", "nM", "RE"),
                           add_flags = NULL,
                           read_scores = T,
                           revcomp_minus_strand = T,
                           lapply_fun = lapply,
                           ...) {

  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }
  if (!requireNamespace("Rsamtools", quietly = T)) {
    BiocManager::install("Rsamtools")
  }


  lapply_fun <- match.fun(lapply_fun)

  if (!missing(add_flags)) {
    if (any(length(genomic_ranges) != sapply(add_flags, length))) {
      stop(paste0("All add_flags need to have the length of genomic_ranges, which is ", length(genomic_ranges), "."))
    }
  }

  print("Reading BAM file.")
  params <- Rsamtools::ScanBamParam(which = genomic_ranges, what = Rsamtools::scanBamWhat(), tag = add_tags, reverseComplement = revcomp_minus_strand) # ... #reverseComplement = FALSE --> all seqs refer to the (+)Strand, reads are (-)Strand are provided as rev.comp
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
    temp <- data.frame(
      readName = reads[[x]][["qname"]],
      start = reads[[x]][["pos"]],
      length = reads[[x]][["qwidth"]],
      strand = reads[[x]][["strand"]],
      seq = reads[[x]][["seq"]], # as.character implicit
      qual = reads[[x]][["qual"]],
      mapq = reads[[x]][["mapq"]],
      cigar = reads[[x]][["cigar"]],
      stringsAsFactors = FALSE
    )

    for (i in add_tags) {
      temp[, i] <- reads[[x]][["tag"]][[i]]
    }

    if (!is.null(add_flags)) {
      for (k in seq_along(reads)) {
        for (p in seq_along(add_flags)) {
          temp[, names(add_flags)[p]] <- add_flags[[p]][x]
        }
      }
    }
    return(temp)
  }))

  if (read_scores) {
    print("Calculating read score.")
    xx <- methods::as(Biostrings::PhredQuality(reads$qual), "IntegerList")
    reads$readQualNum <- unlist(lapply_fun(xx, paste, collapse = ".", ...))
    reads$minQual <- unlist(lapply_fun(xx, min, ...))
    reads$meanQual <- unlist(lapply_fun(xx, mean, ...))
    reads$n_belowQ30 <- unlist(lapply_fun(xx, function(x) sum(x < 30), ...))
  }
  return(reads)
}
