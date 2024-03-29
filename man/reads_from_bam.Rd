% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reads_from_bam.R
\name{reads_from_bam}
\alias{reads_from_bam}
\title{Get reads from a bam file}
\usage{
reads_from_bam(
  file_path,
  genomic_ranges,
  add_tags = c("CR", "CB", "CY", "AS", "UR", "UY", "HI", "NH", "nM", "RE"),
  add_flags = NULL,
  read_scores = T,
  revcomp_minus_strand = T,
  lapply_fun = lapply,
  ...
)
}
\arguments{
\item{file_path}{path to a position-sorted bam file; index file (.bai) has to be in the same directory}

\item{genomic_ranges}{GRanges object with n genomic ranges}

\item{add_tags}{tags to extract from bam file, passed to Rsamtools::ScanBamParam(); character(0) for nothing; missing tags do not seem to matter}

\item{add_flags}{a named list of columns to add; each element must have n elements (length of genomic_ranges)}

\item{read_scores}{calculate read scores from PhredQuality}

\item{revcomp_minus_strand}{calculate reverse complement of reads on minus strand; passed as reverseComplement to Rsamtools::ScanBamParam}

\item{lapply_fun}{lapply function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested}

\item{...}{additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen}
}
\value{
a data frame of reads
}
\description{
This is basically a wrapper around Rsamtools::ScanBamParam and Rsamtools::scanBam. The output from scanBam is processed to a data frame and additional columns
are attached. Providing an exact range has been found to not always work as expected. E.g. there were reads in chr6 outside the exonic regions of HLA-A
that could be mapped to HLA-A. This may be an individual problem of the underlying BAM file (mapping). In order to not miss any relevant reads, one may pass
a wider genomic range for reads to return (e.g. whole chr6 if HLA loci are of interest, see example).
}
\details{
Read scores: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
CellRanger tags: Cell barcode (CR), error-corrected Cell barcode (CB), Cell barcode read quality (CY), Alignment score (AS), UMI (UR),
UMI read quality (UY), Query hit index (HI), Number of reported alignments for query (NH), Number of mismatches per pair (nM),
Region type (E = exonic, N = intronic, I = intergenic) (RE)
}
\examples{
\dontrun{
The seqnames of different BAM file may require different formats:
GenomicRanges::GRanges(seqnames = "chr1", strand = "+", ranges = IRanges::IRanges(start = start, end = end)) or
GenomicRanges::GRanges(seqnames = "1", strand = "+", ranges = IRanges::IRanges(start = start, end = end))

# genomic range over part of chromosome 6 (or whole)
chr6 <- GenomicRanges::GRanges(seqnames = "6", strand = "+",
ranges = IRanges::IRanges(start = 29000000, end = 35000000))
# ranges = IRanges::IRanges(start = 1, end = 536870912))

# alternatively multiple regions of HLA-A exons (in hg19)
# these may have to be obtained from the BAM file, e.g. IGV browser; or the reference genome
hlaa <- GenomicRanges::GRanges(
  seqnames = "6", strand = "+",
  ranges = IRanges::IRanges(
    start = c(29910247, 29910534, 29911045, 29911899, 29912277, 29912836, 29913011, 29913228),
    end = c(29910403, 29910803, 29911320, 29912174, 29912393, 29912868, 29913058, 29913661)
  )
)

reads <- scexpr::reads_from_bam(
  file_path = "my_bam_path",
  genomic_ranges = chr6,
  lapply_fun = parallel::mclapply, mc.cores = parallel::detectCores()
)

# passing multiple regions may return reads twice or multiple times
# if these reads overlap two or more of the regions (see ?scanBam and ?ScanBamParam)
reads <- scexpr::reads_from_bam(
  file_path = "my_bam_path",
  genomic_ranges = hlaa,
  lapply_fun = parallel::mclapply, mc.cores = parallel::detectCores(),
  add_flags = list(exon = c(1:8), gene = rep("A", 8))
)

# filter and process reads
reads <- reads[which(reads$minQual >= 27), ]
reads <- reads[which(reads$n_belowQ30 <= 3), ]
# filter duplicate reads
# if additional flags like exons have been passed
# these columns will prevent dplyr::distinct from filtering
reads <- dplyr::distinct(reads, start, seq, .keep_all = T)
# only reads with standard nucleotides
reads <- reads[which(!grepl("[^ACTGU]", reads[,"seq",drop=T])),]
# readNames were found to be not unique in any case
# (same name for reads with different start and different seq)
reads$readName <- make.unique(reads$readName)
}
}
