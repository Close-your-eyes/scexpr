#' Core cellular machinery gene patterns and gene selection
#'
#' Returns curated regex patterns representing conserved cellular
#' machinery programs (e.g. ribosome, transcription, translation,
#' RNA processing, DNA replication, proteostasis, oxidative
#' phosphorylation, and cytoskeletal transport). Optionally extracts
#' matching genes from an object containing gene features.
#'
#' The function provides two levels of annotation:
#'
#' * `cellular_machinery_patterns`: compact high-level categories
#' * `core_cell_biology_patterns`: more granular biological processes
#'
#' When `obj` is supplied, gene features are extracted using
#' `scexpr:::get_gene_features()` and all genes matching any
#' machinery pattern are returned.
#'
#' @param obj Optional object containing gene features. If `NULL`,
#'   the pattern definitions are returned. Otherwise matching genes
#'   are extracted from the object's feature set.
#'
#' @returns
#' If `obj = NULL`, a named list containing:
#'
#' * `cellular_machinery_patterns`
#' * `core_cell_biology_patterns`
#'
#' If `obj` is supplied, a character vector of unique genes matching
#' one or more cellular machinery patterns.
#'
#' @details
#' Pattern categories include:
#'
#' * Ribosome (`RPL`, `RPS`)
#' * Translation initiation and elongation (`EIF`, `EEF`)
#' * Transcription (`POLR`, `TAF`, `GTF`, `MED`)
#' * RNA splicing and processing (`SRSF`, `HNRNP`, `SF3B`, `PRPF`)
#' * mRNA export (`NXF`, `THOC`)
#' * DNA replication (`MCM`, `RFC`, `RPA`, `ORC`)
#' * Chromatin regulation (`SMARC`, `CHD`, `HDAC`)
#' * Protein folding (`HSP`, `DNAJ`)
#' * Proteasome (`PSMA`, `PSMB`, `PSMC`, `PSMD`)
#' * Oxidative phosphorylation (`ATP5`, `NDUF`, `COX`, `UQCR`)
#' * Cytoskeleton and intracellular transport (`ACT`, `TUB`, `KIF`)
#' * Mitochondrial transcripts (`MT-`)
#'
#' These signatures are intended for exploratory analysis of
#' housekeeping and cellular-state programs in single-cell and
#' bulk transcriptomic datasets.
#'
#' @export
#'
#' @examples
#' # Return pattern definitions
#' patterns <- core_cellular_machinery_genes()
#'
#' # Extract all matching genes from a Seurat object
#' genes <- core_cellular_machinery_genes(seurat_obj)
#'
#' # Access a specific pattern set
#' patterns[[1]]$ribosome
core_cellular_machinery_genes <- function(obj = NULL) {

  core_cell_biology_patterns <- list(

    ribosome = c(
      "^RPL",
      "^RPS"
    ),

    translation_initiation = c(
      "^EIF1",
      "^EIF2",
      "^EIF3",
      "^EIF4",
      "^EIF5",
      "^EIF6"
    ),

    translation_elongation = c(
      "^EEF",
      "^ETF",
      "^GSPT"
    ),

    transcription = c(
      "^POLR",
      "^TAF",
      "^GTF2",
      "^MED"
    ),

    rna_splicing = c(
      "^SRSF",
      "^HNRNP",
      "^SF3B",
      "^PRPF",
      "^U2AF",
      "^DDX",
      "^SNRNP"
    ),

    mrna_export = c(
      "^NXF",
      "^THOC",
      "^ALYREF"
    ),

    dna_replication = c(
      "^MCM",
      "^RFC",
      "^RPA",
      "^POLA",
      "^POLD",
      "^POLE",
      "^ORC"
    ),

    chromatin = c(
      "^SMARC",
      "^CHD",
      "^HIST1H",
      "^HIST2H",
      "^HAT1",
      "^HDAC"
    ),

    protein_folding = c(
      "^HSP",
      "^DNAJ",
      "^BAG"
    ),

    proteasome = c(
      "^PSMA",
      "^PSMB",
      "^PSMC",
      "^PSMD",
      "^PSME"
    ),

    mitochondria_oxphos = c(
      "^ATP5",
      "^NDUF",
      "^COX",
      "^UQCR",
      "^SDH"
    ),

    cytoskeleton = c(
      "^ACT",
      "^TUBA",
      "^TUBB",
      "^KIF",
      "^DYNC",
      "^MYH"
    ),

    mitochondrial = c(
      "^MT-"
    )
  )

  cellular_machinery_patterns <- list(

    ribosome = c("^RPL","^RPS"),

    translation = c(
      "^EIF",
      "^EEF",
      "^ETF",
      "^GSPT"
    ),

    transcription = c(
      "^POLR",
      "^TAF",
      "^GTF",
      "^MED"
    ),

    rna_processing = c(
      "^SRSF",
      "^HNRNP",
      "^SF3B",
      "^PRPF",
      "^U2AF",
      "^SNRNP"
    ),

    dna_replication = c(
      "^MCM",
      "^RFC",
      "^RPA",
      "^ORC",
      "^POLA",
      "^POLD",
      "^POLE"
    ),

    proteostasis = c(
      "^HSP",
      "^DNAJ",
      "^BAG",
      "^PSMA",
      "^PSMB",
      "^PSMC",
      "^PSMD"
    ),

    oxidative_phosphorylation = c(
      "^ATP5",
      "^NDUF",
      "^COX",
      "^UQCR",
      "^SDH"
    ),

    cytoskeleton_transport = c(
      "^ACT",
      "^TUB",
      "^KIF",
      "^DYNC",
      "^MYH"
    )
  )

  if (is.null(obj)) {
    return(list(
      cellular_machinery_patterns = cellular_machinery_patterns,
      core_cell_biology_patterns = core_cell_biology_patterns
    ))
  } else {
    genes <- scexpr:::get_gene_features(obj)
    out1 <- grep_genes(genes, cellular_machinery_patterns)
    out2 <- grep_genes(genes, core_cell_biology_patterns)
    return(unique(c(out1, out2)))
  }

}

grep_genes <- function(genes, patterns) {
  patterns <- unique(unlist(patterns))
  unique(unlist(lapply(patterns, function(p) grep(p, genes, value = TRUE))))
}
