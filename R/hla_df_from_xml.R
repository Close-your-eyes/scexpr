#' Prepare a data frame from a xml file containing hla allele information
#'
#' The xml file required for this function may be downloaded from the GitHub repository of the IMGT/HLA organization:
#' https://github.com/ANHIG/IMGTHLA/blob/Latest/xml/hla.xml.zip. The function then creates a dataframe from this xml
#' which can be used subsequently to match reads of HLA-loci from bulk or single cell RNA seq. By default the sequence
#' of the 2nd + 3rd exon are returned in an extra column as these will most likely be used to infer the HLA type in case of MHC-I.
#'
#' HLA nomenclature: http://hla.alleles.org/nomenclature/naming.html
#'
#' @param file_path path to the xml file, may be in zipped format
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param replace_none_pg replace "None" in p_group and g_group with allele_protein and allele_coding, respectively
#' @param ... additional argument to the lapply function; e.g. mc.cores may be passed when parallel::mclapply is chosen above
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' hla_df_from_xml(file_path = "path/hla.xml", lapply_fun = parallel::mclapply,
#' mc.cores = parallel::detectCores())
#' }
hla_df_from_xml <- function(file_path,
                            lapply_fun = lapply,
                            replace_none_pg = T,
                            ...) {

  if (!requireNamespace("xml2", quietly = T)) {
    utils::install.packages("xml2")
  }

  if (grepl("zip$", file_path)) {
    utils::unzip(file_path, exdir = tempdir())
    file_path <- file.path(tempdir(), "hla.xml")
  }

  lapply_fun <- match.fun(lapply_fun)
  children <- xml2::xml_children(xml2::read_xml(file_path))
  df <- dplyr::bind_rows(lapply_fun(children, read_child, ...))
  rownames(df) <- NULL

  df$allele <- unlist(lapply(df$allele, function(x) {
    while(stringr::str_count(x, ":") < 3) {
      x <- paste0(x, ":01")
    }
    return(x)
  }))

  df$prefix <- substr(df$allele, 1,3)
  df$gene <- stringr::str_sub(stringr::str_extract(df$allele, "-.{1,}\\*"), 2, -2)
  df$gene[which(is.na(df$gene))] <- substr(gsub("-", "", df$allele[which(is.na(df$gene))]), 4,4) # for MIC

  df$allele_coding <- gsub(":[[:digit:]]{2}[[:alpha:]]{0,}$", "", df$allele)
  df$allele_coding <- gsub("MIC", "", gsub("HLA-", "", df$allele_coding))
  df$allele_group <- sapply(strsplit(df$allele_coding, ":"), "[", 1)
  df$allele_protein <- paste0(df$allele_group, ":", sapply(strsplit(df$allele_coding, ":"), "[", 2))
  df$seq_length <- nchar(df$seq)

  df$seq_Exon2_3 <- substr(df$seq, as.numeric(sapply(strsplit(df$Exon2, "_"), "[", 1)), as.numeric(sapply(strsplit(df$Exon3, "_"), "[", 2)))
  df$seq_Exon2_3_length <- nchar(df$seq_Exon2_3)

  if (replace_none_pg) {
    df$g_group <- ifelse(df$g_group == "None", df$allele_coding, df$g_group)
    df$p_group <- ifelse(df$p_group == "None", df$allele_protein, df$p_group)
  }

  return(df)
}

read_child <- function(x) {
  cc <- xml2::xml_children(x)
  attrs <- xml2::xml_attrs(x)
  tryCatch({
    seq <- as.character(xml2::xml_contents(xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]]))[1])
    g_group <- xml2::xml_attrs(cc[[which(xml2::xml_name(cc) == "hla_g_group")]])
    p_group <- xml2::xml_attrs(cc[[which(xml2::xml_name(cc) == "hla_p_group")]])

    # features
    ff <- xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]])[which(xml2::xml_name(xml2::xml_children(cc[[which(xml2::xml_name(cc) == "sequence")]])) == "feature")]
    type <- sapply(xml2::xml_attrs(ff), function(x) x[["name"]])
    bounds <- sapply(ff, function(x) xml2::xml_attrs(xml2::xml_child(x, 1)))[which(type != "Translation")]
    bounds <- sapply(bounds, function(x) paste(x, collapse = "_"))
    type <- gsub(" ", "", (type[which(type != "Translation")]))
    df <- data.frame(allele = attrs[["name"]], seq = seq, g_group = g_group, p_group = p_group)
    df <- cbind(df, tidyr::pivot_wider(data.frame(type, bounds), names_from = type, values_from = bounds))
    return(df)
  }, error = function(e) {
    return(NULL)
  })
}


