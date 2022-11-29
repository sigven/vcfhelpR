#' Write BED data frame records
#'
#' @description
#' Sort and write BED records to file
#'
#'
#' @param bed_records data frame with column names and data similar to
#' mandatory BED columns (CHROM, START, END etc)
#' @param output_dir output directory for VCF file
#' @param bed_fname_prefix VCF file name prefix
#' @param genome_build human genome build (grch37/grch38)
#' @param fsep separator symbol to use in filenames
#' @param chrom_col name of chromosome column
#' @param start_col name of start column
#' @param end_col name of end column
#' @param keep_uncompressed logical indicating if uncompressed BED file should
#' be kept
#'
#' @export
#'
#'
write_bed_records <- function(
    bed_records = NULL,
    output_dir = NA,
    bed_fname_prefix = NA,
    genome_build = "grch37",
    fsep = "_",
    chrom_col = "chrom",
    start_col = "start",
    end_col = "end",
    keep_uncompressed = F) {

  invisible(assertthat::assert_that(
    !is.null(bed_records),
    msg = "Argument 'bed_records' must be a non-NULL object"))
  invisible(assertthat::assert_that(
    is.data.frame(bed_records),
    msg = paste0("Argument 'bed_records' must be of type data.frame, not ",
                 class(bed_records))))
  assertable::assert_colnames(
    bed_records, c(chrom_col, start_col, end_col),
    only_colnames = F, quiet = T)
  if (nrow(bed_records) == 0) {
    return(bed_records)
  }

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(
      fmt = "%t [%L]  %m %j",
      timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(!is.na(bed_fname_prefix))
  stopifnot(!is.na(output_dir))
  stopifnot(dir.exists(output_dir))
  stopifnot(
    genome_build %in% c("grch37","grch38"))
  stopifnot(NROW(bed_records) > 0)

  bed_records <- order_bed_records(
    bed_records = bed_records,
    chrom_col = chrom_col,
    start_col = start_col,
    end_col = end_col)

  bed_fname <- list()
  ## Specify full path for BED records file and BED file
  bed_fname[['final']] <-
    file.path(
      output_dir,
      paste0(
        bed_fname_prefix,
        fsep,
        genome_build, ".bed"
      ))


  utils::write.table(bed_records,
                     file = bed_fname[['final']],
                     sep = "\t",
                     col.names = F,
                     quote = F,
                     row.names = F)

  system(paste0("bgzip -f -c ",
                bed_fname[['final']],
                " > ",
                bed_fname[['final']],
                ".gz"))

  system(paste0("tabix -p bed ",
                bed_fname[['final']],
                ".gz"))

  lgr::lgr$info(paste0(
    "Output_directory: ",
    output_dir
  ))
  lgr::lgr$info(paste0(
    "Output BED file: ",
    basename(bed_fname[['final']]), ".gz"))
  lgr::lgr$info(paste0(
    "Output index file: ",
    basename(bed_fname[['final']]), ".gz.tbi"))

  if (keep_uncompressed == F) {
    system(
      paste0("rm -f ", bed_fname[['final']]))
  }


}

#' Sort BED segments
#'
#' @param bed_records data frame with BED segment records
#' @param chrom_col name of chromosome column
#' @param start_col name of start column
#' @param end_col name of end column
#'
#' @export
#'
order_bed_records <- function(
    bed_records,
    chrom_col = "chrom",
    start_col = "start",
    end_col = "end") {

  invisible(assertthat::assert_that(
    !is.null(bed_records),
    msg = "Argument 'bed_records' must be a non-NULL object"))
  invisible(assertthat::assert_that(
    is.data.frame(bed_records),
    msg = paste0("Argument 'bed_records' must be of type data.frame, not ",
                 class(bed_records))))
  assertable::assert_colnames(
    bed_records, c(chrom_col, start_col, end_col),
    only_colnames = F, quiet = T)
  if (nrow(bed_records) == 0) {
    return(bed_records)
  }

  bed_records[, start_col] <-
    as.integer(bed_records[, start_col])
  bed_records[, end_col] <-
    as.integer(bed_records[, end_col])
  bed_records_sorted <- bed_records

  chr_prefix <- FALSE
  chromosome_names <- unique(bed_records[, chrom_col])
  for (m in chromosome_names) {
    if (startsWith(m, "chr")) {
      chr_prefix <- TRUE
    }
  }

  chr_order <- c(as.character(
    paste0("chr", c(1:22))), "chrX", "chrY")
  if (chr_prefix == FALSE) {
    chr_order <- c(as.character(c(1:22)), "X", "Y")
  }
  bed_records_sorted[, chrom_col] <-
    factor(bed_records_sorted[, chrom_col], levels = chr_order)
  bed_records_sorted <-
    bed_records_sorted[order(bed_records_sorted[, chrom_col]), ]

  bed_records_final <- NULL
  for (chrom in chr_order) {
    if (nrow(bed_records_sorted[!is.na(bed_records_sorted[, chrom_col]) &
                       bed_records_sorted[, chrom_col] == chrom, ]) > 0) {
      chrom_regions <- bed_records_sorted[bed_records_sorted[, chrom_col] == chrom, ]
      chrom_regions_sorted <-
        chrom_regions[with(chrom_regions,
                           order(chrom_regions[, start_col],
                                 chrom_regions[, end_col])), ]
      bed_records_final <-
        dplyr::bind_rows(
          bed_records_final, chrom_regions_sorted)
    }
  }
  return(bed_records_final)
}

#' Order VCF data frame records
#'
#' @description
#' Function that orders VCF records according to order
#' of chromosomes and chromosomal position
#'
#' @param vcf_records data frame with VCF records
#' @param chrom_var variable name in vcf_df for chromosome
#' @param pos_var variable name in vcf_df for chromosomal position
#'
#' @return vcf_df data frame with ordered VCF records
#'
#' @export
#'
order_vcf_records <- function(
    vcf_records = NULL,
    chrom_var = "CHROM",
    pos_var = "POS") {

  stopifnot(
    is.data.frame(vcf_records) &
      chrom_var %in% colnames(vcf_records) &
      pos_var %in% colnames(vcf_records))
  if (NROW(vcf_records) == 0)return(vcf_records)

  vcf_records |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) :=
        factor(
          !!rlang::sym(chrom_var),
          ordered = T,
          levels = c(as.character(seq(1:22)), "X", "Y"))) |>
    dplyr::arrange(
      !!rlang::sym(chrom_var),
      !!rlang::sym(pos_var)) |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) :=
        as.character(!!rlang::sym(chrom_var)))
}



#' Write VCF data frame records
#'
#' @description
#' Validate, sort and write VCF records to file
#'
#'
#' @param vcf_records data frame with column names and data similar to mandatory VCF columns (CHROM, POS, ID etc)
#' @param output_dir output directory for VCF file
#' @param header_lines character vector with VCF header lines
#' @param vcf_fname_prefix VCF file name prefix
#' @param genome_build human genome build (grch37/grch38)
#' @param fsep separator symbol to use in filenames
#' @param sample_names sample names for genotype columns
#' @param keep_uncompressed logical indicating if uncompressed output is kept
#' @param validate logical indicating if vcf-validator should be run on final VCF file
#' @export
#'

write_vcf_records <- function(
    vcf_records = NULL,
    output_dir = NULL,
    header_lines = NULL,
    vcf_fname_prefix = NULL,
    fsep = "_",
    genome_build = "grch38",
    sample_names = NULL,
    keep_uncompressed = F,
    validate = F) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(
      fmt = "%t [%L]  %m %j",
      timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(!is.null(header_lines))
  stopifnot(!is.null(vcf_records))
  stopifnot(!is.null(vcf_fname_prefix))
  stopifnot(is.data.frame(vcf_records))
  stopifnot(dir.exists(output_dir))
  stopifnot(
    genome_build %in% c("grch37","grch38"))
  stopifnot(NROW(vcf_records) > 0)

  samples_found <- F
  if(!is.null(sample_names)){
    stopifnot(is.character(sample_names))
    stopifnot(length(sample_names) >= 1)

    if(!("FORMAT" %in% colnames(vcf_records))){
      lgr::lgr$error(
        paste0("Missing FORMAT column in 'vcf_records' - ",
               "required when 'sample_names' is not NULL"))
      return(0)
    }

    samples_present <- unique(sample_names %in%
      colnames(vcf_records))

    samples_found <- T
    if(length(samples_present) > 1){
      samples_found <- F
    }else{
      if(samples_present == F){
        samples_found <- F
      }
    }
    if(samples_found == F){
      lgr::lgr$error(
        paste0("Could not find all sample names (",
               paste(sample_names, collapse=", "), ") ",
               "as columns in 'vcf_records' data frame"))
      return(0)
    }
  }

  options(scipen = 999)

  assertable::assert_colnames(
    vcf_records,
    c("CHROM", "POS", "ID",
      "REF", "ALT", "QUAL",
      "FILTER", "INFO"),
    quiet = T,
    only_colnames = F
  )

  vcf_records <- vcf_records |>
    dplyr::select(.data$CHROM,
                  .data$POS,
                  .data$ID,
                  .data$REF,
                  .data$ALT,
                  .data$QUAL,
                  .data$FILTER,
                  .data$INFO)

  if(samples_found == T){
    vcf_records <- vcf_records |>
      dplyr::select(.data$CHROM,
                    .data$POS,
                    .data$ID,
                    .data$REF,
                    .data$ALT,
                    .data$QUAL,
                    .data$FILTER,
                    .data$INFO,
                    .data$FORMAT,
                    dplyr::any_of(sample_names))
  }

  human_chromosomes <- data.frame(
    'CHROM' = c(as.character(seq(1:22)), "X", "Y"),
    stringsAsFactors = F)

  ## Check VCF records
  vcf_records <- vcf_records |>
    dplyr::mutate(
      REF = as.character(.data$REF)) |>
    dplyr::mutate(
      ALT = as.character(.data$ALT)) |>
    dplyr::mutate(REF = dplyr::if_else(
      .data$REF == "TRUE",
      as.character("T"),
      as.character(.data$REF)
    )) |>
    dplyr::mutate(ALT = dplyr::if_else(
      .data$ALT == "TRUE",
      as.character("T"),
      as.character(.data$ALT)
    )) |>
    dplyr::mutate(CHROM = as.character(
      stringr::str_replace(
        .data$CHROM, "^chr", ""
      ))) |>
    dplyr::mutate(POS = as.numeric(
      .data$POS
    )) |>
    dplyr::filter(
      !is.na(.data$CHROM) &
        !is.na(.data$ID) &
        !is.na(.data$POS) &
        !is.na(.data$REF) &
        !is.na(.data$ALT) &
        !is.na(.data$FILTER) &
        !is.na(.data$INFO)
    )

  if(samples_found == T){
    vcf_records <- vcf_records |>
      dplyr::filter(
        !is.na(.data$FORMAT)
      )
  }

  ## Check for the existence of duplicate VCF records
  duplicate_records <- vcf_records |>
    dplyr::group_by(
      .data$CHROM,
      .data$POS,
      .data$REF,
      .data$ALT) |>
    dplyr::summarise(
      num_records = dplyr::n(),
      .groups = "drop") |>
    dplyr::filter(
      .data$num_records > 1)

  if(NROW(duplicate_records) > 0){

    lgr::lgr$error(
      paste0("Found n = ", NROW(duplicate_records),
             " duplicate records in VCF file: ")
    )
    for(i in 1:NROW(duplicate_records)){
      lgr::lgr$error(
        paste0(
          paste(
            duplicate_records[i,]$CHROM,
            duplicate_records[i,]$POS,
            duplicate_records[i,]$REF,
            duplicate_records[i,]$ALT,
            sep = "\t"),
          " (n = ",duplicate_records[i,]$num_records,")"
        )
      )
    }
    return(0)
  }


  ## Check for the existence of records with
  ## non-DNA content in REF or ALT (should be excluded)
  non_dna_alleles <-
    vcf_records |>
    dplyr::filter(
      !stringr::str_detect(
        .data$REF,"^(A|T|C|G){1,}$") |
        !stringr::str_detect(
          .data$ALT,"^(A|T|C|G){1,}$")
    )

  if(NROW(non_dna_alleles) > 0){
    lgr::lgr$warn(
      paste0("Found n = ",
      NROW(non_dna_alleles),
      " VCF records with REF/ALT columns containing non-DNA bases")
    )

    vcf_records <- vcf_records |>
      dplyr::filter(
        stringr::str_detect(.data$REF,"^(A|C|T|G){1,}$") &
          stringr::str_detect(.data$ALT,"^(A|C|T|G){1,}$"))

    lgr::lgr$info(
      paste0("Excluding n = ", NROW(non_dna_alleles),
             " VCF records with non-DNA bases")
    )

  }

  ## Check for the existence of records with
  ## identical REF and ALT (should be excluded)
  duplicate_ref_alts <-
    vcf_records |>
    dplyr::filter(.data$REF == .data$ALT)

  if(NROW(duplicate_ref_alts) > 0){
    lgr::lgr$warn(
      paste0("Found n = ",
             NROW(duplicate_ref_alts),
             "VCF records with identical REF and ALT values")
    )
    vcf_records <- vcf_records |>
      dplyr::filter(
        .data$REF != .data$ALT)

    lgr::lgr$info(
      paste0("Excluded n = ", NROW(duplicate_ref_alts),
             " VCF records with identical REF and ALT values")
    )

  }

  vars_other_chrom <-
    vcf_records |>
    dplyr::anti_join(
      human_chromosomes, by = "CHROM")

  if(NROW(vars_other_chrom) > 0){
    lgr::lgr$warn(
      paste0("Found n = ", nrow(vars_other_chrom),
             " VCF records located on the mitochondrial chromosome ('M/MT'),",
             " unlocalized, or unplaced sequences")
    )

    vcf_records <- vcf_records |>
      dplyr::inner_join(human_chromosomes, by = "CHROM")

    lgr::lgr$info(
      paste0("Exluded n = ", NROW(vars_other_chrom),
             " VCF records located on non-ordinary chromosomes")
    )
  }

  stopifnot(NROW(vcf_records) > 0)

  vartype_stats <- list()
  vartype_stats[['n_snv']] <-
    vcf_records |>
    dplyr::filter(
      nchar(.data$REF) == 1 &
        nchar(.data$ALT) == 1) |>
    nrow()

  vartype_stats[['n_mnv']] <-
    vcf_records |>
    dplyr::filter(
      nchar(.data$REF) > 1 &
        nchar(.data$ALT) >1 &
        nchar(.data$REF) == nchar(.data$ALT)) |>
    nrow()

  vartype_stats[['n_del']] <-
    vcf_records |>
    dplyr::filter(
      nchar(.data$REF) > 1 &
        nchar(.data$ALT) < nchar(.data$REF)) |>
    nrow()

  vartype_stats[['n_ins']] <- vcf_records |>
    dplyr::filter(
      nchar(.data$ALT) > 1 &
        nchar(.data$ALT) > nchar(.data$REF)) |>
    nrow()


  lgr::lgr$info(
    paste0("Found n = ",
            NROW(vcf_records), " VCF records in input data frame")
  )
  lgr::lgr$info(
    paste0("Number of SVNs: ", vartype_stats[['n_snv']])
  )
  lgr::lgr$info(
    paste0("Number of MNVs: ", vartype_stats[['n_mnv']])
  )
  lgr::lgr$info(
    paste0("Number of deletions: ", vartype_stats[['n_del']])
  )
  lgr::lgr$info(
    paste0("Number of insertions: ", vartype_stats[['n_ins']])
  )

  ## Order VCFF records according to chrom and pos
  lgr::lgr$info(
    paste0("Ordering VCF records according to 'CHROM' and 'POS'")
  )
  vcf_records <- order_vcf_records(
    vcf_records = vcf_records)

  vcf_fname <- list()

  ## Specify full path for VCF records file and VCF file
  vcf_fname[['final']] <-
    file.path(output_dir, paste0(
      vcf_fname_prefix,
      fsep,
      genome_build, ".vcf"
    ))


  ## Check VCF header lines
  header_lines_filtered <- c()
  for(n in header_lines){
    if(!startsWith(n, "##")){
      lgr::lgr$warn(
        "Found VCF header lines not starting with '##' - ignoring"
      )
    }else{
      header_lines_filtered <-
        c(header_lines_filtered, n)
    }
  }

  record_header_line <-
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

  if(samples_found == T){
    record_header_line <-
      paste0(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
        paste(sample_names, collapse = "\t")
      )
  }

  header_lines_filtered <-
    c(header_lines_filtered, record_header_line)

  ## Write VCF header lines
  write(
    header_lines_filtered,
    file = vcf_fname[['final']],
    sep = "\n"
  )

  ## Write VCF records to records file
  utils::write.table(vcf_records,
              file = vcf_fname[['final']],
              sep = "\t",
              col.names = F,
              quote = F,
              append = T,
              row.names = F)

  system(paste0("bgzip -f -c ",
                vcf_fname[['final']],
                " > ",
                vcf_fname[['final']],
                ".gz"))

  system(paste0("tabix -p vcf ",
                vcf_fname[['final']],
                ".gz"))

  lgr::lgr$info(paste0(
    "Output_directory: ",
    output_dir
  ))
  lgr::lgr$info(paste0(
    "Output VCF file: ",
    basename(vcf_fname[['final']]), ".gz"))
  lgr::lgr$info(paste0(
    "Output index file: ",
    basename(vcf_fname[['final']]), ".gz.tbi"))

  if (keep_uncompressed == F){
    system(
      paste0("rm -f ", vcf_fname[['final']]))
  }


}

#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

utils::globalVariables(c("."))
