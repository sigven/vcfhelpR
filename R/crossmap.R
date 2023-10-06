#' Crossmap VCF file to different assembly version
#'
#' @description
#' Function that "lifts" the VCF records to a different
#' assembly version of the human genome, using CrossMap.py
#'
#' @param target_genome_file fasta file of target genome (e.g. ucsc.hg19.fa)
#' @param direction hg38Tohg19 or hg19Tohg38
#' @param crossmap_cmd_path Full path to CrossMap.py
#' @param chain_file_directory Directory containing chain files
#' @param fsep separator to use in filenames (default '_')
#' @param source_vcf VCF file to be lifted
#' @param target_vcf VCF file with lifted coordinates
#' @param verbose logical indicating if function should log extensively
#'
#' @export
#'
crossmap_vcf <- function(
    target_genome_file = '/Users/sigven/research/DB/hg19/ucsc.hg19.fa',
    direction = "hg38Tohg19",
    crossmap_cmd_path = "/Users/sigven/miniconda3/envs/py36/bin/CrossMap.py",
    chain_file_directory = '/Users/sigven/research/DB/chainFiles',
    fsep = '_',
    source_vcf = NULL,
    target_vcf = NULL,
    verbose = T) {

  if (is.null(target_vcf)) {
    lgr::lgr$error("Please specifiy 'target_vcf'")
    return(0)
  }
  if (is.null(source_vcf)) {
    lgr::lgr$error("Please specifiy 'source_vcf'")
    return(0)
  }
  if (direction != 'hg38Tohg19' && direction != 'hg19Tohg38') {
    lgr::lgr$error(paste0('parameter \'direction\' (',
                         direction,') must be either hg38Tohg19 or hg19Tohg38'))
    return(0)
  }
  chainFile <- paste0(chain_file_directory,'/',direction,'.over.chain.gz')
  if (!file.exists(chainFile)) {
    lgr::lgr$error(
      paste0('Directory with crossmap chain files (parameter \'chain_file_directory\') (',
             chain_file_directory,') does not exist'))
    return(0)
  }
  if (!file.exists(target_genome_file)) {
    lgr::lgr$error(paste0('Target genome file (',
                         target_genome_file,') does not exist'))
    return(0)
  }

  if (!file.exists(crossmap_cmd_path)) {
    lgr::lgr$error(paste0('Path to CrossMap (',
                         crossmap_cmd_path,') does not exist'))
    return(0)
  }

  if (!file.exists(source_vcf)) {
    lgr::lgr$error(paste0('Source VCF file (',
                          source_vcf,') does not exist'))
  }

  if(direction == "hg19Tohg38"){
    target_genome_file = '/Users/sigven/research/DB/hg38/hg38.fa'
  }

  ## make temporary VCF file (with random index)
  rand_index <- floor(stats::runif(1, min = 0, max = 100000000))
  tmp_vcf_file <- paste0(
    paste('tmp','crossmapr',rand_index, sep = fsep),
    '.vcf')

  cmd <- paste0(crossmap_cmd_path,' vcf ',
                chainFile,' ',source_vcf,' ',
                target_genome_file,' ',
                tmp_vcf_file)
  if (verbose == T) {
    lgr::lgr$info(paste0("Command: ", cmd))
  }
  system(cmd)

  con <- file(tmp_vcf_file)
  vcf_header_lines <-
    readLines(con = con)
  close(con)

  vcf_header_lines <-
    vcf_header_lines[
      stringr::str_detect(vcf_header_lines,
                          "^##") &
        !stringr::str_detect(
          vcf_header_lines,
          "^##(liftOver|contig|originalFile|targetRefGenome)")]

  variants <- as.data.frame(
    readr::read_tsv(
      tmp_vcf_file,
      comment = "#", show_col_types = F,
      col_names = F)) |>

    ## ignore indels that are not properly left-aligned
    dplyr::filter(
      nchar(X4) == nchar(X5) |
        (nchar(X4) != nchar(X5) &
           substr(X4,1,1) == substr(X5,1,1)))


  variant_set <- variants |>
    dplyr::mutate(X1 = stringr::str_replace(X1,"chr","")) |>
    dplyr::filter(!stringr::str_detect(X1,"Un_|_random|_alt|_hap[0-9]")) |>
    dplyr::mutate(mut_id = paste(X1,X2,X4,X5, sep="_"))

  duplicated_variants <-
    plyr::count(variant_set$mut_id) |>
    dplyr::filter(freq > 1)

  variant_set_clean <- variant_set |>
    dplyr::anti_join(
      dplyr::select(duplicated_variants, x),
      by = c("mut_id" = "x")) |>
    dplyr::select(-mut_id) |>
    dplyr::rename(
      CHROM = X1, POS = X2, ID = X3,
      REF = X4, ALT = X5, QUAL = X6,
      FILTER = X7,
      INFO = X8)

  write_vcf_records(
    vcf_records = variant_set_clean,
    header_lines = vcf_header_lines,
    output_dir = dirname(target_vcf),
    genome_build_in_fname = F,
    vcf_fname_prefix = stringr::str_replace(
      basename(target_vcf),"\\.vcf","")
    )

  system(paste0('rm -f ',tmp_vcf_file,'*'))

  lgr::lgr$info(paste0("Crossmapped VCF file: ", target_vcf,".gz"))

}

#' Crossmap BED file to different assembly version
#'
#' @description
#' Function that "lifts" the BED records to a different
#' assembly version of the human genome, using CrossMap.py region
#'
#' @param direction hg38Tohg19 or hg19Tohg38
#' @param crossmap_cmd_path Full path to CrossMap.py file
#' @param chain_file_directory Directory containing chain files
#' @param wdir working directory
#' @param debug logical indicating if debug mode is on or off
#' @param fsep file name separator (default '_')
#' @param remap_ratio Minimum ratio of bases that must remap
#' (CrossMap.py region -r <remap_ratio>)
#' @param source_bed BED file to be lifted
#' @param target_bed BED file with lifted coordinates
#'
#' @export
#'
crossmap_bed <- function(
    direction = "hg38Tohg19",
    crossmap_cmd_path = "/Users/sigven/miniconda3/envs/py36/bin/CrossMap.py",
    chain_file_directory = '/Users/sigven/research/DB/chainFiles',
    wdir = NULL,
    debug = FALSE,
    fsep = '_',
    remap_ratio = 0.9,
    source_bed = NULL,
    target_bed = NULL) {

  if (is.null(target_bed)) {
    lgr::lgr$error('Please specifiy target_bed')
    return(0)
  }
  if (direction != 'hg38Tohg19' &&
     direction != 'hg19Tohg38' &&
     direction != 'hg18Tohg38') {
    lgr::lgr$error(
      paste0('parameter \'direction\' (',
             direction,
             ') must be either hg38Tohg19 or hg19Tohg38 or hg18Tohg38'))
    return(0)
  }
  chainFile <- paste0(chain_file_directory,'/',direction,'.over.chain.gz')
  if (!file.exists(chainFile)) {
    lgr::lgr$error(paste0(
      'Directory with crossmap chain files (parameter \'chain_file_directory\') (',
      chain_file_directory,') does not exist'))
    return(0)
  }


  if (!file.exists(source_bed)) {
    lgr::lgr$error(paste0('Source BED file (',source_bed,') does not exist'))
    return(0)
  }
  ## make temporary BED file (with random index)
  tmp_bed_file <- list()
  rand_index <- floor(stats::runif(3, min = 0, max = 100000000))
  for (i in 1:3) {
    tmp_bed_file[[as.character(i)]] <-
      paste0(paste('tmp','crossmapr',rand_index[i],
                   sep = fsep),'.bed')
    if (!is.null(wdir)) {
      if (dir.exists(wdir)) {
        tmp_bed_file[[as.character(i)]] <-
          file.path(
            wdir,
            paste0(paste('tmp','crossmapr',rand_index[i],
                         sep = fsep),'.bed')
          )
      }
    }
  }

  source_records <- utils::read.table(
    file = source_bed, sep = "\t", quote = "")

  ncol_bed <- ncol(source_records)

  if (direction == "hg19Tohg38") {
    source_records$V1 <- stringr::str_replace_all(
      paste0('chr', source_records$V1), 'chrchr', 'chr')
    utils::write.table(
      source_records,
      file = tmp_bed_file[['3']],
      quote = F,
      row.names = F,
      col.names = F)
  }else{
    source_records$V1 <- stringr::str_replace_all(
      source_records$V1, 'chr', '')
    utils::write.table(
      source_records,
      file = tmp_bed_file[['3']],
      quote = F,
      row.names = F,
      col.names = F)
  }

  system(paste0(crossmap_cmd_path,
                ' region -r ',
                remap_ratio, ' ',
                chainFile,' ',
                tmp_bed_file[['3']],
                ' ', tmp_bed_file[['1']]))

  #tmp_bed_file[['2']]
  system(paste0('egrep -v \'^##(lift|new)\' ',
                tmp_bed_file[['1']],
                ' | sed \'s/^chr//\' | ',
                'egrep -v \'Un_|_random|_alt|_hap[0-9]\' > ',
                tmp_bed_file[['2']]))
  system(paste0('cat ',tmp_bed_file[['2']],
                ' | egrep -v \"^#\" | egrep -v \"^[XYM]\" ',
                '| sort -k1,1n -k2,2n -k3,3n ',
                '| cut -f1-', ncol_bed, ' > ',
                target_bed))
  system(paste0('cat ',tmp_bed_file[['2']],
                ' | egrep -v \"^#\" | egrep \"^[XYM]\" ',
                '| sort -k1,1 -k2,2n -k3,3n ',
                '| cut -f1-', ncol_bed, ' >> ',
                target_bed))
  system(paste0('bgzip -f -c ',target_bed,' > ',target_bed,'.gz'))
  system(paste0('tabix -f -p bed ',target_bed,'.gz'))
  if (debug == F) {
    system(paste0('rm -f ',tmp_bed_file[['1']],'*'))
    system(paste0('rm -f ',tmp_bed_file[['2']],'*'))
    system(paste0('rm -f ',tmp_bed_file[['3']],'*'))
  }

  lgr::lgr$info(paste0("Crossmapped BED file: ", target_bed,".gz"))

}
