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

  ## make temporary VCF file (with random index)
  tmp_vcf_file <- list()
  rand_index <- floor(stats::runif(3, min = 0, max = 1000000))
  for (i in 1:3) {
    tmp_vcf_file[[as.character(i)]] <-
      paste0(paste('tmp','crossmapr',rand_index[i],
                   sep = fsep),'.vcf')
  }

  cmd <- paste0(crossmap_cmd_path,' vcf ',
                chainFile,' ',source_vcf,' ',
                target_genome_file,' ',
                tmp_vcf_file[['1']])
  if (verbose == T) {
    lgr::lgr$info(paste0("Command: ", cmd))
  }
  system(cmd)

  system(paste0('egrep -v \'^##(lift|new|contig)\' ',tmp_vcf_file[['1']],
                ' | sed \'s/^chr//\' | egrep -v \'Un_|_random|_alt|_hap[0-9]\' > ',
                tmp_vcf_file[['2']]))
  system(paste0('egrep \'^#\' ',tmp_vcf_file[['2']],' > ',target_vcf))
  system(paste0('cat ',tmp_vcf_file[['2']],
                ' | egrep -v \"^#\" | egrep -v \"^[XYM]\" ',
                '| sort -k1,1n -k2,2n -k4,4 -k5,5 >> ',
                tmp_vcf_file[['3']]))
  system(paste0('cat ',tmp_vcf_file[['2']],
                ' | egrep -v \"^#\" | egrep \"^[XYM]\" ',
                '| sort -k1,1 -k2,2n -k4,4 -k5,5 >> ',
                tmp_vcf_file[['3']]))
  system(paste0('cat ',tmp_vcf_file[['3']], ' >> ',
                target_vcf))
  # system(paste0("awk 'BEGIN{FS=\"\t\"}{if (NF == 8)print;}' ",
  #               tmp_vcf_file[['3']],
  #               " | awk '!seen[$2,$4,$5]++' >> ",
  #               target_vcf))
  system(paste0('bgzip -f -c ',
                target_vcf,' > ',
                target_vcf,'.gz'))
  system(paste0('tabix -f -p vcf ',target_vcf,'.gz'))
  system(paste0('rm -f ', target_vcf))
  system(paste0('rm -f ',tmp_vcf_file[['1']],'*'))
  system(paste0('rm -f ',tmp_vcf_file[['2']],'*'))
  system(paste0('rm -f ',tmp_vcf_file[['3']],'*'))

  lgr::lgr$info(paste0("Crossmapped VCF file: ", target_vcf,".gz"))

}

#' Crossmap BED file to different assembly version
#'
#' @description
#' Function that "lifts" the BED records to a different
#' assembly version of the human genome, using CrossMap.py
#'
#' @param direction hg38Tohg19 or hg19Tohg38
#' @param crossmap_cmd_path Full path to CrossMap.py file
#' @param chain_file_directory Directory containing chain files
#' @param wdir working directory
#' @param debug logical indicating if debug mode is on or off
#' @param fsep file name separator (default '_')
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

  source_records <- utils::read.table(file = source_bed, sep = "\t", quote = "")

  if (direction == "hg19Tohg38") {
    source_records$V1 <- stringr::str_replace_all(
      paste0('chr', source_records$V1), 'chrchr', 'chr')
    utils::write.table(source_records, file = tmp_bed_file[['3']],
                quote = F, row.names = F, col.names = F)
  }else{
    source_records$V1 <- stringr::str_replace_all(
      source_records$V1, 'chr', '')
    utils::write.table(source_records, file = tmp_bed_file[['3']],
                quote = F, row.names = F, col.names = F)
  }

  system(paste0(crossmap_cmd_path,' region -r 0.8 ',
                chainFile,' ', tmp_bed_file[['3']],
                ' ', tmp_bed_file[['1']]))

  #tmp_bed_file[['2']]
  system(paste0('egrep -v \'^##(lift|new)\' ',
                tmp_bed_file[['1']],
                ' | sed \'s/^chr//\' | ',
                'egrep -v \'Un|_alt|random\' > ',
                tmp_bed_file[['2']]))
  system(paste0('egrep \'^#\' ',tmp_bed_file[['2']],' > ',target_bed))
  system(paste0('cat ',tmp_bed_file[['2']],
                ' | egrep -v \"^#\" | egrep -v \"^[XYM]\" ',
                '| sort -k1,1n -k2,2n -k3,3n >> ',
                target_bed))
  system(paste0('cat ',tmp_bed_file[['2']],
                ' | egrep -v \"^#\" | egrep \"^[XYM]\" ',
                '| sort -k1,1 -k2,2n -k3,3n >> ',
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
