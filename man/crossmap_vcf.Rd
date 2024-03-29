% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossmap.R
\name{crossmap_vcf}
\alias{crossmap_vcf}
\title{Crossmap VCF file to different assembly version}
\usage{
crossmap_vcf(
  target_genome_file = "/Users/sigven/research/DB/hg19/ucsc.hg19.fa",
  direction = "hg38Tohg19",
  crossmap_cmd_path = "/Users/sigven/miniconda3/envs/py36/bin/CrossMap.py",
  chain_file_directory = "/Users/sigven/research/DB/chainFiles",
  fsep = "_",
  source_vcf = NULL,
  target_vcf = NULL,
  verbose = T
)
}
\arguments{
\item{target_genome_file}{fasta file of target genome (e.g. ucsc.hg19.fa)}

\item{direction}{hg38Tohg19 or hg19Tohg38}

\item{crossmap_cmd_path}{Full path to CrossMap.py}

\item{chain_file_directory}{Directory containing chain files}

\item{fsep}{separator to use in filenames (default '_')}

\item{source_vcf}{VCF file to be lifted}

\item{target_vcf}{VCF file with lifted coordinates}

\item{verbose}{logical indicating if function should log extensively}
}
\description{
Function that "lifts" the VCF records to a different
assembly version of the human genome, using CrossMap.py
}
