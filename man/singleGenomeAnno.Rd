% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getNCBI.R
\name{singleGenomeAnno}
\alias{singleGenomeAnno}
\title{Get single NCBI whole genomic gene annotation}
\usage{
singleGenomeAnno(genomeID, type = "gene", n = 1)
}
\arguments{
\item{genomeID}{Single NCBI genome ID.}

\item{type}{"gene" or "CDS". The KEGG database use "CDS" as the protein gene count.}

\item{n}{The number of CPUs or processors, and the default value is 1.}
}
\value{
A list of annotation.
}
\description{
NCBI Database API - Get single NCBI whole genomic gene annotation
}
\details{
Get whole gene annotation form single NCBI genome ID. The locus tag is used as names for each gene. If one of the gene feature value is missed, a "" (with length of 1) will return. If the genome has no gene featurs, "NULL" will be returned.
This function now supports two feature types, "gene" or "CDS" (coding sequence). Other features such as RNAs ("ncRNA", "rRNA", "tRNA", "tmRNA"), "misc_feature", "rep_origin", "repeat_region" are not supported yet. It is shown in E. coli K-12 MG1655 "genes" features not only includes all "CDS" and RNAs, but some sRNA ("b4714"). "misc_feature" are mainly about cryptic prophage genes, and "repeat_region" are repetitive extragentic palindromic (REP) elements.
}
\examples{
## no gene features
nofeature <- singleGenomeAnno('BA000048')

\dontrun{
aeuGenome <- singleGenomeAnno('CP007715', n = 4)

## missed value is replaced by ''
pseudo <- singleGenomeAnno('AE001826')[54]}
}
\author{
Yulong Niu \email{niuylscu@gmail.com}
}
